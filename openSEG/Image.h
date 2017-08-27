#pragma once
#include <ezi/Image2D.h>
#include <ezi/Table2D.h>
#include <SparseMatrix.h>
#include <assert.h>


// get average sigma_square for smoothness term
double getsmoothnessterm(const Table2D<RGB> &img,vector<PointPair> & pointpairs,int connecttype);
void rgb2indeximg(Table2D<int> & indeximg,const Table2D<RGB> & img,double channelstep);
int getcompactlabel(Table2D<int> & colorlabel,double channelstep,vector<int> & compacthist);
// edge-constrast sensitive smoothness penalty
inline double fn(const double dI, double lambda,double sigma_square);

class Image{
public:
	Image();
	Image(Table2D<RGB> img_, double channelstep, const char * imgname_ = "", int connecttype_ = 16);
	Image(const char * imgpath, const char * imgname_, double channelstep, int connecttype_);
	void computesmoothnesscost();
	void print();                        //打印相关信息

	Table2D<RGB> img;
	const char * imgname;
	int img_w;
	int img_h;
	int img_size;
	double sigma_square;
	vector<PointPair> pointpairs;
	vector<double> smoothnesscosts;
	int colorbinnum;
	Table2D<int> colorlabel;
	vector<int> compacthist; // color bin histogram (whole image) 这个直方图是去掉了为0的bin的
	int connecttype; // can be 4 or 8 or 16

	vector<Trituple<double>> pair_arcs;
};

Image::Image()
{

}

void Image::computesmoothnesscost()
{
	// number of neighboring pairs of pixels
	int numNeighbor = pointpairs.size();
	smoothnesscosts = vector<double>(numNeighbor);
	for(int i=0;i<numNeighbor;i++)
	{
		PointPair pp = pointpairs[i];
		smoothnesscosts[i] = fn(dI(img[pp.p1],img[pp.p2]),1.0,sigma_square)/(pp.p1-pp.p2).norm();
	} 
	//dI(a,b): a与b的颜色rgb差异的平方再求和
	//fn: exp(-...)  
	//(pp.p1-pp.p2).norm():  p1与p2的欧氏距离

	pair_arcs = vector<Trituple<double>>();
	for(int i=0;i<numNeighbor;i++)
	{
		Point p1 = pointpairs[i].p1;
		Point p2 = pointpairs[i].p2;
		pair_arcs.push_back(Trituple<double>(p1.x+p1.y*img_w,p2.x+p2.y*img_w,smoothnesscosts[i]));
		pair_arcs.push_back(Trituple<double>(p2.x+p2.y*img_w,p1.x+p1.y*img_w,smoothnesscosts[i]));
	}

}

//当申请了一个Image对象的时候，就会自动得到很多东西！！！
//图像，图像的名字，长宽，像素数量，几邻接，二阶项中的sigma^2
//pointpair向量（PointPair），smoothnesscosts向量（double），pair_arcs向量（三元组）。 它们长度都是邻居对的数量
//直方图compacthist(没有值为0的bin)，直方图的长度colorbinnum
//大小和图像一样的矩阵colorlabel
Image::Image(Table2D<RGB> img_, double channelstep, const char * imgname_, int connecttype_)
{
	assert((connecttype_==4)||(connecttype_==8)||(connecttype_==16),"wrong connect type!");
	imgname = imgname_;
	img = Table2D<RGB>(img_);
	img_w = img.getWidth();
	img_h = img.getHeight();
	img_size = img_w * img_h;

	connecttype = connecttype_;

	//这个函数不光是计算了sigma^2, 还得到了pointpairs向量
	sigma_square = getsmoothnessterm(img,pointpairs,connecttype);
	
	colorlabel= Table2D<int>(img_w,img_h);

	//这是很关键的一个函数，求解的是colorlabel: 每个像素掉到index为多少的直方图bin里头
	rgb2indeximg(colorlabel,img,channelstep);   //channelstep是直方图的每个bin的宽度

	//将原直方图中为0的bin去掉，新的直方图的bin的数量就是colorbinnum
	//新的直方图叫做compacthist
	//colorlabel也会被修改,存储的是每个像素属于新的直方图中的bin的序号
	colorbinnum = getcompactlabel(colorlabel,channelstep,compacthist);

	//更新了两个向量：长度都是邻居对的数量
	//smoothnesscosts：存的是两个邻居间的二阶项的能量值
	//pair_arcs：存三元组，分别是第一个，第二个 邻居 在图像中的index，和它们对应的smoothnesscosts
	computesmoothnesscost();
}


Image::Image(const char * imgpath, const char * imgname_, double channelstep, int connecttype_)
{
	Table2D<RGB> img_ = loadImage<RGB>(imgpath);
	new (this)Image(img_, channelstep, imgname_, connecttype_);
}
void Image::print()
{
	cout<<"sigma_square: "<<sigma_square<<endl;   //sigma^2,就是在二阶项的地方用的
	cout<<"image width: "<<img_w<<endl;
	cout<<"image height: "<<img_h<<endl;
	cout<<"# of pairs of neighborhoods"<<pointpairs.size()<<endl;    //邻居的数量：二阶项的数量
	cout<<"# of color bins"<<colorbinnum<<endl;                      //rgb三通道的直方图，在一阶项（高阶项）中，以这么多个bin的形式呈现
	cout<<"size of compact hist: "<<compacthist.size()<<endl;
}

//这是一个很关键的函数，因为他拿到了很多信息：
//得到二阶项需要的：sigma^2  和  所有的邻居关系对：pointpairs向量
double getsmoothnessterm(const Table2D<RGB> &img,vector<PointPair> & pointpairs, int connecttype)
{
	int node_id = 0;
	int img_w = img.getWidth();
	int img_h = img.getHeight();
	double sigma_sum = 0;
	double sigma_square_count = 0;
	Point kernelshifts [] = {Point(1,0),Point(0,1),Point(1,1),Point(1,-1),
		Point(2,-1),Point(2,1),Point(1,2),Point(-1,2),};
	for (int y=0; y<img_h; y++) // adding edges (n-links)
	{
		for (int x=0; x<img_w; x++) 
		{ 
			Point p(x,y);
			for(int i=0;i<connecttype/2;i++)
			{
				Point q = p + kernelshifts[i];
				if(img.pointIn(q))
				{
					sigma_sum += dI(img[p],img[q]);
					pointpairs.push_back(PointPair(p,q));
					sigma_square_count ++;
				}
			}
		}
	}
	return sigma_sum/sigma_square_count;
}

//这个函数是为了得到indeximg
//indeximg是一个图像大小的矩阵，每个的值表示的是这个位置的像素掉到的是3通道直方图第几个bin
/*eg:rgb值=（50，10，30） channelstep=16,表示的是bin的宽度，所以： 0~15 16~31 32~47 48~63 ...
那么这个像素掉到rgb分别为 (r_idx,g_idx,b_idx) = (3,0,2) 的bin里头。
在三通道的直方图中，一共有 （256/channelstep） ^3 = 4096个bin，序号分别是 0到4095
那么这个像素就掉到 3 + 0*16 +2*16*16 这个bin中啦
*/
void rgb2indeximg(Table2D<int> & indeximg,const Table2D<RGB> & img,double channelstep)
{
	RGB rgb_v;
	int r_idx =0, g_idx = 0, b_idx = 0, idx =0;
	int channelbin = (int)ceil(256.0/channelstep);  //每个通道的 直方图的 bin数量
	for(unsigned int j=0;j<img.getHeight();j++)
	{
		for(unsigned int i=0;i<img.getWidth();i++)
		{
			rgb_v = img[i][j];
			r_idx = (int)(rgb_v.r/channelstep);
			g_idx = (int)(rgb_v.g/channelstep);
			b_idx = (int)(rgb_v.b/channelstep);
			idx = r_idx + g_idx*channelbin+b_idx*channelbin*channelbin;
			indeximg[i][j] = idx;
		}
	}
}

//将原图的直方图的空的bin去掉，得到了新的直方图compacthist
//函数返回：新的直方图的bin的数量
//colorlabel会在函数中被修改:存储的是在compacthist中的序号
int getcompactlabel(Table2D<int> & colorlabel,double channelstep,vector<int> & compacthist)
{
	int channelbin = (int)ceil(256/channelstep);
	vector<int> colorhist(channelbin*channelbin*channelbin,0);
	for(unsigned int j=0;j<colorlabel.getHeight();j++)
	{
		for(unsigned int i=0;i<colorlabel.getWidth();i++)
		{
			colorhist[colorlabel[i][j]] = colorhist[colorlabel[i][j]]+1;
		}
	}  //得到原图的直方图
	
	vector<int> correspondence(colorhist.size(),-1);
	int compactcount = 0;
	for(unsigned int i=0;i<colorhist.size();i++)
	{
		if(colorhist[i]!=0)
		{
			compacthist.push_back(colorhist[i]);
			correspondence[i] = compactcount;
			compactcount++;
		}
	} //将原来直方图colorhist中为0的bin去掉，得到compacthist； 
	  //新的直方图compacthist的每一个bin，对应的在原直方图中的bin的序号，存在向量conrespondence里头

	for(unsigned int j=0;j<colorlabel.getHeight();j++)
	{
		for(unsigned int i=0;i<colorlabel.getWidth();i++)
		{
			colorlabel[i][j] = correspondence[colorlabel[i][j]];
		}
	}
	return compacthist.size();
}

inline double fn(const double dI, double lambda,double sigma_square) 
{
	return lambda*(exp(-dI/2/sigma_square));
}