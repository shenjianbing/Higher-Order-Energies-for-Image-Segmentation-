#pragma once
#include "Image.h"
#include <ezi/Cstr.h>
#include <entropy.h>

class BCD{
public:
	BCD(Image image_, double lambda_, double w_bits_ = 1, int max_iteration_ = 200, bool saveoption_ = false);
	void initlabeling(const Table2D<Label> & initiallabeling_);
	void updatemodel();    //更新前背景的直方图
	bool updatelabeling( int iter_count ); //建立graph，优化得到结果，成功/失败 返回 true/false
	bool optimize(int max_iteration_ = 0);
	double trivialsolutionenergy();

	int max_iteration;
	Image image;
	int img_w;
	int img_h;

	int OBJsize;
	int BKGsize;
	vector<int> OBJhist;  //存储前景直方图 长度为bin的数量
	vector<int> BKGhist;  //存储背景直方图

	Table2D<Label> hardconstraints;  //标记哪些像素一定为前景或背景

	Table2D<bool> ROI;  //region of interest 图像中要处理的区域，标记为1，不需要处理则标记为0
	double init_e;                   // energy of initial labeling
	double current_e;                // 当前labeling的能量
	Table2D<Label> current_labeling;

	// energy parameter
	double w_bits; //一阶项的weight
	double lambda; // weight of smoothness term. 二阶项的weight
	bool saveoption;
};

BCD::BCD(Image image_, double lambda_, double w_bits_, int max_iteration_, bool saveoption_)
	:image(image_),lambda(lambda_),max_iteration(max_iteration_),saveoption(saveoption_)
{
	img_w = image_.img_w;
	img_h = image_.img_h;
	w_bits = w_bits_;
	ROI = Table2D<bool>(img_w,img_h,true);
	hardconstraints = Table2D<Label>(img_w,img_h,NONE);
}

void BCD::initlabeling(const Table2D<Label> & initiallabeling_)
{
	current_labeling = initiallabeling_;
	ROI.reset(true);
	for(int y=0; y<img_h; y++)
	{
		for(int x=0; x<img_w; x++) 
		{ 
			if(current_labeling[x][y]==NONE)
				ROI[x][y] = false;
		}
	}
	init_e = getgrabcutenergy(image, w_bits, lambda, current_labeling);
	//cout<<"init energy: "<<init_e<<endl;
	current_e = init_e;
	updatemodel();
}
void BCD::updatemodel()
{
	OBJhist = vector<int>(image.colorbinnum,0);
	BKGhist = vector<int>(image.colorbinnum,0);
	OBJsize = 0;
	BKGsize = 0;
	for(int x =0;x<img_w;x++)
	{
		for(int y=0;y<img_h;y++)
		{
			if(current_labeling[x][y]==OBJ)
			{
				OBJsize++;
				OBJhist[image.colorlabel[x][y]]++;
			}
			else if(current_labeling[x][y]==BKG)
			{
				BKGsize++;
				BKGhist[image.colorlabel[x][y]]++;
			}
		}
	}
	//outv(OBJsize);
	//outv(BKGsize);
}

bool BCD::updatelabeling( int iter_count )
{
	//算法中，控制初始步长: lambda_alg
	float lambda_alg = 0.01 ;

	if((OBJsize==0)||(BKGsize==0))
		return false;
	GraphType * g = new GraphType(img_w*img_h, 4*img_w*img_h); 
	g->add_node(img_w*img_h);    // adding nodes


	// add smoothness term
	addsmoothnessterm(g, image, lambda, ROI);


	// add 高阶 term
	double ans , A , B , C , D , ti , In , FH , BH , FS , BS ;
	double step = 1 ;

	In = log(2.71828)/log(2.0) ;  //In = 1/ ln2 ;
	FS = (double)OBJsize ;
    BS = (double)BKGsize ;

	for(int j=0;j<img_h;j++)
	{
		for(int i=0;i<img_w;i++)
		{
			int n = i+j*img_w;
			if(hardconstraints[i][j]==BKG) // hard constraints to background
				g->add_tweights(n,0,INFTY);
			else if(hardconstraints[i][j]==OBJ) // hard constraints to foreground
				g->add_tweights(n,INFTY,0);
			else if(ROI[i][j])
			{
				int idx = image.colorlabel[i][j];

				if (current_labeling[i][j] == OBJ )
					ti=1.0 ;
				if (current_labeling[i][j] == BKG )
					ti=0.0 ;

				FH = (double)OBJhist[idx] ;
				BH = (double)BKGhist[idx] ;

/*				
				A = (1.0-BH)/(BH-ti) * In - log(BH-ti)/log(2.0) ;
				C = FH/(FH-1.0+ti) * In + log(FH-1.0+ti)/log(2.0) ;
				B = (1.0-BS)/(BS-ti) * In - log(BS-ti)/log(2.0) ;
				D = FS/(FS-1.0+ti) * In + log(FS-1.0+ti)/log(2.0) ;
*/				

				A = 0 - In - log(BH)/log(2.0) ;
				C = In + log(FH)/log(2.0) ;
				B = 0 - In - log(BS)/log(2.0) ;
				D = In + log(FS)/log(2.0) ;

				ans = B + D - A - C ;

				//cout<<"  "<<ans;

				g->add_tweights(n , 0 , ans*w_bits) ;//histogram natural log

				//算法中，控制步长的参数: lambda_alg*(inter_count+1)
				g->add_tweights(n, 0 , lambda_alg*(iter_count+1)*(1-2*ti) ) ;//histogram natural log

			}
		}
	}

	//////以上，把能量函数对应的图模型构建好了///////////

	//用maxflow去就得能量函数的解
	double newflow = g -> maxflow();

	// set next labeling， 导出优化的结果，做为下一次优化的constraint
	Table2D<Label> next_labeling(img_w,img_h,NONE);

	for (int y=0; y<img_h; y++) 
	{
		for (int x=0; x<img_w; x++) 
		{ 
			if(ROI[x][y])
			{
				int n = x+y*img_w;
				if(g->what_segment(n) == GraphType::SOURCE)
				{
					next_labeling[x][y]=OBJ;
				}
				else if(g->what_segment(n) == GraphType::SINK)
				{
					next_labeling[x][y]=BKG;
				}
			}
		}
	}
	delete g;
	g = NULL;

	double next_e = getgrabcutenergy(image, w_bits, lambda, next_labeling);
	                //得到了label，求这个结果的能量


	outv(current_e);
	//outv(next_e);


	//如果要增加迭代次数，把返回值就设成true就可以了
	if(next_e<current_e)   //-0.1
	{
		current_e = next_e;
		current_labeling = next_labeling;
		return true;
	}
	else
		return false;
}


bool BCD::optimize(int max_iteration_)
{
	int itr_count = 0;
	bool returnv = false;

	if(max_iteration_ !=0)
		max_iteration = max_iteration_;

	while(updatelabeling(itr_count)&&(itr_count<max_iteration))
	{
		itr_count++;
		//outv(current_e);
		//outv(itr_count);

		updatemodel();
		returnv = true;
		// assert(OBJsize!=0&&BKGsize!=0,"trivial solution for BCD!");
	}
	return returnv;
}

double BCD::trivialsolutionenergy()
{
	Table2D<Label> labeling(img_w,img_h,NONE);
	for (int y=0; y<img_h; y++) 
	{
		for (int x=0; x<img_w; x++) 
		{ 
			if(ROI[x][y])
				labeling[x][y] = OBJ;
		}
	}
	return getgrabcutenergy(image, w_bits, lambda, labeling);
}