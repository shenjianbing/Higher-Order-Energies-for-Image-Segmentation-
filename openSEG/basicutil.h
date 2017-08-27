// Basic utilities
#ifndef _BASICUTIL_H__
#define _BASICUTIL_H__
#include <maxflow/graph.h>
#include <ezi/Basics2D.h>
#include <ezi/cs1037lib-window.h> // for basic keyboard/mouse interface operations 
#include <ezi/cs1037lib-button.h> // for basic buttons, check-boxes, etc...
#include <ezi/Table2D.h>
#include <ezi/Math2D.h>
#include <ezi/Image2D.h>
#include <ezi/Cstr.h>

#include <Image.h>
#include <time.h>
#include <string>

#define outv(A) cout << #A << ": " << (double)(A) << endl; // output a Variable
#define outs(S) cout<<S<<endl; // output a String

#define INFTY 1e+20
#define EPS 1e-20 


enum Label {NONE=0, OBJ=1, BKG=2, UNKNOWN=3}; 
// UNKNOWN is introduced for monotonic parametric maxflow
// In this case, NONE means out of region of interest.
typedef Graph<double,double,double> GraphType;

// count certain element in table
template<typename T>
int countintable(const Table2D<T> & table, T t);
template<typename T>
void replaceintable(Table2D<T> & table, T oldt, T newt);
// count certain element in ROI in table
template<typename T>
int countintableROI(const Table2D<T> & table, T t, const Table2D<bool> & ROI);
// labeling.l ---> OBJ
// any other label ---> BKG
Table2D<Label> replacelabeling(const Table2D<Label> labeling, Label l);
// get labeling according to OBJ color
Table2D<Label> getinitlabeling(const Table2D<int> & initimg, int OBJcolor);
Table2D<Label> getinitlabelingFB(const Table2D<RGB> & initimg, RGB OBJcolor, RGB BKGcolor);
// save binary labeling as B&W image
void savebinarylabelingBW(const Table2D<Label> & labeling,string outname);
// save binary labeling
void savebinarylabeling(const Table2D<RGB> & img, const Table2D<Label> & labeling,string outname,bool BW = false);
// save multi labeling
template<typename T>
void savemultilabeling(Table2D<T> & labeling,const char * savedname,RGB * colors=NULL,Table2D<RGB> rgbimg = Table2D<RGB>());
// get labeling from maxflow instances
bool getlabeling(GraphType * g, Table2D<Label> & m_labeling);
// add smoothness term to the graph
// lambda is the weight of the smoothness term
// ROI is the region of interest
void addsmoothnessterm(GraphType * g, const Image & image, double lambda, 
	const Table2D<bool> & ROI, bool bordersmooth = false);
// error rate
double geterrorrate(Table2D<Label> & m_labeling,Table2D<int> & gtimg, int boxsize, int gtOBJcolor=0);
// smoothness cost - weight not applied
double getsmoothnesscost(const Image & image, const Table2D<Label> & m_labeling, bool bordersmoothness = false);
// distance transform based for L2 distance
// positive for object and negative for background
Table2D<double> getDistanceTransform(Table2D<Label> & labeling);

vector<int> getrandomvector(int n);
vector<Point> getrandomvector2dim(int n);
Table2D<Label> complementlabel(Table2D<Label> table);

Table2D<Label> complementlabel(Table2D<Label> table)
{
	int img_w = table.getWidth();
	int img_h = table.getHeight();
	Table2D<Label> complement(img_w,img_h,NONE);
	for(int i=0;i<img_w;i++)
	{
		for(int j=0;j<img_h;j++)
		{
			if(table[i][j]==OBJ)
				complement[i][j]=BKG;
			else if(table[i][j]==BKG)
				complement[i][j]=OBJ;
		}
	}
	return complement;
}
void savetableasgrayimage(Table2D<double> table, const char * imgname);

void savetableasgrayimage(Table2D<double> table, const char * imgname)
{
	int img_w = table.getWidth();
	int img_h = table.getHeight();

	for(int i=0;i<img_w;i++)
	{
		for(int j=0;j<img_h;j++)
		{
			if(table[i][j]>1000)
				table[i][j]=0;
		}
	}
	double tablemax = table.getMax();
	Table2D<int> tmp(img_w,img_h);
	for(int i=0;i<img_w;i++)
	{
		for(int j=0;j<img_h;j++)
		{
			int grayvalue = (int)(table[i][j]*10);
			if(grayvalue>255) grayvalue=255;
			tmp[i][j] = grayvalue;
		}
	}
	saveImage(tmp, imgname);
	outv(tmp.getMax());
	cout<<"saved into "<<imgname<<endl;
}

template<typename T>
int countintable(const Table2D<T> & table, T t)
{
	int table_w = table.getWidth();
	int table_h = table.getHeight();
	int tsize = 0;
	for(int y=0; y<table_h; y++) 
	{
		for(int x=0; x<table_w; x++) 
		{ 
			if(table[x][y]==t) // certain element t
				tsize++;
		}
	}
	return tsize;
}

template<typename T>
void replaceintable(Table2D<T> & table, T oldt, T newt)
{
	int table_w = table.getWidth();
	int table_h = table.getHeight();
	for(int y=0; y<table_h; y++) 
	{
		for(int x=0; x<table_w; x++) 
		{ 
			if(table[x][y]==oldt) // certain element t
				table[x][y]=newt;
		}
	}
}

template<typename T>
int countintableROI(const Table2D<T> & table, T t, const Table2D<bool> & ROI)
{
	int table_w = table.getWidth();
	int table_h = table.getHeight();
	int tsize = 0;
	for(int y=0; y<table_h; y++) 
	{
		for(int x=0; x<table_w; x++) 
		{ 
			if((table[x][y]==t)&&(ROI[x][y])) // certain element t
				tsize++;
		}
	}
	return tsize;
}

// labeling.l ---> OBJ
// any other label ---> BKG
Table2D<Label> replacelabeling(const Table2D<Label> labeling, Label l)
{
	Table2D<Label> newlabeling(labeling.getWidth(),labeling.getHeight(),BKG);
	for(unsigned int i=0;i<labeling.getWidth();i++)
	{
		for(unsigned int j=0;j<labeling.getHeight();j++)
		{
			if(labeling[i][j]==l)
				newlabeling[i][j]=OBJ;
		}
	}
	return newlabeling;
}

Table2D<Label> getinitlabeling(const Table2D<int> & initimg, int OBJcolor)
{
	int img_w = initimg.getWidth();
	int img_h = initimg.getHeight();
	Table2D<Label> initlabeling(img_w,img_h);
	for(int i=0;i<img_w;i++)
	{
		for(int j=0;j<img_h;j++)
		{
			if(initimg[i][j]==OBJcolor)
				initlabeling[i][j] = OBJ;
			else
				initlabeling[i][j] = BKG;
		}
	}
	return initlabeling;
}

Table2D<Label> getinitlabelingFB(const Table2D<RGB> & initimg, RGB OBJcolor, RGB BKGcolor)
{
	int img_w = initimg.getWidth();
	int img_h = initimg.getHeight();
	Table2D<Label> initlabeling(img_w,img_h);
	for(int i=0;i<img_w;i++)
	{
		for(int j=0;j<img_h;j++)
		{
			if(initimg[i][j]==OBJcolor)
				initlabeling[i][j] = OBJ;
			else if(initimg[i][j]==BKGcolor)
				initlabeling[i][j] = BKG;
			else
				initlabeling[i][j] = NONE;
		}
	}
	return initlabeling;
}

void savebinarylabelingBW(const Table2D<Label> & labeling,string outname)
{
	int img_w = labeling.getWidth();
	int img_h = labeling.getHeight();
	Table2D<RGB> tmp(img_w,img_h);
	for(int i=0;i<img_w;i++)
	{
		for(int j=0;j<img_h;j++)
		{
			if(labeling[i][j]==BKG)
				tmp[i][j] = white;
			else if(labeling[i][j]==OBJ)
				tmp[i][j] = black;
			else
				tmp[i][j] = gray;
		}
	}
	saveImage(tmp, to_Cstr(outname));
	cout<<"saved into: "<<outname<<endl;
}

void savebinarylabeling(const Table2D<RGB> & img, const Table2D<Label> & labeling,string outname,bool BW)
{
	int img_w = labeling.getWidth();
	int img_h = labeling.getHeight();
	Table2D<RGB> tmp(img_w,img_h);
	for(int i=0;i<img_w;i++)
	{
		for(int j=0;j<img_h;j++)
		{
			if(labeling[i][j]==OBJ)
			{
				if(BW) tmp[i][j] = black;
				else tmp[i][j] = img[i][j];
			}
			else
				tmp[i][j] = silver;
		}
	}
	saveImage(tmp, to_Cstr(outname));
	cout<<"saved into: "<<outname<<endl;
}

template<typename T>
void savemultilabeling(Table2D<T> & labeling,const char * outname,RGB * colors,Table2D<RGB> rgbimg)
{
	int img_w = labeling.getWidth();
	int img_h = labeling.getHeight();
	Table2D<RGB> tmp(img_w,img_h);
	if(colors!=NULL)
	for(int i=0;i<img_w;i++)
	{
		for(int j=0;j<img_h;j++)
		{
			tmp[i][j] = colors[labeling[i][j]];
		}
	}
	else
	{
		int numseg = labeling.getMax()+1;
		vector<int> r_vec(numseg,0);
		vector<int> g_vec(numseg,0);
		vector<int> b_vec(numseg,0);
		vector<int> labelhist(numseg,0);
		vector<RGB> averagecolors(numseg,0);
		for(int i=0;i<img_w;i++)
		{
			for(int j=0;j<img_h;j++)
			{
				T l = labeling[i][j];
				labelhist[l]++;
				r_vec[l]+=rgbimg[i][j].r;
				g_vec[l]+=rgbimg[i][j].g;
				b_vec[l]+=rgbimg[i][j].b;
			}
		}
		for(int l=0;l<numseg;l++)
			if(labelhist[l])
			averagecolors[l] = RGB(r_vec[l]/labelhist[l],g_vec[l]/labelhist[l],b_vec[l]/labelhist[l]);
		for(int i=0;i<img_w;i++)
		{
			for(int j=0;j<img_h;j++)
			{
				tmp[i][j] = averagecolors[labeling[i][j]];
			}
		}
	}
	saveImage(tmp, to_Cstr(outname));
	cout<<"saved into: "<<outname<<endl;
}

bool getlabeling(GraphType * g, Table2D<Label> & m_labeling)
{
	int img_w = m_labeling.getWidth();
	int img_h = m_labeling.getHeight();
	int n=0;
	m_labeling.reset(NONE);
	int sumobj =0, sumbkg = 0;
	for (int y=0; y<img_h; y++) 
	{
		for (int x=0; x<img_w; x++) 
		{ 
			n = x+y*img_w;
			if(g->what_segment(n) == GraphType::SOURCE)
			{
				m_labeling[x][y]=OBJ;
				sumobj++;
			}
			else if(g->what_segment(n) == GraphType::SINK)
			{
				m_labeling[x][y]=BKG;
				sumbkg++;
			}
		}
	}
	if((sumobj==0)||(sumbkg==0))
		return false;
	else
		return true;
}

// add smoothness term to the graph
// lambda is the weight of the smoothness term
// ROI is the region of interest
void addsmoothnessterm(GraphType * g, const Image & image, double lambda,
	const Table2D<bool> & ROI, bool bordersmooth)
{
	int img_w = image.img_w;
	int img_h = image.img_h;
	// number of neighboring pairs of pixels
	int numNeighbor = image.pointpairs.size();
	// n-link - smoothness term
	int node_id1 =0, node_id2 =0;
	
	for(int i=0;i<numNeighbor;i++)
	{
		PointPair pp = image.pointpairs[i];
		node_id1 = pp.p1.x+pp.p1.y*img_w;
		node_id2 = pp.p2.x+pp.p2.y*img_w;
		// if the two points are inside region of interest.
		if(ROI[pp.p1]&&ROI[pp.p2])
		{  
			double v = lambda*image.smoothnesscosts[i];
			g->add_edge(node_id1,node_id2,v,v);
		}
	}
	if(bordersmooth)
	{
		for(int i=0;i<image.img_w;i++)
		{
			for(int j=0;j<image.img_h;j++)
				if((i==0)||(i==img_w-1)||(j==0)||(j==img_h-1)) // border
					g->add_tweights(i+j*img_w,0,lambda);
		}
	}
}

double geterrorrate(Table2D<Label> & m_labeling,Table2D<int> & gtimg, int boxsize, int gtOBJcolor)
{
	double errorrate = 0 ;
	int errornum = 0;
	for(unsigned int j=1;j<gtimg.getHeight()-1;j++)
	{
		for(unsigned int i=1;i<gtimg.getWidth()-1;i++)
		{
			if((gtimg[i][j]==gtOBJcolor)&&(m_labeling[i][j]==OBJ))
				errornum++;
			else if((gtimg[i][j]==(255-gtOBJcolor))&&(m_labeling[i][j]==BKG))
				errornum++;
		}
	}
	errorrate = (double)errornum / boxsize;
	return errorrate;
}

double getsmoothnesscost(const Image & image, const Table2D<Label> & m_labeling, bool bordersmoothness)
{
	int img_w = image.img_w;
	int img_h = image.img_h;
	double smoothenergy = 0;
	// number of neighboring pairs of pixels
	int numNeighbor = image.pointpairs.size();
	for(int i=0;i<numNeighbor;i++)
	{
		PointPair pp = image.pointpairs[i];
		if((m_labeling[pp.p1]!=NONE)&&(m_labeling[pp.p2]!=NONE))
		{
			if(m_labeling[pp.p1]!=m_labeling[pp.p2])
				smoothenergy += image.smoothnesscosts[i];
		}
		if((m_labeling[pp.p1]==NONE)||(m_labeling[pp.p2]==NONE))
			exit(-1);
	}
	if(bordersmoothness)
	{
		for(int i=0;i<image.img_w;i++)
		{
			for(int j=0;j<image.img_h;j++)
				if((i==0)||(i==img_w-1)||(j==0)||(j==img_h-1)) // border
					if(m_labeling[i][j]==OBJ)
						smoothenergy = smoothenergy +1;
		}
	}
	return smoothenergy;
}

////////////////////////////////////////////////////////////////////////
// getDistanceTransform is a function that computes distance 
// transform. It implements forward-backward pass algorithm with
// windows that approximate Euclidean Metric.
Table2D<double> getDistanceTransform(Table2D<Label> & labeling)
{
	int img_w = labeling.getWidth();
	int img_h = labeling.getHeight();
	Table2D<double> dt(img_w,img_h,INFTY); // distance transform
	//find all the edges
	for(int i=0;i<img_w;i++)
	{
		for(int j=0;j<img_h;j++)
		{
			if((i==0)||(i==img_w-1)||(j==0)||(j==img_h-1)) // border
			{
				if(labeling[i][j]==OBJ)
					dt[i][j] = 0;
			}
			else // non-border
			{
				if((labeling[i][j]==OBJ)&&
				((labeling[i+1][j]==BKG)||(labeling[i-1][j]==BKG)||(labeling[i][j+1]==BKG)||(labeling[i][j-1]==BKG)))
					dt[i][j] = 0;
			}
		}
	}
	//forward pass
	for(int i=0;i<img_w;i++)
	{
		for(int j=0;j<img_h;j++)
		{
			Point p(i,j);
			Point q;
			q=p+Point(-1,-1);
			if(dt.pointIn(q))
				dt[p]=min(dt[p],(double)1.414+dt[q]);
			q=p+Point(0,-1);
			if(dt.pointIn(q))
				dt[p]=min(dt[p],(double)1.0+dt[q]);
			q=p+Point(1,-1);
			if(dt.pointIn(q))
				dt[p]=min(dt[p],(double)1.414+dt[q]);
			q=p+Point(-1,0);
			if(dt.pointIn(q))
				dt[p]=min(dt[p],(double)1.0+dt[q]);
		}
	}
	//backward pass
	for(int i=img_w-1;i>=0;i--)
	{
		for(int j=img_h-1;j>=0;j--)
		{
			Point p(i,j);
			Point q;
			q=p+Point(1,1);
			if(dt.pointIn(q))
				dt[p]=min(dt[p],(double)1.414+dt[q]);
			q=p+Point(0,1);
			if(dt.pointIn(q))
				dt[p]=min(dt[p],(double)1.0+dt[q]);
			q=p+Point(-1,1);
			if(dt.pointIn(q))
				dt[p]=min(dt[p],(double)1.414+dt[q]);
			q=p+Point(1,0);
			if(dt.pointIn(q))
				dt[p]=min(dt[p],(double)1.0+dt[q]);
		}
	}

	// second time forward and backward

	//forward pass
	for(int j=0;j<img_h;j++)
	{
		for(int i=0;i<img_w;i++)
		{
			Point p(i,j);
			Point q;
			q=p+Point(-1,-1);
			if(dt.pointIn(q))
				dt[p]=min(dt[p],(double)1.414+dt[q]);
			q=p+Point(0,-1);
			if(dt.pointIn(q))
				dt[p]=min(dt[p],(double)1.0+dt[q]);
			q=p+Point(1,-1);
			if(dt.pointIn(q))
				dt[p]=min(dt[p],(double)1.414+dt[q]);
			q=p+Point(-1,0);
			if(dt.pointIn(q))
				dt[p]=min(dt[p],(double)1.0+dt[q]);
		}
	}
	//backward pass
	for(int j=img_h-1;j>=0;j--)
	{
		for(int i=img_w-1;i>=0;i--)
		{
			Point p(i,j);
			Point q;
			q=p+Point(1,1);
			if(dt.pointIn(q))
				dt[p]=min(dt[p],(double)1.414+dt[q]);
			q=p+Point(0,1);
			if(dt.pointIn(q))
				dt[p]=min(dt[p],(double)1.0+dt[q]);
			q=p+Point(-1,1);
			if(dt.pointIn(q))
				dt[p]=min(dt[p],(double)1.414+dt[q]);
			q=p+Point(1,0);
			if(dt.pointIn(q))
				dt[p]=min(dt[p],(double)1.0+dt[q]);
		}
	}

	for(int i=0;i<img_w;i++)
	{
		for(int j=0;j<img_h;j++)
		{
			if(labeling[i][j]==OBJ)
				dt[i][j] = dt[i][j]+0.5;
			else if(labeling[i][j]==BKG)
				dt[i][j] = -(dt[i][j]-0.5);
		}
	}
	return dt;
}

vector<int> getrandomvector(int n)
{
	vector<int> v(n,0);
	for (int i=0;i<n;i++)
		v[i]=i;
	int p;
	int tmp;

	for (int i=n-1;i>0;i--)
	{
		p=rand()%i;
		tmp=v[p];
		v[p]=v[i];
		v[i]=tmp;
	}
	return v;
}

vector<Point> getrandomvector2dim(int n)
{
	vector<Point> v;
	for(int i=0;i<n-1;i++)
		for(int j=i+1;j<n;j++)
			v.push_back(Point(i,j));
	int size_v = v.size();

	for (int i=size_v-1;i>0;i--)
	{
		int p=rand()%i;
		Point tmp=v[p];
		v[p]=v[i];
		v[i]=tmp;
	}
	return v;
}

#endif