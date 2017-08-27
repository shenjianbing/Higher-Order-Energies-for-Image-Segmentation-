#ifndef _PPBCEntropy_H__
#define _PPBCEntropy_H__
// class for parametric maxflow (minimize BCD energy)
#include "PPBCBase.h"
#include <matlabIO.h>
#include <engine.h>
//extern Engine  *ep; // matlab engine

class PPBCEntropy: public PPBCBase{
public:
	PPBCEntropy(const Image & image_, double w_smooth_, double w_bits_)
		:image(image_),w_smooth(w_smooth_),w_bits(w_bits_){};

	virtual double computeenergy(const Table2D<Label> & labeling); //得到现在这个label下的能量值
	virtual void updatemodel(); // collect statistics about current labeling 更新前背景的直方图
	virtual void * parabasegraph(double para_, UnknownRegion * unknownregion_p,
		double * flowoffset_p, vector<Point> * node_corr_p);  //？？
	virtual double getnewpara(ParaInterval paraInterval);
	Table2D<Label> hardconstraints;
private:
	int OBJsize;
	int BKGsize;
	vector<int> OBJhist;
	vector<int> BKGhist;
	double w_smooth;
	double w_bits;
	Image image;
};

double PPBCEntropy::computeenergy(const Table2D<Label> & labeling)
{
	return getgrabcutenergy(image, w_bits, w_smooth, labeling);
}
void PPBCEntropy::updatemodel()
{
	objbkghist(image, OBJhist, BKGhist, OBJsize, BKGsize, initlabeling);
}
void * PPBCEntropy::parabasegraph(double para_, UnknownRegion * unknownregion_p,
	double * flowoffset_p, vector<Point> * node_corr_p)
{
	Table2D<double> capsource(img_w,img_h,0);
	Table2D<double> capsink(img_w,img_h,0);

	for(int j=0;j<img_h;j++)
	{
		for(int i=0;i<img_w;i++)
		{
			int n = i+j*img_w;
			if(hardconstraints.pointIn(i,j)&&hardconstraints[i][j]==BKG) // hard constraints to background
			{
				capsource[i][j]=0;
				capsink[i][j]=INFTY;
			}
			else if(hardconstraints.pointIn(i,j)&&hardconstraints[i][j]==OBJ) // hard constraints to foreground
			{
				capsource[i][j]=INFTY;
				capsink[i][j]=0;
			}
			else if(ROI[i][j])
			{
				int idx = image.colorlabel[i][j];
				if(BKGsize)
					capsource[i][j]  = (-log(max((double)BKGhist[idx],(double)1.0))+log((double)BKGsize))
						/log(2.0)*w_bits;
				else
					capsource[i][j] = -log(1.0/(double)image.img_size)/log(2.0)*w_bits; // uniform model
				if(OBJsize)
					capsink[i][j]  = (-log(max((double)OBJhist[idx],(double)1.0))+log((double)OBJsize))
						/log(2.0)*w_bits;
				else
					capsink[i][j] = -log(1.0/(double)image.img_size)/log(2.0)*w_bits; // uniform model
				capsink[i][j] += para_;
			}
		}
	}
	if(unknownregion_p==NULL){
	GraphType * g = new GraphType(img_w*img_h, 4*img_w*img_h); 
	g->add_node(img_w*img_h);    // adding nodes
	// add smoothness term
	addsmoothnessterm(g, image, w_smooth, ROI);
	// add unary term
	for(int j=0;j<img_h;j++)
	{
		for(int i=0;i<img_w;i++)
		{
			if(ROI[i][j])
			{
				int n = i+j*img_w;
				g->add_tweights(n,capsource[i][j],capsink[i][j]);//histogram natural log
			}
		}
	}
	return g;}

	// monotonic speedup
	// for monotonic mode
	double flowoffset = 0; 
	Table2D<bool> incompactgraph;
	GraphType * compact_g;
	Table2D<int> img_corr;
	vector<PointPair> compactpointpairs;
	vector<double> compactsmoothnesscosts;
	int compactsize = getcompactgraph(image,unknownregion_p->knownlabeling, incompactgraph,
		*node_corr_p,img_corr, compactpointpairs, compactsmoothnesscosts);
	compact_g = new GraphType(compactsize,10*compactsize);
	compact_g->add_node(compactsize);

	// number of neighboring pairs of pixels
	int compactnumNeighbor = compactpointpairs.size();
	// n-link - smoothness term
	for(int i=0;i<compactnumNeighbor;i++)
	{
		PointPair pp = compactpointpairs[i];
		if(incompactgraph[pp.p1]&&incompactgraph[pp.p2])
		{
			double v = w_smooth*compactsmoothnesscosts[i];
			compact_g->add_edge(img_corr[pp.p1],img_corr[pp.p2],v,v);
		}
	}

	for(int i=0;i<img_w;i++)
	{
		for(int j=0;j<img_h;j++)
		{
			int n = img_corr[i][j];
			if(incompactgraph[i][j])
			{
				if(unknownregion_p->knownlabeling[i][j]==OBJ)
					compact_g->add_tweights(n,INFTY,capsink[i][j]);
				else if(unknownregion_p->knownlabeling[i][j]==BKG)
					compact_g->add_tweights(n,capsource[i][j],INFTY);
				else if(unknownregion_p->knownlabeling[i][j]==UNKNOWN)
					compact_g->add_tweights(n,capsource[i][j],capsink[i][j]);
			}
			else
			{
				if(unknownregion_p->knownlabeling[i][j]==OBJ)
					flowoffset += capsink[i][j];
				else if(unknownregion_p->knownlabeling[i][j]==BKG)
					flowoffset += capsource[i][j];
			}
		}
	}
	if(flowoffset_p) *flowoffset_p = flowoffset;
	return compact_g;

}

double PPBCEntropy::getnewpara(ParaInterval interval)
{
	// ballooning cost for lower parameter solution
	int ssize_low = interval.bplow.ssize;
	double ballooncost_low = ssize_low*interval.bplow.para;
	double entropy_low = interval.bplow.flow-ballooncost_low;
	// ballooning cost for upper parameter solution
	int ssize_up = interval.bpup.ssize;
	double ballooncost_up = ssize_up*interval.bpup.para;
	double entropy_up = interval.bpup.flow-ballooncost_up;
	//out(ssize_up - ssize_low);
	double newpara = (entropy_low - entropy_up) / (double)(ssize_up - ssize_low);
	return newpara;
}

#endif