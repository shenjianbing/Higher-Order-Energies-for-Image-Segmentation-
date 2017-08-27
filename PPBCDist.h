#ifndef _PPBCDIST_H__
#define _PPBCDIST_H__
// class for parametric maxflow (matching target distribution)
#include "PPBCBase.h"
#include "boundutil.h"
#include "SparseMatrix.h"

// global paremeter
double epsilon = 1/1000000.0; // for KL divergence

double distmatchingenergy(double w_smooth, double w_dist, const Image & image, 
	const SEGSOLUTION & segsolution, const vector<double> & targetdist, METRIC metric);
double KLdistance(const vector<int> & objhist, int ssize, const vector<double> & targetdist);
double BHAdistance(const vector<int> & objhist, int ssize, const vector<double> & targetdist);
double distmatchingenergy(double w_smooth, double w_dist, const Image & image, 
	const Table2D<Label> & labeling, const vector<double> & targetdist, METRIC metric);

class PPBCDist: public PPBCBase{
public:
	PPBCDist(double w_smooth_, double w_dist_, const Image &image_, 
		const vector<double> targetdist_, METRIC metric_)
		:image(image_),targetdist(targetdist_),metric(metric_)
	{
		w_smooth = w_smooth_;
		w_dist = w_dist_;
	}
	virtual double computeenergy(const Table2D<Label> & labeling);
	virtual void updatemodel(); // collect statistics about current labeling
	virtual void * parabasegraph(double para_, UnknownRegion * unknownregion_p,
		double * flowoffset_p, vector<Point> * node_corr_p);
	virtual double getnewpara(ParaInterval paraInterval);
private:
	double w_smooth;
	double w_dist;
	Image image;

	vector<double> targetdist;
	METRIC metric;
	SEGSOLUTION initsolution;
};

double PPBCDist::computeenergy(const Table2D<Label> & labeling)
{
	return distmatchingenergy(w_smooth, w_dist, image, labeling, 
		targetdist, metric);
}
void PPBCDist::updatemodel()
{
	initsolution = SEGSOLUTION(initlabeling,image);
}
void * PPBCDist::parabasegraph(double para_, UnknownRegion * unknownregion_p,
	double * flowoffset_p, vector<Point> * node_corr_p)
{
	assert(initsolution.ssize!=0);
	int histlen = targetdist.size();
	vector<double> slopes(histlen,0);
	for(int i=0;i<histlen;i++)
	{
		double t_i = targetdist[i];
		int s_i = initsolution.objhist[i];
		if(s_i!=0)
		{
			if(metric == KLDIVERGENCE)
				slopes[i] = -w_dist*t_i*(log(epsilon)-log((double)s_i/initsolution.ssize+epsilon))/s_i;
			else if(metric == BHATTACHARYYA)
				slopes[i] = w_dist*sqrt(t_i*s_i/(double)initsolution.ssize)/s_i;
		}
	}

	Table2D<double> capsource(img_w,img_h,0);
	Table2D<double> capsink(img_w,img_h,0);
	for(int i=0;i<img_w;i++)
	{
		for(int j=0;j<img_h;j++)
		{
			int n = i + j*img_w;
			if(initsolution.labeling[i][j] == BKG)
			{
				capsource[i][j] = 0;
				capsink[i][j] = INFTY;
			}
			else
			{
				capsource[i][j] = slopes[image.colorlabel[i][j]]+para_;
				capsink[i][j] = 0;
			}
		}
	}

	if(unknownregion_p==NULL)
	{
		GraphType * g = new GraphType(img_w*img_h,10*img_w*img_h);
		g->add_node(img_w*img_h);
		addsmoothnessterm(g, image, w_smooth, ROI);
		for(int i=0;i<img_w;i++)
		{
			for(int j=0;j<img_h;j++)
			{
				int n = i + j*img_w;
				g->add_tweights(n,capsource[i][j],capsink[i][j]);
			}
		}
		return g;
	}
	// for monotonic mode
	Table2D<int> img_corr(img_w,img_h,-1);
	double flowoffset = 0;

	flowoffset = 0; 
	Table2D<bool> incompactgraph;
	GraphType * compact_g;
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

double PPBCDist::getnewpara(ParaInterval interval)
{
	// ballooning cost for lower parameter solution
	int ssize_low = interval.bplow.ssize;
	double ballooncost_low = (initsolution.ssize-ssize_low)*interval.bplow.para;
	double Entropy_low = interval.bplow.flow-ballooncost_low;
	// ballooning cost for upper parameter solution
	int ssize_up = interval.bpup.ssize;
	double ballooncost_up = (initsolution.ssize-ssize_up)*interval.bpup.para;
	double Entropy_up = interval.bpup.flow-ballooncost_up;
	double newpara = -(Entropy_low - Entropy_up) / (double)(ssize_up - ssize_low);
	return newpara;
}

double distmatchingenergy(double w_smooth, double w_dist, const Image & image, 
	const SEGSOLUTION & segsolution, const vector<double> & targetdist,METRIC metric)
{
	double dist = 0;
	if(metric == KLDIVERGENCE)
		dist = KLdistance(segsolution.objhist,segsolution.ssize,targetdist);
	else if(metric == BHATTACHARYYA)
		dist = BHAdistance(segsolution.objhist,segsolution.ssize,targetdist);
	double smoothcost = getsmoothnesscost(image,segsolution.labeling);
	return dist*w_dist + w_smooth * smoothcost;
}

double KLdistance(const vector<int> & objhist, int ssize, const vector<double> & targetdist)
{
	if(ssize==0)
		return INFTY;
	double dist = 0;
	int histlen = objhist.size();
	for(int i=0;i<histlen;i++)
		if(targetdist[i]!=0)
			dist += -targetdist[i]*log(objhist[i]/(double)ssize+epsilon);
	return dist;
}

double BHAdistance(const vector<int> & objhist, int ssize, const vector<double> & targetdist)
{
	assert(ssize!=0,"for BHA distance, foreground size should not be zero!");
	double dist = 0;
	int histlen = objhist.size();
	for(int i=0;i<histlen;i++)
			dist += -sqrt((double)targetdist[i]*objhist[i]/(double)ssize);
	return dist;
}

double distmatchingenergy(double w_smooth, double w_dist, const Image & image, 
	const Table2D<Label> & labeling, const vector<double> & targetdist,METRIC metric)
{
	double dist = 0;
	SEGSOLUTION segsolution(labeling,image);
	if(metric == KLDIVERGENCE)
		dist = KLdistance(segsolution.objhist,segsolution.ssize,targetdist);
	else if(metric == BHATTACHARYYA)
		dist = BHAdistance(segsolution.objhist,segsolution.ssize,targetdist);
	double smoothcost = getsmoothnesscost(image,segsolution.labeling);
	return dist*w_dist + w_smooth * smoothcost;
}

#endif