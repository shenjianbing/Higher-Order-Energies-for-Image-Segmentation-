#ifndef _FTRDIST_H__
#define _FTRDIST_H__
#include <FTRBase.h>
#include "boundutil.h"
#include "PPBCDist.h"

class FTRDist: public FTRBase{
public:
	FTRDist(double w_smooth_, double w_dist_, const Image &image_, 
		const vector<double> targetdist_, METRIC metric_)
		:image(image_),targetdist(targetdist_),metric(metric_)
	{
		w_smooth = w_smooth_;
		w_dist = w_dist_;
	}
	virtual void computeenergy(const Table2D<Label> & labeling, 
		double * act_energy_p, double * app_energy_p);
	// collect statistics about current labeling and update approximate model
	virtual void updatemodel();
	virtual GraphType * ftrbasegraph(double lambda_);
private:
	double w_smooth;
	double w_dist;
	Image image;

	vector<double> targetdist;
	METRIC metric;

	// statistics about current labeling
	vector<int> current_obj_hist;
	int current_obj_size;
	vector<double> app_a; // ax + b // weight not applied
	double app_b;
};

void FTRDist::computeenergy(const Table2D<Label> & labeling, 
		double * act_energy_p, double * app_energy_p)
{
	vector<int> obj_hist(image.colorbinnum,0);
	for(int i=0;i<img_w;i++)
	{
		for(int j=0;j<img_h;j++)
			if(labeling[i][j]==OBJ)
				obj_hist[image.colorlabel[i][j]]++;
	}
	double smooth_cost = getsmoothnesscost(image,labeling);
	outv(smooth_cost);
	double dist_cost;
	if(metric==KLDIVERGENCE)
		dist_cost = KLdistance(obj_hist,countintable(labeling,OBJ),targetdist);
	else if(metric == BHATTACHARYYA)
		dist_cost = BHAdistance(obj_hist,countintable(labeling,OBJ),targetdist);
	double act_energy = w_smooth*smooth_cost+w_dist*dist_cost;
	*act_energy_p = act_energy;
	if(app_energy_p)
	{
		double app_dist_cost = 0;
		for(int i=0;i<image.colorbinnum;i++)
		{
			app_dist_cost += app_a[i]*obj_hist[i];
		}
		app_dist_cost += app_b;
		double app_energy = w_smooth*smooth_cost+w_dist*app_dist_cost;
		*app_energy_p = app_energy;
	}
}
void FTRDist::updatemodel()
{
	// distance transform
	dt = getDistanceTransform(current_labeling);
	// update current foreground histogram
	int current_obj_size = countintable(current_labeling,OBJ);
	current_obj_hist = vector<int>(image.colorbinnum,0);
	for(int i=0;i<img_w;i++)
	{
		for(int j=0;j<img_h;j++)
			if(current_labeling[i][j]==OBJ)
				current_obj_hist[image.colorlabel[i][j]]++;
	}
	double common_a = 0;
	int s = current_obj_size;
	for(int i=0;i<image.colorbinnum;i++)
	{
		double t_i = targetdist[i];
		int s_i = current_obj_hist[i];
		if(metric==KLDIVERGENCE)
			common_a += t_i*s_i/((double)s*s)/((double)s_i/s+epsilon);
		else if(metric == BHATTACHARYYA)
			common_a += 0.5*sqrt(t_i*s_i/s/s/s);
	}
	// update a and b ( ax + b)
	app_a = vector<double>(image.colorbinnum,0);
	for(int i=0;i<image.colorbinnum;i++)
	{
		double t_i = targetdist[i];
		int s_i = current_obj_hist[i];
		// update a(i)
		app_a[i] = common_a;
		if(t_i==0)
			continue;
		if(metric==KLDIVERGENCE)
			app_a[i] += -t_i/((double)s_i/s+epsilon)/s;
		else if(metric == BHATTACHARYYA)
			app_a[i] += -0.5*sqrt(t_i/(s*max((double)1.0,(double)s_i)));
	}
	
	// update b
	double tmp = w_smooth * getsmoothnesscost(image,current_labeling);
	for(int i=0;i<image.colorbinnum;i++)
	{
		tmp += w_dist*app_a[i]*current_obj_hist[i];
	}
	app_b = (current_e - tmp)/w_dist;

}
GraphType * FTRDist::ftrbasegraph(double lambda_)
{
	GraphType * g = new GraphType(img_w*img_h,10*img_w*img_h);
	g->add_node(img_w*img_h);
	addsmoothnessterm(g, image, w_smooth, Table2D<bool>(img_w,img_h,true));
	for(int i=0;i<img_w;i++)
	{
		for(int j=0;j<img_h;j++)
		{
			int n = i+j*img_w;
			if(ROI[i][j])
				g->add_tweights(n,dt[i][j]*lambda_,app_a[image.colorlabel[i][j]]*w_dist);
			else // background
				g->add_tweights(n,0,INFTY);
		}
	}
	return g;
}

#endif