#ifndef _BOUNDUTIL_H__
#define _BOUNDUTIL_H__
#include "utilities.h"

enum METRIC { 
	KLDIVERGENCE=0, 
	BHATTACHARYYA=1};

vector<int> gettargethistogram(const Image & image, const Table2D<int> & gtimg);
vector<int> getobjhist(const Image & image, const Table2D<Label> & initlabeling);
int histogramdistance(const vector<int> & objhist, const vector<int> & targethist, int metric);
vector<double> gettargetdistribution(const Image & image, const Table2D<int> & gtimg);

class SEGSOLUTION{
public:
	bool operator<(const SEGSOLUTION& b) const   {return energy < b.energy;}
	SEGSOLUTION()
	{}
	SEGSOLUTION(double energy_)
	{
		energy = energy_;
	}
	SEGSOLUTION(const Table2D<Label> & labeling_,const Image & image)
	{
		labeling = labeling_;
		int img_w = labeling.getWidth();
		int img_h = labeling.getHeight();
		ssize = 0;
		energy = 0;
		objhist = vector<int>(image.colorbinnum,0);
		for(int i=0;i<img_w;i++)
		{
			for(int j=0;j<img_h;j++)
			{
				if(labeling[i][j]==OBJ)
				{
					ssize++;
					objhist[image.colorlabel[i][j]]++;
				}
			}
		}
	}
	void print()
	{
		cout<<"foreground size: "<<ssize<<endl;
		cout<<"energy "<<energy<<endl;
	}
	Table2D<Label> labeling;
	double energy;
	int ssize; // foreground size
	vector<int> objhist;
};

void addhistogramterm(GraphType * g, const vector<int> & targethist, int metric, 
	const Image & image, const SEGSOLUTION & segsolution, double w_hist);

vector<int> gettargethistogram(const Image & image, const Table2D<int> & gtimg)
{
	vector<int> targethist(image.colorbinnum,0);
	for(int i=0;i<image.img_w;i++)
	{
		for(int j=0;j<image.img_h;j++)
		{
			if(gtimg[i][j]==255) // white
				targethist[image.colorlabel[i][j]]++;
		}
	}
	return targethist;
}


vector<int> getobjhist(const Image & image, const Table2D<Label> & initlabeling)
{
	vector<int> objhist(image.colorbinnum,0);
	for(int i=0;i<image.img_w;i++)
	{
		for(int j=0;j<image.img_h;j++)
		{
			if(initlabeling[i][j]==OBJ) // white
				objhist[image.colorlabel[i][j]]++;
		}
	}
	return objhist;
}

int histogramdistance(const vector<int> & objhist, const vector<int> & targethist, int metric)
{
	int dist = 0;
	int histlen = objhist.size();
	for(int i=0;i<histlen;i++)
		dist += (int)pow(abs(objhist[i]-targethist[i]),(double)metric);
	return dist;
}

double histogrammatchingenergy(double w_smooth, double w_hist, const Image & image, 
	const SEGSOLUTION & segsolution, const vector<int> & targethist, int metric)
{
	return histogramdistance(segsolution.objhist,targethist,metric)*w_hist 
		+ w_smooth * getsmoothnesscost(image,segsolution.labeling);
}

void addhistogramterm(GraphType * g, const vector<int> & targethist, int metric,
	const Image & image, const SEGSOLUTION & segsolution, double w_hist)
{
	int histlen = targethist.size();
	int img_w = segsolution.labeling.getWidth();
	int img_h = segsolution.labeling.getHeight();
	vector<double> slopes(histlen,0);
	for(int i=0;i<histlen;i++)
	{
		int t_i = targethist[i];
		int s_i = segsolution.objhist[i];
		if(s_i!=0)
		{
			if(metric==2)
				slopes[i] = 2*t_i-s_i;
			else if(metric==1)
				slopes[i] = (t_i - abs(s_i-t_i)) / (double)s_i;
		}
	}
	double extraslope = 0.5*histogramdistance(segsolution.objhist, targethist, metric)/(double)segsolution.ssize;
	for(int i=0;i<img_w;i++)
	{
		for(int j=0;j<img_h;j++)
		{
			int n = i + j*img_w;
			if(segsolution.labeling[i][j] == BKG)
				g->add_tweights(n,0,INFTY);
			else
				g->add_tweights(n,slopes[image.colorlabel[i][j]]*w_hist+extraslope,0);
		}
	}
}

vector<double> gettargetdistribution(const Image & image, const Table2D<int> & gtimg)
{
	vector<double> targetdist(image.colorbinnum,0);
	int ssize = 0;
	for(int i=0;i<image.img_w;i++)
	{
		for(int j=0;j<image.img_h;j++)
		{
			if(gtimg[i][j]==255) // white
			{
				targetdist[image.colorlabel[i][j]]++;
				ssize++;
			}
		}
	}
	for(unsigned int i=0;i<targetdist.size();i++)
		targetdist[i] = targetdist[i] / (double)ssize;
	return targetdist;
}

#endif