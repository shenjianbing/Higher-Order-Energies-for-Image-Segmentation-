// base class for parametric maxflow
#ifndef _PPBCBASE_H__
#define _PPBCBASE_H__

#include "utilities.h"
#include <time.h>
#include <queue>

// Interval to be explored in parametric maxflow
struct ParaInterval{
	BreakPoint bplow; // breakpoint for lower parameter
	BreakPoint bpup; // breakpoint for upper parameter
};

// UNKNOWN is introduced for monotonic parametric maxflow
// In this case, NONE means out of region of interest.
struct UnknownRegion{
	Table2D<Label> knownlabeling; // fixed labeling; NONE is out of ROI. UNKNOWN is not known.
	int unknownsize;
};

class PPBCBase{
public:
	PPBCBase(double paramin_ = -12, double paramax_ = 4, double paradelta_ = 0.1
		,int pixeldelta_ = 20, bool monotonicflag_ = false);
	void setpara(double paramin_, double paramax_, double paradelta_,int pixeldelta_, bool monotonicflag_)
	{paramin=paramin_;paramax=paramax_;paradelta=paradelta_;pixeldelta=pixeldelta_;monotonicflag=monotonicflag_;};
	virtual double computeenergy(const Table2D<Label> & labeling)=0;
	void setinitlabeling(const Table2D<Label> & initiallabeling_,Label nonlabel = BKG);
	virtual void updatemodel()=0; // collect statistics about current labeling
	virtual void * parabasegraph(double para_, UnknownRegion * unknownregion_p = NULL, 
		double * flowoffset_p = NULL, vector<Point> * node_corr_p = NULL)=0;
	BreakPoint explorepara(double para_, UnknownRegion * unknownregion_p = NULL);
	virtual double getnewpara(ParaInterval paraInterval)=0;
	UnknownRegion getkonwnregion(const ParaInterval & paraInterval);
	void AddBreakPoint(BreakPoint & bp);
	void explore();
	BreakPoint gridsearch(double delta = 0);
	void savesolutions(const Image & image);
	BreakPoint SelectBestBP();
	Point mappoint(int idx){return Point(idx%img_w,idx/img_w);}
	void emplyROImask(Table2D<Label> & labeling);

protected:
	double paramin;
	double paramax;
	double paradelta;
	int pixeldelta;
	bool monotonicflag;

	int img_w;
	int img_h;

	priority_queue<BreakPoint> bps; // breakpoints
	Table2D<Label> initlabeling;
	Table2D<bool> ROI;
	int ROISize;
	double init_e;
private:
	BreakPoint best_bp;// the breakpoint with minimum original energy
};

PPBCBase::PPBCBase(double paramin_ , double paramax_, double paradelta_
	,int pixeldelta_, bool monotonicflag_)
{
	paramin = paramin_;
	paramax = paramax_;
	paradelta = paradelta_;
	pixeldelta = pixeldelta_;
	monotonicflag = monotonicflag_;
	best_bp.original_e = 1e+9;
}

void PPBCBase::setinitlabeling(const Table2D<Label> & initiallabeling_,Label nonlabel)
{
	initlabeling = initiallabeling_;
	img_w = initlabeling.getWidth();
	img_h = initlabeling.getHeight();
	ROI.reset(img_w,img_h,true);
	ROISize = img_w * img_h;
	for(int y=0; y<img_h; y++)
	{
		for(int x=0; x<img_w; x++) 
		{ 
			if(initlabeling[x][y]==NONE)
			{
				ROI[x][y] = false;
				ROISize--;
			}
		}
	}
	outv(ROISize);
	init_e = computeenergy(initlabeling);
	updatemodel();   //用现有的标记去更新前背景直方图。这个函数属于BCD类
	outv(init_e);
	bps = priority_queue<BreakPoint>();
	best_bp.original_e = 1e+9;
}

BreakPoint PPBCBase::explorepara(double para_, UnknownRegion * unknownregion_p)
{
	Table2D<Label> m_labeling(img_w,img_h,NONE);
	BreakPoint breakpoint;
	breakpoint.para = para_;

	//cout<<"explore para "<<para_<<endl;
	double flowoffset = 0;
	vector<Point> node_corr;
	clock_t start = clock();
	GraphType * g = (GraphType *)parabasegraph(para_, unknownregion_p, &flowoffset, &node_corr);
	clock_t finish = clock();
	breakpoint.flow = g->maxflow();
	breakpoint.flow+=flowoffset;
	if((unknownregion_p==NULL))
		getlabeling(g, m_labeling);
	else
		m_labeling = mergelabelingcompact(g, unknownregion_p->knownlabeling, node_corr);
	delete g;

	if(ROISize!=img_w*img_h) emplyROImask(m_labeling);
	breakpoint.solution = m_labeling;
	breakpoint.ssize = countintable(m_labeling,OBJ);
	return breakpoint;
}

void PPBCBase::emplyROImask(Table2D<Label> & labeling)
{
	for(int i=0;i<img_w;i++)
	{
		for(int j=0;j<img_h;j++)
		{
			if(ROI[i][j]==NONE)
				labeling[i][j]=NONE;
		}
	}
}
UnknownRegion PPBCBase::getkonwnregion(const ParaInterval & paraInterval)
{
	UnknownRegion unknownregion;
	unknownregion.knownlabeling.reset(img_w,img_h,NONE);
	int unknownsize = 0;
	for(int i=0;i<img_w;i++)
	{
		for(int j=0;j<img_h;j++)
		{
			if(ROI[i][j])   // If in ROI (Region of interest)
			{
				Label uplabel = paraInterval.bpup.solution[i][j];
				Label lowlabel = paraInterval.bplow.solution[i][j];
				if((uplabel==lowlabel))
					unknownregion.knownlabeling[i][j]=uplabel;
				else
				{
					unknownregion.knownlabeling[i][j]=UNKNOWN;
					unknownsize++;
				}
			}
		}
	}
	unknownregion.unknownsize = unknownsize;
	return unknownregion;
}

inline void PPBCBase::AddBreakPoint(BreakPoint & bp)
{
	bp.original_e = computeenergy(bp.solution);
	//outv(bp.original_e);
	bps.push(bp);
	if(bp.original_e<best_bp.original_e)
		best_bp = bp;
}
void PPBCBase::explore()
{
	vector<ParaInterval> intervals;
	ParaInterval interval;
	BreakPoint bpmin = explorepara(paramin);   //就是pPBC算法中的 \lambda^min
	BreakPoint bpmax = explorepara(paramax);   //                \lamnda^max
	if(bpmin.solution==bpmax.solution)
	{
		AddBreakPoint(bpmin);
		return;
	}
	else
	{
		AddBreakPoint(bpmin);
		AddBreakPoint(bpmax);
		interval.bplow = bpmin;
		interval.bpup = bpmax;
		intervals.push_back(interval);
	}
	while(!intervals.empty())
	{
		interval = intervals[intervals.size()-1];
		intervals.pop_back();
		if(abs(interval.bplow.para-interval.bpup.para)<paradelta)
			continue;
		double newpara = getnewpara(interval);
		if((newpara>interval.bpup.para+paradelta)||((newpara<interval.bplow.para-paradelta)))
		{
			cout<<"Invalid new para!!!!!!!!!!!!!!!!!!! "<<endl;
			continue;
		}
		else if(!((newpara>=(interval.bplow.para+paradelta))&&(newpara<=(interval.bpup.para-paradelta))))
		{
			//cout<<"new para: "<<newpara<<" in "<<interval.bplow.para<<" "<<interval.bpup.para<<" same as bound"<<endl;
			continue;
		}
		//cout<<"new para: "<<newpara<<" in "<<interval.bplow.para<<" "<<interval.bpup.para<<endl;
		BreakPoint newbp;
		if(monotonicflag==false)
			newbp = explorepara(newpara);
		else
		{
			UnknownRegion unknownregion = getkonwnregion(interval);
			newbp = explorepara(newpara,&unknownregion);
		}
		bool equaltolow = (abs(newbp.ssize-interval.bplow.ssize)<pixeldelta);
		bool equaltoup = (abs(newbp.ssize-interval.bpup.ssize)<pixeldelta);
		if(equaltolow||equaltoup)
			continue;
		AddBreakPoint(newbp);
		//if(!(newbp.solution==interval.bplow.solution))
		if(abs(newbp.ssize-interval.bplow.ssize)>pixeldelta)
		{
			ParaInterval newinterval = interval;
			newinterval.bpup = newbp;
			intervals.push_back(newinterval);
		}
		//if(!(newbp.solution==interval.bpup.solution))
		if(abs(newbp.ssize-interval.bpup.ssize)>pixeldelta)
		{
			ParaInterval newinterval = interval;
			newinterval.bplow = newbp;
			intervals.push_back(newinterval);
		}
	}
}

BreakPoint PPBCBase::gridsearch(double delta)
{
	BreakPoint bestbp;
	double best_e = INFTY;
	if(delta==0)delta = paradelta;
	for(double para_ = paramin;para_<=paramax+1e-5;para_=para_+delta)
	{
		BreakPoint bp = explorepara(para_);
		bp.original_e = computeenergy(bp.solution);
		if(bp.original_e<best_e-(1e-4))
		{
			best_e = bp.original_e;
			bestbp = bp;
		}
	}
	return bestbp;
}
void PPBCBase::savesolutions(const Image & image)
{
	priority_queue<BreakPoint> bps_copy(bps);
	int i=0;
	while(!bps_copy.empty())
	{
		BreakPoint bp = bps_copy.top();
		bps_copy.pop();
		string outname = "Entropy/multilabel/expansion/breakpoints/";
		outname +=image.imgname;
		outname +="_result_";
		outname +=tostr(bp.para);
		outname += ".bmp";
		savebinarylabeling(image.img, bp.solution,outname);
	}
}

BreakPoint PPBCBase::SelectBestBP()
{
	cout<<"number of breakpoints: "<<bps.size()<<endl;
	return best_bp;
}

#endif