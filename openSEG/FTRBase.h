#pragma once
// Base class of FTR for image labeling problem

class FTRBreakPoint{
public:
	double lambda;
	Table2D<Label> solution;
	double flow;
	int ssize; // foreground size
	double actual_en; // actual energy
	double approximate_en; // approximate energy
	bool operator<(const FTRBreakPoint& b) const   {return lambda < b.lambda;}
};

class FTRBase{
public:
	FTRBase(double lambda_init_ = 0.001, double lambda_min_ = 1e-3 , double lambda_multiplier_ = 5, 
		double ratio_t_ = 0.25, int max_iteration_ = 200, double energy_delta_ = 0.1);

	// computer actual energy and store in *act_energy_p
	// also can compute approximate energy and store in *app_energy_p
	virtual void computeenergy(const Table2D<Label> & labeling, 
		double * act_energy_p, double * app_energy_p = NULL)=0;
	void setinitlabeling(const Table2D<Label> & initiallabeling_, Label nonlabel = BKG);
	// collect statistics about current labeling and update approximate model
	virtual void updatemodel()=0; 

	virtual GraphType * ftrbasegraph(double lambda_)=0;
	void setpara(double lambda_init_, double lambda_min_, double lambda_multiplier_, 
		double ratio_t_, int max_iteration_, double energy_delta_){lambda_init=lambda_init_;
		lambda_min=lambda_min_;lambda_multiplier=lambda_multiplier_;ratio_t=ratio_t_;
		max_iteration=max_iteration_;energy_delta=energy_delta_;};
	void updatelabeling(const FTRBreakPoint & bp);
	bool istrivial(const FTRBreakPoint & bp);
	double FTRBase::findlambdamax(double lambda_top);
	
	FTRBreakPoint findgivenpoint(double lambda_);
	bool optimize();

	Table2D<Label> initlabeling;
	Table2D<Label> current_labeling;

	double init_e;
	double current_e;

	int gc_count; // counter of graph cut instances
	bool showflag;
protected:
	// parameters of FTR
	double ratio_t; // ratio threshold, commonly set to 0.25
	double lambda_init;
	double lambda_min;
	double lambda_multiplier;
	int max_iteration;
	double energy_delta;

	int img_w;
	int img_h;

	Table2D<bool> ROI;
	int ROISize;

	double lambda_ftr; // current lambda

	// distance transform
	Table2D<double> dt;
};

FTRBase::FTRBase(double lambda_init_, double lambda_min_, double lambda_multiplier_, 
	double ratio_t_, int max_iteration_, double energy_delta_)
	:lambda_init(lambda_init_),lambda_min(lambda_min_),lambda_multiplier(lambda_multiplier_)
	,ratio_t(ratio_t_),max_iteration(max_iteration_),energy_delta(energy_delta_)
{
	lambda_ftr = lambda_init_;
	showflag = true;
}

void FTRBase::setinitlabeling(const Table2D<Label> & initiallabeling_, Label nonlabel)
{
	lambda_ftr = lambda_init;
	initlabeling = initiallabeling_;
	img_w = initlabeling.getWidth();
	img_h = initlabeling.getHeight();
	ROISize = img_w * img_h;
	ROI.reset(img_w,img_h,true);
	for(int y=0; y<img_h; y++)
	{
		for(int x=0; x<img_w; x++) 
		{
			if(initlabeling[x][y]==NONE)
			{
				ROI[x][y] = false;
				initlabeling[x][y] = nonlabel;
				ROISize--;
			}
		}
	}
	current_labeling = initlabeling;
	computeenergy(initlabeling,&init_e);
	current_e = init_e;
	updatemodel();
	if(showflag) outv(init_e);
	gc_count = 0;
}

void FTRBase::updatelabeling(const FTRBreakPoint & bp)
{
	current_labeling = bp.solution;
	current_e = bp.actual_en;
	updatemodel();
}

bool FTRBase::istrivial(const FTRBreakPoint & bp)
{
	return bp.ssize==0||bp.ssize==ROISize;
}

// Find lambdamax in the range [lambda_min lambda_top]
double FTRBase::findlambdamax(double lambda_top)
{
	double lambda_bot = lambda_min;
	while(1)
	{
		double lambda_middle = ( lambda_top + lambda_bot ) / 2.0;
		FTRBreakPoint bp = findgivenpoint(lambda_middle);
		if(bp.solution==current_labeling)
			lambda_top = lambda_middle;
		else
			lambda_bot = lambda_middle;
		if(abs(lambda_top - lambda_bot)<=lambda_min)
			return lambda_bot;
	}
}

FTRBreakPoint FTRBase::findgivenpoint(double lambda_)
{
	GraphType * g = ftrbasegraph(lambda_);
	FTRBreakPoint bp;
	bp.lambda = lambda_;
	bp.flow = g->maxflow();
	gc_count++;
	Table2D<Label> m_labeling(img_w,img_h,NONE);
	getlabeling(g, m_labeling);
	delete g;
	bp.solution = m_labeling;
	bp.ssize = countintableROI(m_labeling,OBJ,ROI);
	return bp;
}

bool FTRBase::optimize()
{
	int itr_id=0;
	while((itr_id++)<max_iteration)
	{
		if(showflag) cout<<"Iteration "<<itr_id<<"---------------------------------"<<endl;
		if(showflag) cout<<"Actual energy of current labeling "<<current_e<<endl;
		if(itr_id>=max_iteration)
			break;
		if(lambda_ftr > 1e+6)
			break;
		FTRBreakPoint bp = findgivenpoint(lambda_ftr);
		if(showflag) cout<<"lambda = "<<lambda_ftr<<" new labeling foreground size "<<bp.ssize<<endl;
		if(istrivial(bp))
		{
			if(showflag) outs("trivial solution!!!");
			lambda_ftr = lambda_ftr*lambda_multiplier;
			if(showflag) outv(lambda_ftr);
			continue;
		}
		
		// lambda is greater than lambda_max
		if(bp.solution == current_labeling)
		{
			if(showflag) outs("S_\lambda = S_0, find min change breakpoint");
			double lambda_max = findlambdamax(lambda_ftr);
			if(showflag) cout<<" lambda_max "<<lambda_max<<endl;
			bp = findgivenpoint(lambda_max);
			if(istrivial(bp))
				break; // get trivial solution
			else if(bp.solution == current_labeling)
				break;
			else
			{
				//if(showflag) outs("next_labeling is not equal to current_labeling");
				computeenergy(bp.solution,&(bp.actual_en),&(bp.approximate_en));
				if(current_e - bp.actual_en > energy_delta)
				{
					lambda_ftr = lambda_max;
					double actual_decrease = current_e- bp.actual_en;
					double predicted_decrease = current_e - bp.approximate_en;
					if(showflag) outv(actual_decrease);
					if(showflag) outv(predicted_decrease);
					if(showflag) outv(actual_decrease/predicted_decrease);
					if(showflag) cout<<"actual reduction / predicted reduction "<<actual_decrease/predicted_decrease<<endl;
					// adjust trust region size
					if(actual_decrease/predicted_decrease>ratio_t)
					{
						lambda_ftr = lambda_ftr/lambda_multiplier;
						if(showflag) cout<<"Increase step size, set lambda = "<<lambda_ftr<<endl;
					}
					else
					{
						lambda_ftr = lambda_ftr*lambda_multiplier;
						if(showflag) cout<<"Reduce step size, set lambda = "<<lambda_ftr<<endl;
					}

					updatelabeling(bp);
					continue;
				}
				else
				{
					if(showflag) outs("local minima!");
					break;
				}
			}
		}

		computeenergy(bp.solution,&(bp.actual_en),&(bp.approximate_en));

		double actual_decrease = current_e- bp.actual_en;
		double predicted_decrease = current_e - bp.approximate_en;

		if(showflag) outv(actual_decrease);
		if(showflag) outv(predicted_decrease);
		if(showflag) outv(actual_decrease/predicted_decrease);
		
		// update labeling
		if(actual_decrease>energy_delta)
		{
			updatelabeling(bp);
		}
		else
		{
			if(showflag) outs("Did not decrease energy/////");
			// update stepsize
			lambda_ftr = lambda_ftr*lambda_multiplier;
			if(showflag) cout<<"Reduce step size, set lambda = "<<lambda_ftr<<endl;
			continue;
		}

		// adjust trust region size
		if(actual_decrease/predicted_decrease>ratio_t)
		{
			lambda_ftr = lambda_ftr/lambda_multiplier;
			if(showflag) cout<<"Increase step size, set lambda = "<<lambda_ftr<<endl;
		}
		else
		{
			lambda_ftr = lambda_ftr*lambda_multiplier;
			if(showflag) cout<<"Reduce step size, set lambda = "<<lambda_ftr<<endl;
		}
	}
	if(current_e<init_e-energy_delta)
		return true;
	else
		return false;
}