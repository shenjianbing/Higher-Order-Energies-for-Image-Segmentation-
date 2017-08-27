#ifndef _FTREntropy_H__
#define _FTREntropy_H__
#include <FTRBase.h>
class FTREntropy: public FTRBase{
public:
	FTREntropy(const Image & image_, double w_smooth_, double w_bits_, int num_comp_);
	virtual void computeenergy(const Table2D<Label> & labeling, 
		double * act_energy_p, double * app_energy_p);
	// collect statistics about current labeling and update approximate model
	virtual void updatemodel(); 
	virtual GraphType * ftrbasegraph(double lambda_);
	double colorseparationcost(const Table2D<Label> & labeling) const;
	double trivialsolutionenergy()const;
private:
	double w_smooth;
	double w_bits;
	// for color separation term
	int num_comp; // number of components for approximating JS color separation term
	vector<double> thetamaxes, slopes;
	Image image;

	vector<int> ROIhist;
	double JSoffset; // n_k*log(n_k) sum over k
	// appriximate model: a*|s| + b
	double app_a;
	double app_b;

	Table2D<double> dtOBJ;
	Table2D<double> dtBKG;
	Table2D<double> dt;
};

FTREntropy::FTREntropy(const Image & image_, double w_smooth_, double w_bits_, int num_comp_)
		:image(image_),w_smooth(w_smooth_),w_bits(w_bits_),num_comp(num_comp_)
{
	cout<<"number of breakpoints: "<<num_comp<<endl;
	char txtfilename[100];
	strcpy(txtfilename,"D:/repositories/openSEG/JS");
	strcat(txtfilename, "/thetamaxslope");
	char anumbrk[3];
	itoa(num_comp, anumbrk, 10);
	strcat(txtfilename, anumbrk);
	strcat(txtfilename, ".txt");
	cout<<txtfilename<<endl;
	Table2D<double> thetamaxslopetable = readtxtfile(txtfilename,num_comp,2);
	thetamaxes = vector<double>(num_comp,0);
	slopes = vector<double>(num_comp,0);
	for(int i=0;i<num_comp;i++)
	{
		thetamaxes[i] = thetamaxslopetable[0][i];
		slopes[i] = thetamaxslopetable[1][i];
	}
}

void FTREntropy::computeenergy(const Table2D<Label> & labeling, 
		double * act_energy_p, double * app_energy_p)
{
	if(ROIhist.size()==0)
	{
		// ROI compacthist
		ROIhist = vector<int>(image.colorbinnum,0);
		for(int x =0; x< img_w;x++)
		{
			for(int y=0;y<img_h;y++)
			{
				if(ROI[x][y])
				{
					ROIhist[image.colorlabel[x][y]]++;
				}
			}
		}
		JSoffset = 0;
		for(int i=0;i<image.colorbinnum;i++)
		{
			if(ROIhist[i]) JSoffset += ROIhist[i]*log((double)ROIhist[i])/log(2.0);
		}
	}
	double smooth_cost = getsmoothnesscost(image,labeling);
	int objsize = countintableROI(labeling,OBJ,ROI);
	int bkgsize = countintableROI(labeling,BKG,ROI);
	double volume_cost = objsize*log(max((double)objsize,1.0))+bkgsize*log(max((double)bkgsize,1.0));
	volume_cost = volume_cost/log(2.0);
	double separation_cost = colorseparationcost(labeling);
	separation_cost -= JSoffset;
	double act_energy = w_smooth*smooth_cost+w_bits*(volume_cost+separation_cost);
	*act_energy_p = act_energy;
	if(app_energy_p)
	{
		double app_volume_cost = app_a*objsize + app_b;
		double app_energy = w_smooth*smooth_cost+w_bits*(app_volume_cost+separation_cost);
		*app_energy_p = app_energy;
	}
}
void FTREntropy::updatemodel()
{
	int S0 = countintableROI(current_labeling,OBJ,ROI);
	S0 = max(1,S0); // make sure S0 is not zero
	S0 = min(ROISize-1,S0); // make sure S0 is not ROI_size
	app_a = log((double)S0/(ROISize-S0))/log(2.0);
	app_b = (S0*log((double)S0)+(ROISize-S0)*log((double)ROISize-S0))/log(2.0)-S0*app_a;

	dt = getDistanceTransform(replacelabeling(current_labeling, OBJ));
	//dtOBJ =  getDistanceTransform(replacelabeling(current_labeling, OBJ));
	//dtBKG =  getDistanceTransform(replacelabeling(current_labeling, BKG));
	//dt = dtOBJ;
	/*for(int i=0;i<img_w;i++)
	{
		for(int j=0;j<img_h;j++)
		{
			if(current_labeling[i][j]==BKG)
				dt[i][j] = - dtBKG[i][j];
		}
	}*/
}
GraphType * FTREntropy::ftrbasegraph(double lambda_)
{
	int colorbinnum = image.colorbinnum;
	int N = img_w*img_h;

	GraphType * g = new GraphType(/*estimated # of nodes*/ N+2*num_comp*colorbinnum, 
		/*estimated # of edges*/ (4+2*num_comp)*N); 
	g->add_node(N+2*num_comp*colorbinnum);    // adding nodes
	addsmoothnessterm(g, image, w_smooth, ROI);
	for(int i=0;i<img_w;i++)
	{
		for(int j=0;j<img_h;j++)
		{
			int n = i+j*img_w;
			g->add_tweights(n,dt[i][j]*lambda_,app_a);
		}
	}
	addJSseparationterm(g, image.colorlabel,w_bits*2, ROI, num_comp, ROIhist, thetamaxes, slopes);

	return g;
}

double FTREntropy::colorseparationcost(const Table2D<Label> & labeling) const
{	
	int colorbinnum = image.colorbinnum;
	vector<int> OBJhist(colorbinnum,0);
	vector<int> BKGhist(colorbinnum,0);
	for(int i=0;i<img_w;i++)
	{
		for(int j=0;j<img_h;j++)
		{
			if((labeling[i][j]==OBJ)&&(ROI[i][j]))
				OBJhist[image.colorlabel[i][j]]++;
			else if((labeling[i][j]==BKG)&&(ROI[i][j]))
				BKGhist[image.colorlabel[i][j]]++;
		}
	}

	double returnv = 0.0;
	double theta=0,slope = 0;
	for(int k=0;k<colorbinnum;k++)
	{
		int nk = ROIhist[k];
		if(nk)
		{
			int nsk = OBJhist[k];
			int nsbark = BKGhist[k];
			for(int i=0;i<num_comp;i++)
			{
				theta = thetamaxes[i];
				slope = slopes[i];
				double v = min(theta*nk,min(slope*nsk*2,(slope*nsbark*2)));;
				returnv += v;
			}
		}
	}
	return returnv;
}

double FTREntropy::trivialsolutionenergy()const
{
	return ROISize * log(max(1.0,(double)ROISize))/log(2.0);
}

#endif