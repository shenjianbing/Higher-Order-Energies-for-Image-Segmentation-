#ifndef _MATCHDIST_H__
#define _MATCHDIST_H__

#include "Image.h"
#include "utilities.h"
#include "boundutil.h"
#include "PPBCDist.h"
#include "FTRDist.h"

void MatchDist(int argc, char * argv[]);

// matching target distribution
void MatchDist(int argc, char * argv[])
{
	char * rootdir = "D:/repositories/PBO_public_code/";
	//energy parameter
	char * imgname = argv[1];
	int numbin = atoi(argv[2]);
	METRIC metric;
	if(0==strcmp(argv[3],"KL")) metric= KLDIVERGENCE;
	else if(0==strcmp(argv[3],"BHA")) metric = BHATTACHARYYA;
	double w_smooth = atof(argv[4]); // weight of smoothness term

	double w_dist = atof(argv[5]);
	outv(w_smooth);
	outv(w_dist);
	double errorrate;

	//initialization - Image, target distribution, initial solution
	Image image = Image(to_Cstr(rootdir<<"images/"<<imgname<<".bmp"),imgname,256/numbin,8);
	image.print();
	int img_w = image.img_w;
	int img_h = image.img_h;

	Table2D<int> gtimg = loadImage<RGB>(to_Cstr(rootdir<<"images/"<<imgname<<"_groundtruth.bmp"));
	Table2D<int> initimg = loadImage<RGB>(to_Cstr(rootdir<<"images/"<<imgname<<"_box.bmp"));

	vector<double> targetdist = gettargetdistribution(image,gtimg);
	//initimg.reset(0); // can also start from trivial solution (whole image as foreground)
	Table2D<Label> initlabeling = getinitlabeling(initimg,0);

	
	// FTR
	FTRDist ftr = FTRDist(w_smooth,w_dist,image,targetdist,metric);
	ftr.setinitlabeling(initlabeling);
	ftr.optimize();
	savebinarylabeling(image.img, ftr.current_labeling,
		to_Cstr(rootdir<<"images/"<<imgname<<"_"<<(metric==KLDIVERGENCE?"KL":"BHA")<<"_FTR.bmp"));
	errorrate = geterrorrate(ftr.current_labeling, gtimg, countintable(initimg,0));
	cout<<"FTR error rate: %"<<errorrate*100<<"%"<<endl;
	cout<<"FTR energy: $"<<ftr.current_e<<"$"<<endl;

	// Auxiliary Cut
	Table2D<Label> currentlabeling = initlabeling;
	PPBCDist PPBCdist(w_smooth,w_dist,image,targetdist,metric);
	PPBCdist.setpara(-4,4,0.01,10,true);
	double best_e = INFTY;
	while(1)
	{
		PPBCdist.setinitlabeling(currentlabeling);
		BreakPoint bp = PPBCdist.explorepara(0);
		bp.original_e = PPBCdist.computeenergy(bp.solution);
		if(bp.original_e<best_e-1)
		{
			currentlabeling = bp.solution;
			best_e = bp.original_e;
			cout<<"best breakpoint para: "<<bp.para<<endl;
		}
		else
			break;
	}
	savebinarylabeling(image.img, currentlabeling,
		to_Cstr(rootdir<<"images/"<<imgname<<"_"<<(metric==KLDIVERGENCE?"KL":"BHA")<<"_AUX.bmp"));
	errorrate = geterrorrate(currentlabeling, gtimg, countintable(initimg,0));
	cout<<"Auxiliary Cut error rate: %"<<errorrate*100<<"% / 100"<<endl;
	cout<<"Auxiliary Cut energy: $"<<best_e<<"$"<<endl;

	// PPBC
	currentlabeling = initlabeling;
	best_e = INFTY;
	while(1)
	{
		PPBCdist.setinitlabeling(currentlabeling);
		PPBCdist.explore();
		BreakPoint bp = PPBCdist.SelectBestBP();
		if(bp.original_e<best_e-1)
		{
			currentlabeling = bp.solution;
			best_e = bp.original_e;
			cout<<"best breakpoint para: "<<bp.para<<endl;
		}
		else
			break;
	}
	savebinarylabeling(image.img, currentlabeling,
		to_Cstr(rootdir<<"images/"<<imgname<<"_"<<(metric==KLDIVERGENCE?"KL":"BHA")<<"_PBO.bmp"));

	errorrate = geterrorrate(currentlabeling, gtimg, countintable(initimg,0));
	cout<<"PBO error rate: %"<<errorrate*100<<"% / 100"<<endl;
	cout<<"PBO energy: $"<<best_e<<"$"<<endl;
}

#endif