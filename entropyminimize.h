#ifndef _EntropyMINIMIZE_H__
#define _EntropyMINIMIZE_H__
#include <BCD.h>
//#include "PPBCEntropy.h"

void EntropyMinimize(int argc, char * argv[]);

// for Entropy-base image segmentation 
void EntropyMinimize(int argc, char * argv[])
{
	char * rootdir = "E:/eccv14_PBO/PBO_code/test/";
	// parse arguments and set energy parameters
	char * imgname = argv[1];
	int numbin = atoi(argv[2]);
	double w_smooth = atof(argv[3]); // weight of smoothness term

	// load original image
	Image image(to_Cstr(rootdir<<"images/"<<imgname<<".bmp"),
		imgname,256/numbin,8); // 8 connect smoothness term
	image.print();
	int img_w = image.img_w;
	int img_h = image.img_h;
	srand( (unsigned int) time( NULL ) ); // RANDOM NUMBER INITIALIZER
	clock_t start,finish; // Timing

	// load bounding box
	Table2D<int> boximg = loadImage<RGB>(to_Cstr(rootdir<<"images/"<<imgname<<"_box.bmp"));
	int boxsize = countintable(boximg,0);
	outv(boxsize);
	Table2D<Label> initlabeling(img_w,img_h,NONE);
	Table2D<Label> hardconstraints(img_w,img_h,NONE);
	for(int x=0;x<img_w;x++){
		for(int y=0;y<img_h;y++){
			if(boximg[x][y]==0)//inside the box
				initlabeling[x][y] = OBJ;
			else{//outside the box
				initlabeling[x][y] = BKG;
				hardconstraints[x][y] = BKG;
			}                                      //这里，用输入的框框图去给initlabeling和hardconstraints赋值
		}                                          //现在得到的结果是：initlabeling只有OBJ和BKG; hardconstraints有BKG和NONE
	}
	Table2D<int> gt = loadImage<RGB>(to_Cstr(rootdir<<"images/"<<imgname<<"_groundtruth.bmp")); // ground truth


	////////////////////////////////////////////////////////////////////
	//用grabcut的方法做分割，没有用GMM，用的是直方图
	
	cout<<"\n\n -----------\n"<<"optimize with Block-coordinate-descent (GrabCut)"<<endl;
	start = clock();
	BCD bcd(image,w_smooth);
	bcd.initlabeling(initlabeling);
	bcd.hardconstraints = hardconstraints;
	bcd.optimize();
	finish = clock();
	double grabcuttime = (double)(finish-start)/CLOCKS_PER_SEC;


	savebinarylabeling(image.img, bcd.current_labeling,
		to_Cstr(rootdir<<"images/"<<imgname<<"_BCD_withbox.bmp"));

	// energies
	cout<<"\n\n*********\n";
	cout<<"BCD energy: "<<bcd.current_e<<endl;

	// error rate
	double errorrate = geterrorrate(bcd.current_labeling, gt, boxsize);
	cout<<"BCD error rate: %"<<errorrate<<"%"<<endl;

	// timing
	cout<<"BCD takes "<<grabcuttime<<" seconds!"<<endl;

	return;
}

#endif