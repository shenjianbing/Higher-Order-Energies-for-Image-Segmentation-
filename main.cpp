//higher-order image segmentation
//paper: A New Energy Minimization Method for Higher-Order Binary Energy Functions and Its Application to Better Image Segmentation

#include "Image.h"
#include "utilities.h"

#include <BCD.h>

int main(int argc, char * argv[])
{
	clock_t start = clock() ,finish ;

	char * rootdir = "./test/images/"; 

	char * imgname = argv[1] ;

	int numbin = atoi(argv[2]);  //16 bins
	double w_smooth = atof(argv[3]); // weight of smoothness term
	w_smooth =   15 ;


	// load original image
	Image image( to_Cstr(rootdir<<imgname<<".bmp"),
		         imgname , 256/numbin , 8 ); // 8 connect smoothness term

	//image.print();
	int img_w = image.img_w;
	int img_h = image.img_h;
	srand( (unsigned int) time( NULL ) ); // RANDOM NUMBER INITIALIZER


	// load bounding box
	Table2D<int> boximg = loadImage<RGB>(to_Cstr(rootdir<<imgname<<"_box.bmp"));

	int boxsize = countintable(boximg,128); 
	//outv(boxsize);

	//initialize two matrixes
	Table2D<Label> initlabeling(img_w,img_h,NONE);
	Table2D<Label> hardconstraints(img_w,img_h,NONE);
	for(int x=0;x<img_w;x++){
		for(int y=0;y<img_h;y++){
			if(boximg[x][y]>=100)//inside the box
				initlabeling[x][y] = OBJ;
			else{//outside the box
				initlabeling[x][y] = BKG;
				hardconstraints[x][y] = BKG;
			}                                      
		}                                          
	}
	Table2D<int> gt = loadImage<RGB>(to_Cstr(rootdir<<imgname<<"_groundtruth.bmp")); // ground truth


	////////////////////////////////////////////////////////////////////
	start = clock();
	BCD bcd(image,w_smooth);
	bcd.initlabeling(initlabeling);  

	bcd.hardconstraints = hardconstraints;  
	
	bcd.optimize();

	savebinarylabeling(image.img, bcd.current_labeling,
		to_Cstr(rootdir<<imgname<<"_OUR.bmp"));


	//output£º energy, error , time 
	finish = clock();
	double grabcuttime = (double)(finish-start)/CLOCKS_PER_SEC;
	double errorrate = geterrorrate(bcd.current_labeling, gt, boxsize);


	cout<<"Error rate:  "<<errorrate*100<<"%"<<endl;
	cout<<"Computational time: "<<grabcuttime<<" seconds!"<<endl;

	system("pause");
	return -1;
}