// high-order consistency (color separation)
#pragma once

// add L1 color separation term to the graph
// ROI is the region of interest
// separation_w is the weight of the color separation term
void addl1separationterm(GraphType * g, const Table2D<int> &colorlabel,double separation_w, const Table2D<bool> & ROI);
int getl1penalty(Table2D<int> & colorlabel,Table2D<Label> & labeling);
Table2D<double> readtxtfile(char * txtfilename, int row,int col);
void addJSseparationterm(GraphType * g, const Table2D<int> &colorlabel,double separation_w, 
	const Table2D<bool> & ROI, int numbrk, const vector<int> & ROIhist, const vector<double> & thetamaxes, const vector<double> & slopes);

void addl1separationterm(GraphType * g, const Table2D<int> &colorlabel,double separation_w, const Table2D<bool> & ROI)
{
	int node_id = 0;
	int img_w = colorlabel.getWidth();
	int img_h = colorlabel.getHeight();
	for(int y=0; y<img_h; y++) // adding links to auxiliary nodes
	{
		for(int x=0; x<img_w; x++) 
		{ 
			if(ROI[x][y])
			{
				node_id = x+y*img_w;
				g->add_edge( node_id, colorlabel[x][y]+img_w*img_h,separation_w, separation_w);
			}
		}
	}
}

// Color histogram overlap of box and its outside region based on L1 metric
int getl1penalty(Table2D<int> & colorlabel,Table2D<Label> & labeling)
{
	int returnv = 0;
	int bin_num = colorlabel.getMax()+1;
	vector<int> obj_vector(bin_num,0);
	vector<int> bkg_vector(bin_num,0);
	for(unsigned int j=0;j<labeling.getHeight();j++)
	{
		for(unsigned int i=0;i<labeling.getWidth();i++)
		{
			if(labeling[i][j] == OBJ)
			{
				obj_vector[colorlabel[i][j]]++;
			}
			else
				bkg_vector[colorlabel[i][j]]++;
		}
	}
	for(int i=0;i<bin_num;i++)
	{
		returnv+=min(obj_vector[i],bkg_vector[i]);
	}
	return returnv;
}

Table2D<double> readtxtfile(char * txtfilename,int row, int col)
{
	Table2D<double> data(col,row,0);
	FILE *f = fopen(txtfilename, "r");
	float tmp;
	if(f) 
	{
	  for (int y=0; y<row; y++) 
	  {
		  for (int x=0; x<col; x++) 
		  {
			  fscanf(f, "%f ", &tmp);
			  cout<<tmp<<' ';
			  data[x][y] = tmp;
		  }
		  cout<<'\n';
	  }
	  cout<<"Reading data - done"<<endl;
	}  
	else 
	{
		cout << "Unable to open data file"<<endl;
		exit(-1);
	}
	fclose(f);
	return data;
}

void addJSseparationterm(GraphType * g, const Table2D<int> &colorlabel,double separation_w, 
	const Table2D<bool> & ROI, int numbrk, const vector<int> & ROIhist, const vector<double> & thetamaxes, const vector<double> & slopes)
{
	int node_id = 0;
	int img_w = colorlabel.getWidth();
	int img_h = colorlabel.getHeight();
	double theta =0, slope = 0;
	int colorbinnum = ROIhist.size();
	for(int y=0; y<img_h; y++) // adding links to auxiliary nodes
	{
		for(int x=0; x<img_w; x++) 
		{ 
			if(ROI[x][y])
			{
				node_id = x+y*img_w;
				for(int n=0;n<numbrk;n++)
				{
					slope = slopes[n];
					g->add_edge( node_id, img_w*img_h+colorlabel[x][y]+n*2*colorbinnum,0, separation_w *slope);
					g->add_edge( node_id, img_w*img_h+colorlabel[x][y]+colorbinnum+n*2*colorbinnum,separation_w*slope,0);
				}
			}
		}
	}
	
	for(int i=0;i<colorbinnum;i++)
	{
		for(int n=0;n<numbrk;n++)
		{
			theta = thetamaxes[n];
			if(ROIhist[i]!=0)
				g->add_edge(img_w*img_h+i+n*2*colorbinnum,img_w*img_h+colorbinnum+i+n*2*colorbinnum,0,(ROIhist[i])*separation_w/2.0*theta);
		}
	}
	//cout<<"add multi finished"<<endl;
}