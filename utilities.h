#ifndef _UTILITIES_H__
#define _UTILITIES_H__
#include <basicutil.h>
#include <consistency.h>
#include <entropy.h>
#include <iostream>  
#include <sstream>  
#include <fstream>   

// k_labeling is the known labeling
Table2D<Label> mergelabeling(GraphType * g, const Table2D<Label> & k_labeling);
// k_labeling is the known labeling
Table2D<Label> mergelabelingcompact(GraphType * g, const Table2D<Label> & k_labeling, 
	vector<Point> node_corr);
int getcompactgraph(const Image & image,const Table2D<Label> & k_labeling, Table2D<bool> & incompactgraph,
	vector<Point> & node_corr, Table2D<int> & img_corr,vector<PointPair> & compactpointpairs, 
	vector<double> & compactsmoothnesscosts);
Table2D<bool> addsmoothnessterm(GraphType * g, const Image & image, double lambda,const Table2D<Label> & k_labeling);
GraphType * getbasegraph(const Image & image, double lambda, double beta, const Table2D<int> & box);
template <typename T> string tostr(const T& t);
void objbkghist(const Image & image, vector<int> & OBJhist, vector<int> & BKGhist, 
	int & OBJsize, int & BKGsize, const Table2D<Label> & labeling);

class BreakPoint{
public:
	bool operator<(const BreakPoint& b) const   {return para < b.para;}
	double para;
	Table2D<Label> solution;
	double flow;
	int ssize; // foreground size
	double original_e; // original energy
};

// k_labeling is the known labeling
Table2D<Label> mergelabeling(GraphType * g, const Table2D<Label> & k_labeling)
{
	int img_w = k_labeling.getWidth();
	int img_h = k_labeling.getHeight();
	int n=0;
	Table2D<Label> m_labeling = k_labeling;
	for (int y=0; y<img_h; y++) 
	{
		for (int x=0; x<img_w; x++) 
		{ 
			n = x+y*img_w;
			if(k_labeling[x][y]==NONE&&g->what_segment(n) == GraphType::SOURCE)
			{
				m_labeling[x][y]=OBJ;
			}
			else if(k_labeling[x][y]==NONE&&g->what_segment(n) == GraphType::SINK)
			{
				m_labeling[x][y]=BKG;
			}
			//else
				//exit(-1);
		}
	}
	return m_labeling;
}

// k_labeling is the known labeling
Table2D<Label> mergelabelingcompact(GraphType * g, const Table2D<Label> & k_labeling, 
	vector<Point> node_corr)
{
	int img_w = k_labeling.getWidth();
	int img_h = k_labeling.getHeight();
	int n=0;
	Table2D<Label> m_labeling = k_labeling;
	int nodenum = node_corr.size();
	for (int n=0; n<nodenum; n++) 
	{
		if(g->what_segment(n) == GraphType::SOURCE)
		{
			if(m_labeling[node_corr[n]]==UNKNOWN)
				m_labeling[node_corr[n]]=OBJ;
		}
		else if(g->what_segment(n) == GraphType::SINK)
		{
			if(m_labeling[node_corr[n]]==UNKNOWN)
				m_labeling[node_corr[n]]=BKG;
		}
	}
	return m_labeling;
}

// get compact graph flag
// k_labeling is the known labeling
// return compact graph size
int getcompactgraph(const Image & image,const Table2D<Label> & k_labeling, Table2D<bool> & incompactgraph,
	vector<Point> & node_corr, Table2D<int> & img_corr,vector<PointPair> & compactpointpairs, 
	vector<double> & compactsmoothnesscosts)
{
	int img_w = image.img_w;
	int img_h = image.img_h;
	// number of neighboring pairs of pixels
	int numNeighbor = image.pointpairs.size();
	// n-link - smoothness term
	int node_id1 =0, node_id2 =0;
	incompactgraph.reset(img_w,img_h,false);
	img_corr.reset(img_w,img_h,-1);
	for(int i=0;i<numNeighbor;i++)
	{
		PointPair pp = image.pointpairs[i];
		Label l1 = k_labeling[pp.p1];
		Label l2 = k_labeling[pp.p2];
		if(!(((l1==OBJ)&&(l2==OBJ))||((l1==BKG)&&(l2==BKG))||((l1==NONE)&&(l2==NONE))))
		{
			compactpointpairs.push_back(pp);
			compactsmoothnesscosts.push_back(image.smoothnesscosts[i]);
			incompactgraph[pp.p1] = true;
			incompactgraph[pp.p2] = true;
		}
	}
	int compactsize = countintable(incompactgraph,true);
	node_corr = vector<Point>(compactsize,Point());
	int idx = 0;
	for(int j=0;j<img_h;j++)
	{
		for(int i=0;i<img_w;i++)
			if(incompactgraph[i][j])
			{
				img_corr[i][j] = idx;
				node_corr[idx++] = Point(i,j);
			}
	}
	return compactsize;
}

// add smoothness term to the graph
// lambda is the weight of the smoothness term
// k_labeling is the known labeling
Table2D<bool> addsmoothnessterm(GraphType * g, const Image & image, double lambda,const Table2D<Label> & k_labeling)
{
	int img_w = image.img_w;
	int img_h = image.img_h;
	// number of neighboring pairs of pixels
	int numNeighbor = image.pointpairs.size();
	// n-link - smoothness term
	int node_id1 =0, node_id2 =0;
	Table2D<bool> incompactgraph(img_w,img_h,false);
	for(int i=0;i<numNeighbor;i++)
	{
		PointPair pp = image.pointpairs[i];
		node_id1 = pp.p1.x+pp.p1.y*img_w;
		node_id2 = pp.p2.x+pp.p2.y*img_w;
		// if the two nodes are not given labeling
		Label l1 = k_labeling[pp.p1];
		Label l2 = k_labeling[pp.p2];
		if(!(((l1==OBJ)&&(l2==OBJ))||((l1==BKG)&&(l2==BKG))))
		{
			double v = lambda*image.smoothnesscosts[i];
			g->add_edge(node_id1,node_id2,v,v);
			incompactgraph[pp.p1] = true;
			incompactgraph[pp.p2] = true;
		}
	}
	return incompactgraph;
}

// base graph: smothness + L1 color separation + hard constraints
GraphType * getbasegraph(const Image & image, double lambda, double beta, const Table2D<int> & box)
{
	
	GraphType * g = new GraphType(/*estimated # of nodes*/ image.img_size+image.colorbinnum, 
		/*estimated # of edges*/ 5*image.img_size); 
	g->add_node(image.img_size+image.colorbinnum);    // adding nodes
	Table2D<bool> ROI(image.img_w,image.img_h,true); // Region of interest
	addsmoothnessterm(g, image, lambda,ROI);
	addl1separationterm(g, image.colorlabel, beta, ROI);
	for(int x=0;x<image.img_w;x++)  
	{
		for(int y=0;y<image.img_h;y++)
			if(box[x][y]==255)
				g->add_tweights(x+y*image.img_w,0,INFTY);// Hard constraint outside the bounding box
	}
	return g;
}

template <typename T> string tostr(const T& t) 
{ 
	ostringstream os; 
	os<<t; 
	return os.str(); 
}


void objbkghist(const Image & image, vector<int> & OBJhist, vector<int> & BKGhist, 
	int & OBJsize, int & BKGsize, const Table2D<Label> & labeling)
{
	OBJhist = vector<int>(image.colorbinnum,0);
	BKGhist = vector<int>(image.colorbinnum,0);
	OBJsize = 0;
	BKGsize = 0;
	for(int x =0;x<image.img_w;x++)
	{
		for(int y=0;y<image.img_h;y++)
		{
			if(labeling[x][y]==OBJ)
			{
				OBJsize++;
				OBJhist[image.colorlabel[x][y]]++;
			}
			else if(labeling[x][y]==BKG)
			{
				BKGsize++;
				BKGhist[image.colorlabel[x][y]]++;
			}
		}
	}
}

#endif