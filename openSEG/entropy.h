#ifndef _ENTROPY_H__
#define _ENTROPY_H__

double getgrabcutenergy(const Image & image, double w_bits, double w_smooth, const Table2D<Label> & m_labeling);
double getmultienergy(const Image & image, double w_bits, double lambda, 
	const Table2D<unsigned char> & labeling,int maxnumSeg,double labelcost = 0);

// 得到grabcut能量函数的能量值!
//                                     ￥￥￥***这就是我们要优化的公式***￥￥￥
// smoothness cost + number of bits
// ignore smoothness penalty between OBJ and NONE or BKG and NONE.
double getgrabcutenergy(const Image & image, double w_bits, double w_smooth, const Table2D<Label> & m_labeling)
{
	int img_w = image.img_w;
	int img_h = image.img_h;
	double smoothenergy = 0;
	// number of neighboring pairs of pixels
	int numNeighbor = image.pointpairs.size();
	// n-link - smoothness term
	int node_id1 =0, node_id2 =0;
	for(int i=0;i<numNeighbor;i++)
	{
		PointPair pp = image.pointpairs[i];
		node_id1 = pp.p1.x+pp.p1.y*img_w;
		node_id2 = pp.p2.x+pp.p2.y*img_w;
		if((m_labeling[pp.p1]!=NONE)&&(m_labeling[pp.p2]!=NONE))
		{
			//double v = fn(dI(image.img[pp.p1],image.img[pp.p2]),lambda,image.sigma_square)/(pp.p1-pp.p2).norm();  
			if(m_labeling[pp.p1]!=m_labeling[pp.p2])
			{
				double v = w_smooth*image.smoothnesscosts[i];
				smoothenergy += v;
			}
		}
	}
	//printf("smoothenergy = %f\n", smoothenergy);

	//estimate appearance model
	vector<int> OBJHIST(image.colorbinnum,0);
	vector<int> BKGHIST(image.colorbinnum,0);
	int OBJsize=0, BKGsize = 0;
	for(int j=0;j<img_h;j++)
	{
		for(int i=0;i<img_w;i++)
		{
			if(m_labeling[i][j]==OBJ)
			{
				OBJsize++;
				OBJHIST[image.colorlabel[i][j]] ++;
			}
			else if(m_labeling[i][j]==BKG)
			{
				BKGsize++;
				BKGHIST[image.colorlabel[i][j]] ++;
			}
		}
	}

	double bits =0;
	for(int i=0;i<image.colorbinnum;i++)
	{
		if(OBJHIST[i]!=0)
			bits+=-(double)OBJHIST[i]*log(OBJHIST[i]/(double)OBJsize)/log(2.0);
		if(BKGHIST[i]!=0)
			bits+=-(double)BKGHIST[i]*log(BKGHIST[i]/(double)BKGsize)/log(2.0);
	}
	//printf("smoth + bits = %f\n",smoothenergy+bits);
	return smoothenergy+bits*w_bits;
}

// smoothness cost + number of bits for multilabel
double getmultienergy(const Image & image, double w_bits, double lambda, 
	const Table2D<unsigned char> & labeling,int maxnumSeg,double labelcost)
{
	int img_w = image.img_w;
	int img_h = image.img_h;
	double smoothenergy = 0;
	// number of neighboring pairs of pixels
	int numNeighbor = image.pointpairs.size();
	// n-link - smoothness term
	int node_id1 =0, node_id2 =0;
	for(int i=0;i<numNeighbor;i++)
	{
		PointPair pp = image.pointpairs[i];
		node_id1 = pp.p1.x+pp.p1.y*img_w;
		node_id2 = pp.p2.x+pp.p2.y*img_w;
		if(labeling[pp.p1]!=labeling[pp.p2])
			smoothenergy += fn(dI(image.img[pp.p1],image.img[pp.p2]),lambda,image.sigma_square)/(pp.p1-pp.p2).norm();   
	}
	//printf("smoothenergy = %f\n", smoothenergy);
	vector<int> seghist(maxnumSeg,0); // histogram for each segment
	vector<vector<int>> segcolorhist(maxnumSeg,vector<int>(image.colorbinnum,0)); // histogram for each segment and color bin
	
	for(int i=0;i<img_w;i++)
	{
		for(int j=0;j<img_h;j++)
		{
			int l = labeling[i][j];
			int c = image.colorlabel[i][j];
			seghist[l]++;
			segcolorhist[l][c]++;
		}
	}

	double bits =0;
	for(int l=0;l<maxnumSeg;l++)
	{
		if(seghist[l]==0)
			continue;
		for(int c=0;c<image.colorbinnum;c++)
		{
			if(segcolorhist[l][c]!=0)
				bits+=-(double)segcolorhist[l][c]*log(segcolorhist[l][c]/(double)seghist[l])/log(2.0);
		}
	}
	//printf("bits = %f\n", bits);
	//printf("smoth + bits = %f\n",smoothenergy+bits);
	double costs = 0;
	if(labelcost!=0)
	{
		for(int l=0;l<maxnumSeg;l++)
		{
			if(seghist[l]!=0)
				costs +=labelcost;
		}
	}
	if((1191944==(int)(smoothenergy+bits*w_bits+costs))
		||(1286300==(int)(smoothenergy+bits*w_bits+costs)))
	{
		outv(smoothenergy);
		outv(bits*w_bits);
		outv(costs);
	}
	return smoothenergy+bits*w_bits+costs;
}

#endif