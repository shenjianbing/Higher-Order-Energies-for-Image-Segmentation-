#ifndef _PPBCQP_H__
#define _PPBCQP_H__
// class for quadratic programming
#include "PPBCBase.h"
#include "SparseMatrix.h"

enum PPBCMODE {PPBCB=0, PPBCL=1, PPBCT=2}; 

class PPBCQP: public PPBCBase{
public:
	PPBCQP()
	{
	}
	PPBCQP(const SparseMatrix<double> & M_, 
		const SparseMatrix<double> & U_, double C_, PPBCMODE mode_=PPBCT)
		:M(M_),U(U_),C(C_),mode(mode_){}
	virtual double computeenergy(const Table2D<Label> & labeling);
	virtual void updatemodel(); // collect statistics about current labeling
	virtual void * parabasegraph(double para_, UnknownRegion * unknownregion_p,
		double * flowoffset_p, vector<Point> * node_corr_p);
	virtual double getnewpara(ParaInterval paraInterval);
	double getlambdaproduct(const Table2D<Label> & labeling);
	void setmode(PPBCMODE mode_){mode=mode_;};
	PPBCMODE mode;
private:
	SparseMatrix<double> M;
	SparseMatrix<double> U; // unary term
	double C; // const
	
	SparseMatrix<double> M_neg; // submodular pairwise
	SparseMatrix<double> M_pos; // nonsubmodular pairwise
	Table2D<double> lambda_multipliers; // + lambda * (x_i-x_i^0)
};

double PPBCQP::computeenergy(const Table2D<Label> & labeling)
{
	double en = 0;
	for(int i=0;i<M.getsize();i++)
	{
		Trituple<double> elem = M[i];
		Point p1 = mappoint(elem.col);
		Point p2 = mappoint(elem.row);
		en += elem.val*(labeling[p1]==OBJ)*(labeling[p2]==OBJ);
	}
	for(int i=0;i<U.getsize();i++)
	{
		Trituple<double> elem = U[i];
		Point p = mappoint(elem.col);
		en += elem.val*(labeling[p]==OBJ);
	}
	en += C;
	return en;
}
void PPBCQP::updatemodel()
{
	M_neg = SparseMatrix<double>(M.width,M.height);
	M_pos = SparseMatrix<double>(M.width,M.height);
	for(int i=0;i<M.getsize();i++)
	{
		Trituple<double> elem = M[i];
		if(elem.val>0)
			M_pos.add(elem);
		else
			M_neg.add(elem);
	}

	if(PPBCB==mode){
		lambda_multipliers = Table2D<double>(img_w,img_h,1);
		return;
	}
	else
		lambda_multipliers = Table2D<double>(img_w,img_h,0);
	for(int i=0;i<M_pos.getsize();i++) // nonsubmodluar pairwise term
	{
		Trituple<double> elem = M_pos[i];
		assert(elem.val>1e-10);
		Point p1 = mappoint(elem.col);
		Point p2 = mappoint(elem.row);
		bool x = (initlabeling[p1]==OBJ);
		bool y = (initlabeling[p2]==OBJ);
		double val = elem.val;
		if(PPBCL==mode){
			lambda_multipliers[p1] +=(-x + 0.5) * val;
			lambda_multipliers[p2] +=(-y + 0.5) * val;
			continue;
		}
		// for PPBCT
		if(x==0&&y==0)
		{
			lambda_multipliers[p1] +=1;
			lambda_multipliers[p2] +=1;
		}
		else if(x==1&&y==1)
		{
			lambda_multipliers[p1] +=1;
			lambda_multipliers[p2] +=1;
		}
		else if(x==0&&y==1)
		{
			lambda_multipliers[p1] +=1;
		}
		else if(x==1&&y==0)
		{
			lambda_multipliers[p2] +=1;
		}
	}
}
void * PPBCQP::parabasegraph(double para_, UnknownRegion * unknownregion_p,
	double * flowoffset_p, vector<Point> * node_corr_p)
{
	Table2D<double> capsource(img_w,img_h,0);
	Table2D<double> capsink(img_w,img_h,0);
	if(PPBCB==mode) capsink = Table2D<double>(img_w,img_h,para_);
	SparseMatrix<double> pair_arcs(M.width,M.height);
	for(int i=0;i<U.getsize();i++) // unary term
	{
		Trituple<double> elem = U[i];
		capsink[mappoint(elem.col)] += elem.val;
	}
	for(int i=0;i<M_neg.getsize();i++) // submodluar pairwise term
	{
		Trituple<double> elem = M_neg[i];
		Point p1 = mappoint(elem.col);
		Point p2 = mappoint(elem.row);
		assert(elem.val<1e-5);
		capsource[p1] += -elem.val;
		pair_arcs.add(Trituple<double>(elem.col,elem.row,-elem.val));
	}
	for(int i=0;i<M_pos.getsize();i++) // nonsubmodluar pairwise term
	{
		Trituple<double> elem = M_pos[i];
		assert(elem.val>1e-10);
		Point p1 = mappoint(elem.col);
		Point p2 = mappoint(elem.row);
		bool x = (initlabeling[p1]==OBJ);
		bool y = (initlabeling[p2]==OBJ);
		double val = elem.val;
		if(PPBCL==mode){  // Laplacian parametric representation
			capsink[p1] +=(-para_*x + y + para_ / 2.0) * val;
			capsink[p2] +=(x- para_*y + para_ / 2.0) * val;
		}
		else if((PPBCT==mode)||(PPBCB==mode)){
			if(x==0&&y==0){
				capsink[p1] += val/2.0+para_*(PPBCT==mode);
				capsink[p2] += val/2.0+para_*(PPBCT==mode);
			}
			else if(x==1&&y==1){
				capsink[p1] += val/2.0+para_*(PPBCT==mode);
				capsink[p2] += val/2.0+para_*(PPBCT==mode);
			}
			else if(x==0&&y==1){
				capsink[p1] += val+para_*(PPBCT==mode);
			}
			else if(x==1&&y==0){
				capsink[p2] += val+para_*(PPBCT==mode);
			}
		}
	}
	if(unknownregion_p==NULL)
	{
		GraphType * g = new GraphType(img_w*img_h,img_w*img_h*2+M_neg.getsize());
		g->add_node(img_w*img_h);
		for(int i=0;i<img_w;i++)
		{
			for(int j=0;j<img_h;j++)
			{
				if((capsource[i][j]!=0)||(capsink[i][j]!=0))
					g->add_tweights(i+j*img_w,capsource[i][j],capsink[i][j]);
			}
		}
		for(int i=0;i<pair_arcs.getsize();i++)
		{
			Trituple<double> arc = pair_arcs[i];
			g->add_edge(arc.col,arc.row,arc.val,0);
		}
		return g;
	}
	else // monotonic
	{
		int unknownsize = unknownregion_p->unknownsize;
		double flowoffset = 0;
		GraphType * g = new GraphType(unknownsize,unknownsize*2+M_neg.getsize());
		g->add_node(unknownsize);

		Table2D<int> img_corr(img_w,img_h,-1);
		node_corr_p->resize(unknownsize,Point());
		int idx = 0;
		for(int i=0;i<img_w;i++)
		{
			for(int j=0;j<img_h;j++)
				if(unknownregion_p->knownlabeling[i][j]==UNKNOWN)
				{
					img_corr[i][j] = idx;
					(*node_corr_p)[idx++] = Point(i,j);
				}
		}

		for(int i=0;i<img_w;i++)
		{
			for(int j=0;j<img_h;j++)
			{
				if((capsource[i][j]!=0)||(capsink[i][j]!=0))
				{
					if(unknownregion_p->knownlabeling[i][j]==OBJ)
						flowoffset += capsink[i][j];
					else if(unknownregion_p->knownlabeling[i][j]==BKG)
						flowoffset += capsource[i][j];
					else if(unknownregion_p->knownlabeling[i][j]==UNKNOWN)
						g->add_tweights(img_corr[i][j],capsource[i][j],capsink[i][j]);
				}
			}
		}

		for(int i=0;i<pair_arcs.getsize();i++)
		{
			Trituple<double> arc = pair_arcs[i];
			Point p1 = mappoint(arc.col);
			Label l1 = unknownregion_p->knownlabeling[p1];
			Point p2 = mappoint(arc.row);
			Label l2 = unknownregion_p->knownlabeling[p2];
			if((l1==UNKNOWN)&&(l2==UNKNOWN))
				g->add_edge(img_corr[p1],img_corr[p2],arc.val,0);
			else if((l1!=UNKNOWN)&&(l2!=UNKNOWN))
				flowoffset += arc.val*(l1==OBJ)*(l2==BKG);
			else if((l1==UNKNOWN)&&(l2==BKG))
				g->add_tweights(img_corr[p1],0,arc.val);
			else if((l1==OBJ)&&(l2==UNKNOWN))
				g->add_tweights(img_corr[p2],arc.val,0);
		}
		*flowoffset_p = flowoffset;
		return g;
	}
}

double PPBCQP::getlambdaproduct(const Table2D<Label> & labeling)
{
	double product = 0;
	for(int i=0;i<img_w;i++)
	{
		for(int j=0;j<img_h;j++)
		{
			if(labeling[i][j]==OBJ)
				product += lambda_multipliers[i][j];
		}
	}
	return product;
}

double PPBCQP::getnewpara(ParaInterval interval)
{
	double a_low = getlambdaproduct(interval.bplow.solution);
	double b_low = interval.bplow.flow-a_low*interval.bplow.para;
	double a_up = getlambdaproduct(interval.bpup.solution);
	double b_up = interval.bpup.flow-a_up*interval.bpup.para;
	double newpara = (b_up - b_low) /(a_low - a_up);
	return newpara;
}

#endif