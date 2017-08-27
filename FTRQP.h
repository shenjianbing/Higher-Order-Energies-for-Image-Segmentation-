#include <FTRBase.h>
// class for quadratic programming using FTR
#include "SparseMatrix.h"

class FTRQP: public FTRBase{
public:
	FTRQP(){}
	FTRQP(const SparseMatrix<double> & M_, 
		const SparseMatrix<double> & U_, double C_):M(M_),U(U_),C(C_)
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
	}
	virtual void computeenergy(const Table2D<Label> & labeling, 
		double * act_energy_p, double * app_energy_p);
	// collect statistics about current labeling and update approximate model
	virtual void updatemodel(); 
	virtual GraphType * ftrbasegraph(double lambda_);
	Point mappoint(int idx){return Point(idx%img_w,idx/img_w);}
private:
	SparseMatrix<double> M;
	SparseMatrix<double> U; // unary term
	double C; // const
	
	SparseMatrix<double> M_neg; // submodular pairwise
	SparseMatrix<double> M_pos; // nonsubmodular pairwise

};

void FTRQP::computeenergy(const Table2D<Label> & labeling, 
		double * act_energy_p, double * app_energy_p)
{
	double act_M = 0;
	for(int i=0;i<M.getsize();i++)
	{
		Trituple<double> elem = M[i];
		Point p1 = mappoint(elem.col);
		Point p2 = mappoint(elem.row);
		act_M += elem.val*(labeling[p1]==OBJ)*(labeling[p2]==OBJ);
	}
	double act_U = 0;
	for(int i=0;i<U.getsize();i++)
	{
		Trituple<double> elem = U[i];
		Point p = mappoint(elem.col);
		act_U += elem.val*(labeling[p]==OBJ);
	}
	*act_energy_p = act_M + act_U + C;
	if(app_energy_p)
	{
		double app_M = 0;
		for(int i=0;i<M_neg.getsize();i++)
		{
			Trituple<double> elem = M_neg[i];
			Point p1 = mappoint(elem.col);
			Point p2 = mappoint(elem.row);
			app_M += elem.val*(labeling[p1]==OBJ)*(labeling[p2]==OBJ);
		}
		for(int i=0;i<M_pos.getsize();i++)
		{
			Trituple<double> elem = M_pos[i];
			Point p1 = mappoint(elem.col);
			Point p2 = mappoint(elem.row);
			int x0=0, y0=0, x=0, y=0;
			if(current_labeling[p1]==OBJ) x0=1;
			if(current_labeling[p2]==OBJ) y0=1;
			if(labeling[p1]==OBJ) x=1;
			if(labeling[p2]==OBJ) y=1;
			app_M += (y0*(x-x0)+x0*(y-y0)+x0*y0)*elem.val;
		}
		*app_energy_p = app_M + act_U + C;
	}
}
void FTRQP::updatemodel()
{
	dt = getDistanceTransform(current_labeling);
}

GraphType * FTRQP::ftrbasegraph(double lambda_)
{
	GraphType * g = new GraphType(img_w*img_h,img_w*img_h*2+M_neg.getsize());
	g->add_node(img_w*img_h);
	Table2D<double> capsource(img_w,img_h,0);
	Table2D<double> capsink(img_w,img_h,0);
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
		assert(elem.val<1e-10);
		capsource[p1] += -elem.val;
		g->add_edge(elem.col,elem.row,-elem.val,0);
	}
	for(int i=0;i<M_pos.getsize();i++) // nonsubmodluar pairwise term
	{
		Trituple<double> elem = M_pos[i];
		assert(elem.val>1e-10);
		Point p1 = mappoint(elem.col);
		Point p2 = mappoint(elem.row);
		int x0=0, y0=0;
		if(current_labeling[p1]==OBJ) x0=1;
		if(current_labeling[p2]==OBJ) y0=1;
		double val = elem.val;

		capsink[p1] += y0 * val;
		capsink[p2] += x0 * val;
	}
	for(int i=0;i<img_w;i++)
	{
		for(int j=0;j<img_h;j++)
		{
			// hamming distance
			if(current_labeling[i][j]==OBJ)
				capsource[i][j] += lambda_;
			else if(current_labeling[i][j]==BKG)
				capsink[i][j] += lambda_;
			//if((capsource[i][j]!=0)||(capsink[i][j]!=0))
			g->add_tweights(i+j*img_w,capsource[i][j],capsink[i][j]);
			// L2 distance
			//g->add_tweights(i+j*img_w,capsource[i][j]+dt[i][j]*lambda_,capsink[i][j]);
		}
	}
	return g;
}