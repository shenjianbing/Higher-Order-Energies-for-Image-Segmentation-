#ifndef _QPPROGRAM_H__
#define _QPPROGRAM_H__
//Quadratic Programming
#include "PPBCQP.h"
#include "FTRQP.h"
#include <matlabIO.h>

void QPprogram(int argc, char * argv[]);

void QPprogram(int argc, char * argv[])
{
	char * rootdir = "D:/repositories/PBO_public_code/pairwise/";
	
	// read data file
	char * qpfilepath = argv[1];
	Table2D<double> c_table;
	Table2D<double> u_table;
	Table2D<double> m_table;
	readmatintotable(c_table,qpfilepath,"C");
	readmatintotable(u_table,qpfilepath,"U");
	readmatintotable(m_table,qpfilepath,"M");
	cout<<m_table[0][0]<<m_table[1][0]<<m_table[2][0]<<endl;
	cout<<m_table[0][1]<<m_table[1][1]<<m_table[2][1]<<endl;
	// C U M
	int img_w = atoi(argv[2]);
	int img_h = atoi(argv[3]);
	int img_size = img_w * img_h;
	SparseMatrix<double> M(img_size,img_size);
	for(unsigned int i=0;i<m_table.getHeight();i++)
	{
		int col_id = ((int)m_table[0][i]-1)/img_h;
		int row_id = ((int)m_table[0][i]-1)%img_h;
		assert(col_id>=0,"col id should be positive");
		assert(row_id>=0,"col id should be positive");
		int id1 = col_id + row_id * img_w;
		col_id = ((int)m_table[1][i]-1)/img_h;
		row_id = ((int)m_table[1][i]-1)%img_h;
		int id2 = col_id + row_id * img_w;
		M.add(Trituple<double>(id1,id2,m_table[2][i]));
	}
	cout<<M.getsize()<<endl;
	SparseMatrix<double> U(img_size,img_size);
	for(unsigned int i=0;i<u_table.getHeight();i++)
	{
		int col_id = i/img_h;
		int row_id = i%img_h;
		int id = col_id + row_id * img_w;
		U.add(Trituple<double>(id,id,u_table[0][i]));
	}
	cout<<u_table.getMax()<<endl;
	cout<<u_table.getMin()<<endl;
	cout<<"c "<<c_table[0][0]<<endl;
	cout<<"Quadratic Programming :)"<<endl;

	double best_e = INFTY;
	int itr_count=0;

	Table2D<Label> initlabeling = Table2D<Label>(img_w,img_h,OBJ);
	PPBCQP PPBCqp(M,U,c_table[0][0]);
	Table2D<Label> currlabeling;

	/*cout<<"Auxiliary Cut"<<endl;
	currlabeling = initlabeling;
	best_e = INFTY;
	itr_count=0;
	
	while(++itr_count)
	{
		cout<<itr_count<<" th iteration"<<endl;
		PPBCqp.setinitlabeling(currlabeling);
		BreakPoint bp = PPBCqp.explorepara(0);
		bp.original_e = PPBCqp.computeenergy(bp.solution);
		cout<<"breakpoint parameter: "<<bp.para<<endl;
		cout<<"breakpoint energy: "<<bp.original_e<<endl;
		if(bp.original_e<best_e-0.001)
		{
			currlabeling = bp.solution;
			best_e = bp.original_e;
		}
		else
			break;
	}
	savebinarylabelingBW(currlabeling,to_Cstr(rootdir<<"qp_auxcut.bmp"));
	cout<<"Auxiliary Cut energy: %"<<best_e<<"%"<<endl;
	return;*/

	cout<<"PBO"<<endl;
	PPBCqp.setmode(PPBCT);
	// parameters for deconvolution
	/*if(PPBCL == PPBCqp.mode)
		PPBCqp.setpara(0,2,0.01,1,true);
	else if((PPBCB == PPBCqp.mode)||(PPBCT == PPBCqp.mode))
		PPBCqp.setpara(-0.2,0.2,0.001,1,true);*/
	// parameters for stenosis
	if(PPBCB==PPBCqp.mode)
		PPBCqp.setpara(-40,40,0.1,1,true);
	else if(PPBCT==PPBCqp.mode)
		PPBCqp.setpara(-20,20,0.1,1,true);
	else if(PPBCL==PPBCqp.mode)
		PPBCqp.setpara(-20,20,0.1,1,true);
	currlabeling= initlabeling;
	best_e = INFTY;
	itr_count=0;
	while(++itr_count)
	{
		cout<<itr_count<<" th iteration"<<endl;
		PPBCqp.setinitlabeling(currlabeling);
		PPBCqp.explore();
		BreakPoint bp = PPBCqp.SelectBestBP();
		cout<<"breakpoint parameter: "<<bp.para<<endl;
		cout<<"breakpoint energy: "<<bp.original_e<<endl;
		if(bp.original_e<best_e-0.001)
		{
			currlabeling = bp.solution;
			best_e = bp.original_e;
		}
		else
			break;
	}
	savebinarylabelingBW(currlabeling,to_Cstr(rootdir<<"qp_pbo.bmp"));
	cout<<"PBO energy: %"<<best_e<<"%"<<endl;

	/*// FTR
	FTRQP ftr(M,U,c_table[0][0]);
	ftr.setpara(0.01, 1e-3, 2, 0.25, 200, 0.1);
	ftr.setinitlabeling(initlabeling);
	ftr.optimize();
	savebinarylabelingBW(ftr.current_labeling,to_Cstr(rootdir<<"qp_"<<"_ftr.bmp"));
	cout<<"FTR energy: %"<<ftr.current_e<<"%"<<endl;
	cout<<"Foreground size: sss"<<countintable(ftr.current_labeling,OBJ)<<"sss"<<endl;
	//return;*/
	return;
}

#endif