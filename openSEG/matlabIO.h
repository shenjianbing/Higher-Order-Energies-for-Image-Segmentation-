// Input and Output with matlab
#ifndef _MATLABIO_H__
#define _MATLABIO_H__

// Interface with Matlab
#ifdef _CHAR16T
#define CHAR16_T
#endif

#include <mat.h>
#pragma comment(lib,"libmat.lib")
#pragma comment(lib,"libmx.lib")


template <class T>
bool readmatintotable(Table2D<T> & table,char * filename,char * variablename);
template <class T>
bool savetalbeasmat(Table2D<T> & table,const char * filename,char * variablename);

template <class T>
bool readmatintotable(Table2D<T> & table,char * filename,char * variablename)
{
	MATFile *pmat; // pointer to matfile
	pmat=matOpen(filename, "r");//打开文件，返回指向文件指针
	if (pmat == NULL)
	{
		cout<<"Error opening file:" <<filename<<endl;
		return false;
	}
	const char **dir;
	int ndir; //ndir 表示mat文件中含有矩阵数目
	dir = (const char **)matGetDir(pmat, &ndir);
	if(dir == NULL)
	{
		printf("Error reading directory of file:");
		return false;
	}
	else
	{
		printf("Directory of %s\n",filename);
		for(int i=0; i < ndir; i++)
			printf("%s\n",dir[i]);//输出所含矩阵数目
	}
	mxArray *pa;
	pa = matGetVariable(pmat,variablename);
	int r = mxGetM(pa);
	int c = mxGetN(pa);
	table.resize(c,r);
	cout<<r<<' '<<c<<endl;
 
	double * dMat = new double[r*c];
	dMat = (double*)mxGetData(pa);

	table.resize(c,r);
	for(int i=0;i<r;i++)
	{
		for(int j=0;j<c;j++)
		{
			table[j][i] = dMat[j*r+i];
		}
	}
	
	matClose(pmat);
	cout<<"mat closed"<<endl;
	return true;
}

template <class T>
bool savetalbeasmat(Table2D<T> & table,const char * filename,char * variablename)
{
	int height = table.getHeight();
	int width = table.getWidth();
	MATFile * poutmatFile = matOpen(filename,"w");  
	mxArray *pouta = mxCreateDoubleMatrix(height,width,mxREAL);
	if (pouta == NULL) {
        printf("Unable to create mxArray with mxCreateDoubleMatrix\n");
        return false;
    }
    double *outfij = mxGetPr(pouta);
	for(int j=0;j<height;j++)
	{
		for(int i=0;i<width;i++)
		{
			outfij[i*height+j] = table[i][j];
		}
	}  
    matPutVariable(poutmatFile, variablename, pouta); 
	matClose(poutmatFile);
	return true;
}

#endif