// TDMatrix.cpp: implementation of the CTDMatrix class.
//
//////////////////////////////////////////////////////////////////////

#include "TDMatrix.h"
#include <iostream>

using namespace std;
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CTDMatrix::CTDMatrix(int m)
{
	numrows = m;
	numcols = 3;
	mat = new CVector* [numrows];
	
	for (int i = 0;i<numrows; ++i)
	{	
		mat[i] = new CVector(numcols);
	}
}

CTDMatrix::CTDMatrix(CMatrix &m)
{
	numrows = m.getnumrows();
	numcols = 3;
	mat = new CVector* [numrows];
	
	for (int i=0;i<numrows; ++i) mat[i] = new CVector(numcols);
	
}


CTDMatrix::~CTDMatrix()
{
	for (int i = numrows; i>0; --i)
		delete mat[i-1];
	delete mat;
}


CVector& CTDMatrix::operator[](int i)
{
	return *mat[i];
}

int CTDMatrix::getnumrows() {return numrows;};

void TDMA(CTDMatrix& m,CVector& v, CVector& x)
{
	int j;
	double bet;
	int n = v.getsize();

	CVector gam(n);
	

	x[0] = v[0]/(bet=m[0][1]);
	for (j=1;j<=n-1;j++)	{
		gam[j] = m[j-1][2]/bet;
		bet = m[j][1] - m[j][0]*gam[j];
		x[j] = (v[j] - m[j][0]*x[j-1])/bet;
	}

	for (j=(n-2); j>=0; j--)
		x[j]-=gam[j+1]*x[j+1];
	
}


CVector operator/(CVector &V, CTDMatrix &M)
{
	int n = M.getnumrows();
	CVector b(n);
	TDMA(M, V, b);
	return b;
}

