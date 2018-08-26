// SpareMatrix.h: interface for the CSpareMatrix class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_SPAREMATRIX_H__B31AD1FB_C6FC_44D1_8D7E_96824F117C46__INCLUDED_)
#define AFX_SPAREMATRIX_H__B31AD1FB_C6FC_44D1_8D7E_96824F117C46__INCLUDED_


#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class CVector;
class CMatrix;
class CSpareMatrix  
{
public:
	CSpareMatrix();
	CSpareMatrix(int n_n, int n_nnz);
	CSpareMatrix(int n_n);
	CSpareMatrix(CSpareMatrix&);
	virtual ~CSpareMatrix();
	int *rowind;
	double *value;
	int *colsub;
	int nnz;
	int n;
	void CSpareMatrix::setval(int ii,int ji, double v);
	double CSpareMatrix::getval(int ii,int ji);
	CVector CSpareMatrix::getcol(int i);
	CVector CSpareMatrix::getrow(int i);

	int CSpareMatrix::getnumrows();
	int CSpareMatrix::getnumcols();
	CSpareMatrix& CSpareMatrix::operator=(const CSpareMatrix&);
	CSpareMatrix& CSpareMatrix::operator+=(const CSpareMatrix&);
	friend CSpareMatrix mult(CSpareMatrix, CSpareMatrix);
	friend CVector diag(CSpareMatrix);
	int CSpareMatrix::colass(int i);
	CSpareMatrix CSpareMatrix::Transpose();
	void CSpareMatrix::addtoval(int i,int j, double val);
	double CSpareMatrix::sumrows(int i);
	double CSpareMatrix::sumcols(int i);
	friend CMatrix inverse(CSpareMatrix A);
	
};

CSpareMatrix operator+(CSpareMatrix, CSpareMatrix);
CSpareMatrix operator*(CSpareMatrix, CSpareMatrix);
CVector operator*(CSpareMatrix, CVector);
CVector operator/(CVector, CSpareMatrix);
CSpareMatrix Transpose(CSpareMatrix M1);
CVector SpareSolve(CSpareMatrix, CVector);
CSpareMatrix mult(CSpareMatrix A, CSpareMatrix B);
CVector mult(CSpareMatrix A, CVector B);
CVector RichardsonIteration(CSpareMatrix A, CVector B, double tol);
CVector BiCG(CSpareMatrix A, CVector B, double tol);



#endif // !defined(AFX_SPAREMATRIX_H__B31AD1FB_C6FC_44D1_8D7E_96824F117C46__INCLUDED_)
