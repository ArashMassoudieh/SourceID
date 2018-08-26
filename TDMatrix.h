// TDMatrix.h: interface for the CTDMatrix class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_TDMATRIX_H__E79AF454_AE6F_4181_91C1_22F120D7BE1D__INCLUDED_)
#define AFX_TDMATRIX_H__E79AF454_AE6F_4181_91C1_22F120D7BE1D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#include "vector.h"
#include "matrix.h"

class CVector;
class CMatrix;
class CTDMatrix  
{
private:
	int numcols;
	int numrows;
	CVector **mat;
	int range(int);
public:
	CTDMatrix(int);
	CTDMatrix(CMatrix&);
	virtual ~CTDMatrix();
	CVector& operator[](int);
	CTDMatrix& CTDMatrix::operator=(const CTDMatrix&);
	int CTDMatrix::getnumrows();
};

void TDMA(CTDMatrix&,CVector&, CVector&);
CVector operator/(CVector&, CTDMatrix&);

#endif // !defined(AFX_TDMATRIX_H__E79AF454_AE6F_4181_91C1_22F120D7BE1D__INCLUDED_)
