// Vector.h: interface for the CVector class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_VECTOR_H__27C038B4_05C9_4727_B9D2_0274D7801E3D__INCLUDED_)
#define AFX_VECTOR_H__27C038B4_05C9_4727_B9D2_0274D7801E3D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <iostream>
#include <vector>

using namespace std;

class CMatrix;
class SizeDist;
class CVector  
{
private:
	
	
public:
	vector<double> vec;
	CVector();
	CVector(int);
	CVector(const vector<double>, int);
	CVector(const vector<double> &v);
	CVector(const double x, int n);
	CVector(const CVector&);
	double& CVector::operator[](int);
	virtual ~CVector();
	int num;
	int CVector::range(int);
	CVector& CVector::operator=(const CVector&);
	CVector& CVector::operator=(const vector<double>&);
	CVector& CVector::operator=(const double &v);

	CVector& CVector::operator+();
	void CVector::swap(int , int );
	int CVector::getsize();
	CVector& CVector::operator*=(double);
	CVector& CVector::operator/=(double);
	CVector& CVector::operator+=(const CVector&);
	CVector& CVector::operator-=(const CVector&);
	CVector& CVector::operator*=(const CVector&);
	friend double dotproduct(CVector, CVector);
	friend CVector mult(CMatrix&, CVector&);
	friend CVector mult(CVector&, CMatrix&);
	friend double norm(CVector);
	friend double dotproduct(CVector v1, CVector v2);
	bool CVector::operator==(double v);
	double CVector::max();
	double CVector::min();
	double CVector::norm2();
	double CVector::sum();
	double CVector::abs_max();
	CMatrix CVector::T();
	CVector CVector::Log();
	void CVector::writetofile(FILE *f);
	CVector CVector::Exp();
	CMatrix CVector::diagmat();
};

CVector Log(CVector &V);
CVector Exp(CVector &V);
CVector operator+(CVector, CVector);
CVector operator+(double, CVector);
CVector operator+(CVector, double);
CVector operator-(CVector, CVector);
CVector operator-(double, CVector&);
CVector operator*(CVector, CVector);
CVector operator*(double, CVector);
CVector operator/(CVector, double); 
CVector operator/(CVector, CVector);
CVector operator/(double, CVector);
CVector zeros(int i);

#endif // !defined(AFX_VECTOR_H__27C038B4_05C9_4727_B9D2_0274D7801E3D__INCLUDED_)
