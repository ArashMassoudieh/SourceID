// SIMPLEX.h: interface for the CSIMPLEX class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_SIMPLEX_H__0F086C8D_C14E_4BAD_91A5_4A653F35F4A1__INCLUDED_)
#define AFX_SIMPLEX_H__0F086C8D_C14E_4BAD_91A5_4A653F35F4A1__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "vector.h"
#include "Individual.h"
#include "SourceID.h"


class CSIMPLEX  
{
public:
	CSIMPLEX();
	virtual ~CSIMPLEX();
	CVector **p;
	int ndim;	
	CSIMPLEX(int n);
	CSIMPLEX(int n,CVector X,double fac);
	CVector *f;
	double CSIMPLEX::onestepsimplex(int i);
	void CSIMPLEX::nstepsimplex(int n, double fac);
	double funk(CVector &X);
	int CSIMPLEX::wst();
	int CSIMPLEX::bst();
	int CSIMPLEX::reflect(int ww, int bb);
	int CSIMPLEX::expand(int ww, int bb);
	int CSIMPLEX::contract(int ww, int bb);
	int CSIMPLEX::multcontract(int ww, int bb);
	CVector *pref;
	CVector *pext;
	CVector *pcon;
	CVector *pmcn;
	double fref,fext,fcon,fmcn;
	double fopt;
	CVector *popt;
	double tol;
	int w,b;
	bool *fixedinput;
	double *fixedvalues;
	CIndividual *Ind;
	CSourceID *SID;

};

#endif // !defined(AFX_SIMPLEX_H__0F086C8D_C14E_4BAD_91A5_4A653F35F4A1__INCLUDED_)

