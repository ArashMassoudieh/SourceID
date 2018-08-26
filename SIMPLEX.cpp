// SIMPLEX.cpp: implementation of the CSIMPLEX class.
//
//////////////////////////////////////////////////////////////////////

#include "SIMPLEX.h"
#include "math.h"


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CSIMPLEX::CSIMPLEX()
{

}

CSIMPLEX::~CSIMPLEX()
{
	delete f;
	delete pref;
	delete pext;
	delete pcon;
	delete pmcn;
	delete popt;
	delete SID;
}

CSIMPLEX::CSIMPLEX(int n)
{
	ndim = n;
	p = new CVector *[ndim+1];
	for (int i=0;i<ndim+1;i++)
		p[i] = new CVector(ndim);
	f = new CVector;
	*f = CVector(ndim+1);
	pref= new CVector;
	pext= new CVector;
	pcon= new CVector;
	pmcn= new CVector;
	Ind = new CIndividual;
	*pref= CVector(ndim);
	*popt= CVector(ndim);
	*pext= CVector(ndim);
	*pcon= CVector(ndim);
	*pmcn= CVector(ndim);
	tol=1E12;


}

CSIMPLEX::CSIMPLEX(int n,CVector X, double fac)
{
	ndim = n;
	p = new CVector *[ndim+1];
	for (int i=0;i<ndim+1;i++)
		p[i] = new CVector(ndim);

	*p[0] = X;

		
	for (int i=1; i<ndim; i++)
	{
		*p[i]=X;
		p[i]->vec[i-1] = p[i]->vec[i-1]*fac;
	}
	f = new CVector;
	*f = CVector(ndim+1);
	Ind = new CIndividual;
	pref= new CVector(ndim);
	pext= new CVector(ndim);
	pcon= new CVector(ndim);
	popt= new CVector(ndim);
	pmcn= new CVector(ndim);
	tol=1E12;
}



double CSIMPLEX::funk(CVector &X)
{
	int j=0;

	double sumfitness = 0;

	for (int i=0; i<SID->numSources; i++)
	{	SID->sourcefracs[i] = pow(10,X[i]);
		SID->inc[i] = true;
	}
		
	sumfitness = SID->calcerr();

	return sumfitness;
}


void CSIMPLEX::nstepsimplex(int n, double fac)
{
	CVector X(ndim);
	for (int i=0; i<ndim; i++)
		X[i] = (*Ind).x[i];

	*p[0] = X;
		
	for (int i=1; i<ndim+1; i++)
	{
		*p[i]=X;
		p[i]->vec[i-1] = p[i]->vec[i-1]*fac;
	}
	
	for (int j = 0; j<n; j++)
		onestepsimplex(j);

	b = bst();
	fopt= f->vec[b];
	*popt= (*p[b]);

}

double CSIMPLEX::onestepsimplex(int ii)
{


	if (ii==0)
	{
		for (int i=0; i<ndim+1;i++)
			f->vec[i]=funk(*p[i]);
	}
	else
		f->vec[w]=funk(*p[w]);
	
	w=wst();
	b=bst();

	tol = fabs(f->vec[w]-f->vec[b]);
	int k=reflect(w,b);
	int l;
	if (k==0)
		*p[w]=(*pref);
	if (k==1)
	{	l=expand(w,b);
		if (l==0) *p[w]=(*pext);
		if (l==1) *p[w]=(*pref);
	}
	else if (k==2)
	{	
		l=contract(w,b);
		if (l==0) *p[w]=(*pcon);
		if (l==1) 
			int l1 = multcontract(w,b);
	}

	return fopt;

}

int CSIMPLEX::wst()
{
	int w;
	double ws=-1E12;
	for (int i=0; i<ndim+1; i++)
	{
		double a = f->vec[i];
		if (a>ws)
		{
			w=i;
			ws=a;
		}
	}
	return w;
}

int CSIMPLEX::bst()
{
	int b;
	double bs=+1E12;
	for (int i=0; i<ndim+1; i++)
	{
		double a = f->vec[i];
		if (a<bs)
		{
			b=i;
			bs=a;
		}
	}
	return b;
}

int CSIMPLEX::reflect(int ww, int bb)
{
	CVector sump(ndim);
	for (int i=0; i<ndim+1; i++)
		if (i!=ww)
			sump+=*p[i];
	CVector Cent=sump/ndim;
	*pref=2.0*Cent - (*p[ww]);
	fref=funk(*pref);
	if (fref<f->vec[bb]) return 1;
	if (fref>f->vec[ww]) return 2;
	return 0;

}

int CSIMPLEX::expand(int ww, int bb)
{
	CVector sump(ndim);
	for (int i=0; i<ndim+1; i++)
		if (i!=ww)
			sump+=*p[i];
	CVector Cent=sump/ndim;
	*pext=3.0*Cent - 2.0*(*p[ww]);
	fext=funk(*pext);
	if (fext<fref)
		return 0;
	else
		return 1;
	
}

int CSIMPLEX::contract(int ww, int bb)
{
	CVector sump(ndim);
	for (int i=0; i<ndim+1; i++)
		if (i!=ww)
			sump+=*p[i];
	CVector Cent=sump/ndim;
	*pcon=3.0*Cent - 2.0*(*p[ww]);
	fcon=funk(*pcon);
	if (fcon<f->vec[ww])
		return 0;
	else
		return 1;
}

int CSIMPLEX::multcontract(int ww, int bb)
{
	for (int i=0; i<ndim+1; i++)
		if (i!=bb)
			*p[i] = 0.5*((*p[bb])+(*p[i]));
	return 1;
}

