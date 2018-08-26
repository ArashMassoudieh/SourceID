// Distribution.cpp: implementation of the CDistribution class.
//
//////////////////////////////////////////////////////////////////////

#include "Distribution.h"
#include "DistributionNUnif.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CDistribution::CDistribution()
{
	n = 1;
	e = new double[n];
	s = new double[n];
	
}

CDistribution::CDistribution(int nn)
{
	n = nn;
	e = new double[n];
	s = new double[n];
	

}

CDistribution::~CDistribution()
{
	delete[] e;
	delete[] s;
}

CDistribution::CDistribution(const CDistribution &C)
{
	n = C.n;
	e = new double[n];
	s = new double[n];
	for (int i=0; i<n; i++)
	{
		e[i] = C.e[i];
		s[i] = C.s[i];
	}


}

CDistribution CDistribution::operator = (const CDistribution &C)
{
	n = C.n;
	e = new double[n];
	s = new double[n];
	for (int i=0; i<n; i++)
	{
		e[i] = C.e[i];
		s[i] = C.s[i];
	}

	return *this;

}

int CDistribution::GetRand()
{
	double x = GetRndUniF(0,1);
	int ii = 0;
	for (int i=0; i<n-1; i++)
	{	
		if (x<e[i] && x>s[i])
			ii = i;
	}
	return ii;

}