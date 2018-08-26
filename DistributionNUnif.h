// Distribution.h: interface for the CDistribution class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_DISTRIBUTION_H__96ABBAC0_5401_4D37_A842_2FEC7C889B95__INCLUDED_)
#define AFX_DISTRIBUTION_H__96ABBAC0_5401_4D37_A842_2FEC7C889B95__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000



class CDistributionNUnif  
{
public:
	CDistributionNUnif();
	CDistributionNUnif(int n);
	virtual ~CDistributionNUnif();
	int n;
	double *x;
	double *y;
	double CDistributionNUnif::GetRndNorm(double mean, double std);
	double CDistributionNUnif::GetRndGamma();
	void CDistributionNUnif::initializeNormal(double dx0, double dxmult,int nint);
	void CDistributionNUnif::initializeGamma(double dx0, double dxmult, int nint, double r, double lambda);
	bool set;
	bool symetrical;
	CDistributionNUnif(const CDistributionNUnif &D);
	CDistributionNUnif CDistributionNUnif::operator=(const CDistributionNUnif &D);
	
};

double Gammapdf(double x, double r, double lambda);
double NormalStdpdf(double x);
double gamma(double x);
double GetRndUniF(double xmin, double xmax);

#endif // !defined(AFX_DISTRIBUTION_H__96ABBAC0_5401_4D37_A842_2FEC7C889B95__INCLUDED_)
