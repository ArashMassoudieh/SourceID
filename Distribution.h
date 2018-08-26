// Distribution.h: interface for the CDistribution class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_DISTRIBUTION_H__D12484DB_60B5_462D_949A_51A2B88BF5B7__INCLUDED_)
#define AFX_DISTRIBUTION_H__D12484DB_60B5_462D_949A_51A2B88BF5B7__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class CDistribution  
{
public:
	CDistribution();
	virtual ~CDistribution();
	int n;
	double *s;
	double *e;
	CDistribution::CDistribution(int nn);
	CDistribution(const CDistribution &C);
	CDistribution CDistribution::operator = (const CDistribution &C);
	int CDistribution::GetRand();

};

#endif // !defined(AFX_DISTRIBUTION_H__D12484DB_60B5_462D_949A_51A2B88BF5B7__INCLUDED_)
