// GA.h: interface for the CGA class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_GA_H__6D71FAC1_C9F2_400A_A150_475CB1945E5E__INCLUDED_)
#define AFX_GA_H__6D71FAC1_C9F2_400A_A150_475CB1945E5E__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "individual.h"
#include <stdio.h>
#include "Simplex.h"
#include "Distribution.h"
#include "SourceID.h"
#include <vector>

using namespace std;

class CGA  
{
public:
	CGA();
	CGA(int n);
	CGA(int n, int nParams);
	CGA(const CGA &C);
	CGA CGA::operator=(CGA &C);
	virtual ~CGA();
	int maxpop;
	vector<CIndividual> Ind;
	vector<CIndividual> Ind_old;
	void CGA::initialize();
	void CGA::Setminmax(int a, double minrange, double maxrange, int prec);
	bool fixedinputs[15];
	double fixedinputvale[15];
	double sumfitness;
	void CGA::fitnessdistini();
	CDistribution fitdist;
	void CGA::crossover();
	double CGA::avgfitness();
	void CGA::mutate(double mu);
	void CGA::assignfitnesses();
	int CGA::maxfitness();
	double CGA::variancefitness();
	double N;
	double pcross;
	double CGA::stdfitness();
	double CGA::avg_actual_fitness();
	char fitnesstype;
	double exponentcoeff;
	int CGA::optimize(int nGens);
	int cross_over_type;
	void CGA::setnumpop(int n);
	double pmute;
	double MaxFitness;
	double CGA::avg_inv_actual_fitness();
	int CGA::optimize(int nGens, string RunFileName);

	void CGA::setnparams(int n);
	double pertscale, nonpurt;

	void CGA::assignrank();
	void CGA::assignfitness_rank(double N);
	double biasfact;
	CSourceID *SID;
	double shakescale, shakescalered;
	void CGA::shake();
	void CGA::enhance(CSIMPLEX &CSimp, CIndividual &InD);
	void CGA::enhance(int t_rank, CSIMPLEX &CSimp);
	int numenhancements, num_enh;
	bool RCGA;
	void CGA::crossoverRC();
	void CGA::fillfitdist();
};

#endif // !defined(AFX_GA_H__6D71FAC1_C9F2_400A_A150_475CB1945E5E__INCLUDED_)
