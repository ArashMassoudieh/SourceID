#pragma once
#include "Vector.h"
#include "Matrix.h"
#include "NormalDist.h"


class CNormalDist
{
public:
	CNormalDist(void);
	~CNormalDist(void);
	double CNormalDist::unitrandom();
	double CNormalDist::getstdnormalrand();
	double CNormalDist::getnormalrand(double mu, double std);
	CVector CNormalDist::getnormal(CVector &mu, CMatrix &sigma);
	CMatrix CNormalDist::getnormal(int m, int n, double mu, double std);
	CVector CNormalDist::getnormal(int m, double mu, double std);
	double CNormalDist::getlognormalrand(double mu, double std);
	CVector CNormalDist::getlognormal(CVector &mu, CMatrix &sigma);
	CMatrix CNormalDist::getlognormal(int m, int n, double mu, double std);
	CVector CNormalDist::getlognormal(int m, double mu, double std);
	

};

double getpdflognormal(CVector &X, CVector &mu, CMatrix &std);
double getpdfnormal(CVector &X, CVector &mu, CMatrix &std);