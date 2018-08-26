#include "StdAfx.h"
#include "Vector.h"
#include "Matrix.h"
#include "NormalDist.h"
#include <cstdlib>
#include "math.h"

CNormalDist::CNormalDist(void)
{
	
}

CNormalDist::~CNormalDist(void)
{

}

double CNormalDist::unitrandom()
{
	double x1 = rand();
	double x2 = rand();
	return (x1+x2*(double(RAND_MAX)-1))/(double(RAND_MAX)*double(RAND_MAX));
}

double CNormalDist::getstdnormalrand()
{
	double x1 = unitrandom();
	double x2 = unitrandom();
	double pi = atan(1.0)*4;
	double y1 = sqrt(-2*log(x1))*cos(2*pi*x2);
	return y1;
}

double CNormalDist::getnormalrand(double mu, double std)
{
	return getstdnormalrand()*std+mu;
}

CVector CNormalDist::getnormal(CVector &mu, CMatrix &sigma)
{
	CMatrix L = sigma.Cholesky_factor();
	CVector V = getnormal(L.getnumrows(),0,1);
	return L*V + mu;
}

CMatrix CNormalDist::getnormal(int m, int n, double mu, double std)
{
	CMatrix M(m,n);
	for (int i=0; i<m; i++)
		for (int j=0; j<n; j++)
			M[i][j] = getstdnormalrand()*std+mu;

	return M;
}

CVector CNormalDist::getnormal(int m, double mu, double std)
{
	CVector M(m);
	for (int i=0; i<m; i++)
		M[i] = getstdnormalrand()*std+mu;

	return M;
}

double CNormalDist::getlognormalrand(double mu, double std)
{
	return exp(getstdnormalrand()*std+mu);
}

CVector CNormalDist::getlognormal(CVector &mu, CMatrix &sigma)
{
	return Exp(getnormal(mu,sigma));
}

CMatrix CNormalDist::getlognormal(int m, int n, double mu, double std)
{
	CMatrix M(m,n);
	for (int i=0; i<m; i++)
		for (int j=0; j<n; j++)
			M[i][j] = exp(getstdnormalrand()*std+mu);

	return M;
}

CVector CNormalDist::getlognormal(int m, double mu, double std)
{
	CVector M(m);
	for (int i=0; i<m; i++)
		M[i] = exp(getstdnormalrand()*std+mu);

	return M;
}

double getpdfnormal(CVector &X, CVector &mu, CMatrix &std)
{
	int k = X.num;
	double pi = atan(1.0)*4;
	CMatrix exparg = 0.5*(X.T()-mu.T())*Invert(std)*(X-mu);

	double pdf = 1/pow(2*pi,k)/sqrt(std.det())*exp(-exparg[0][0]);
	return pdf;

}

double getpdflognormal(CVector &X, CVector &mu, CMatrix &std)
{
	int k = X.num;
	double pi = atan(1.0)*4;
	CMatrix exparg = 0.5*(Log(X.T())-mu.T())*Invert(std)*(Log(X)-mu);

	double pdf = 1/pow(2*pi,k)/sqrt(fabs(std.det()))*exp(-exparg[0][0]);
	return pdf;

}