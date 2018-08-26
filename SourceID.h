#pragma once
#include "vector.h"
#include "Source.h"
#include <iostream>
#include "matrix.h"
#include <vector>

using namespace std;

class CSourceID
{
public:
	CSourceID(void);
	CSourceID(int num_Sources, int num_Constts, int num_Samples);
	CSourceID(int num_Sources, int num_Constts, int num_Samples,int num_iso);
	int numSources, numConstts, numSamples, numisotopes;
	vector<CSource> obsconsts;
	vector<double> obsconst_mean;
	vector<double> obs_iso_conts_mean;
	double obsconst_std;
	double obs_iso_std;
	string path;
	vector<CSource> sourceconsts;
	CVector iniXm;
	double purtfact;
	vector<bool> inc;
	vector<int> majorisotopes;
	vector<double> sourcefracs;
	void CSourceID::getfromfile(string filename);
	void CSourceID::getfromfile_isotope(string filename);
	void CSourceID::getfromfile_v2(char filename[]);
	double CSourceID::calcerr();
	double CSourceID::sumfracs();
	double CSourceID::maxcnst(int i);
	CSourceID& CSourceID::operator = (const CSourceID&);
	CSourceID(const CSourceID&);
	void CSourceID::generate_sample();
	double CSourceID::fluctuate();
	double CSourceID::fluctuate(double x);
	double CSourceID::calcPostXY(int n);
	double CSourceID::calcPostXYMCMC(int n);
	double CSourceID::calcPostXYMCMC(int n, int interval, int nchains);
	CMatrix CSourceID::createY();
	CMatrix CSourceID::createY(double p);
	CMatrix CSourceID::createY_iso(double p);
	CVector CSourceID::createX();
	CVector CSourceID::createX(double p);
	double CSourceID::calcp1(CVector &X, CMatrix &Y, double obsstd);
	double CSourceID::calcp1_iso(CVector &X, CMatrix &Y, CMatrix &Y1, double obsstd_iso);
	double CSourceID::calcp2(CMatrix &Y);
	double CSourceID::calcp2_iso(CMatrix &Y);
	void CSourceID::GetIniYPrior();
	//void CSourceID::GetIniYPrior_v2();
	CVector CSourceID::getX(CVector Xl, double p);
	vector<CVector> Coutl, Coutt, Coutl_iso, Coutt_iso;
	vector<int> elements;
	
public:
	~CSourceID(void);
	
};
double normalrand(double m, double std);
//vector<string> getline(ifstream& file);
//vector<string> split(string s, char del);
string trim(string s);

