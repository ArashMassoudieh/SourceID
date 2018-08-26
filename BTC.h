
#pragma once

#include <string>
#include <vector>

using namespace std;


class CBTC  
{
public:
	CBTC();
	CBTC(int n);
	virtual ~CBTC();
	int n;
	vector<double> t;
	vector<double> C;
	double CBTC::interpol(double x);
	CBTC(const CBTC &C);
	CBTC::CBTC(string Filename);
	//CBTC CBTC::operator = (const CBTC &C);
	CBTC& CBTC::operator = (const CBTC &C);
	void CBTC::readfile(string);
	//void CBTC::writefile(CString Filename);
	void CBTC::writefile(string Filename);
	//void CBTC::readfile(CString);
	double CBTC::tm();
	double CBTC::Vm();
	double CBTC::std();
	double CBTC::GetS0(CBTC &M);
	double CBTC::EMC(CBTC &M);
	double CBTC::maxC();
	double CBTC::Calculate_load(CBTC &M);
	double CBTC::GetS0(CBTC &M, CBTC &Q);
	void CBTC::setnumpoints(int);
	CBTC CBTC::Log();
	CBTC CBTC::Log(double min);
	double CBTC::sum2();
	double CBTC::percentile(double x);
	double CBTC::percentile(double x, int limit);
};

double diff(CBTC &BTC_p, CBTC &BTC_d);
double diff(CBTC BTC_p, CBTC BTC_d, int scale);
double diff(CBTC BTC_p, CBTC BTC_d, CBTC Q);
double diff2(CBTC BTC_p, CBTC BTC_d);
double R2(CBTC BTC_p, CBTC BTC_d);
CBTC operator*(double, CBTC);
CBTC operator*(CBTC, double);
double XYbar(CBTC BTC_p, CBTC BTC_d);
double X2bar(CBTC BTC_p, CBTC BTC_d);
double Y2bar(CBTC BTC_p, CBTC BTC_d);
double Ybar(CBTC BTC_p, CBTC BTC_d);
double Xbar(CBTC BTC_p, CBTC BTC_d);
double diff_n_norm(CBTC &BTC_p, CBTC &BTC_d, double std);
vector<double> QSort(vector<double> V);
