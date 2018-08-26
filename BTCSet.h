#pragma once
#include "BTC.h"
#include <vector>

class CBTCSet
{
public:
	CBTCSet(void);
	CBTCSet(int n);
	CBTCSet(const CBTCSet &BTC);
	int nvars;
	vector<CBTC> BTC;
	void CBTCSet::writetofile(char outputfile[]);
	int CBTCSet::maxnumpoints();
	CBTCSet& CBTCSet::operator = (const CBTCSet &C);
	void CBTCSet::writetofile(string outputfile);
	CBTCSet::CBTCSet(const string filename);
	void CBTCSet::getfromfile(const string filename);
	vector<double> CBTCSet::getrandom();
	bool solved;
	vector<double> CBTCSet::percentile(double x);
	vector<double> CBTCSet::percentile(double x, int limit);
public:
	~CBTCSet(void);
};

vector<string> getline(ifstream& file);
vector<string> split(string s, char del);
