#pragma once
#include <vector>

using namespace std;

class CSource
{
public:
	CSource(void);
	CSource(int nn);
	CSource(int nn, int nn_iso);
	CSource& CSource::operator=(const CSource&);
	CSource(const CSource&);
	void CSource::SetnumConstts(int nn);
	void CSource::SetnumConstts(int nn, int nn_iso);
	int n, n_iso;
	vector<double> constts;
	vector<double> constts_var;
	vector<double> constts_mean;
	vector<double> iso;
	vector<double> iso_var;
	vector<double> iso_mean;
public:
	~CSource(void);
};
