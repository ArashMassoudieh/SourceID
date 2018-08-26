#include "BTCSet.h"
#include "string.h"
#include <fstream>
#include <string>
#include "DistributionNUnif.h"
#include "stdafx.h"

using namespace std;

CBTCSet::CBTCSet(void)
{
	nvars = 0;
	BTC.resize(nvars);
	solved = true;
}

CBTCSet::~CBTCSet(void)
{

}

CBTCSet::CBTCSet(int n)
{
	nvars = n;
	BTC.resize(nvars);
	solved = true;
}

void CBTCSet::writetofile(char outputfile[])
{
	FILE *Fil;
	Fil = fopen(outputfile, "w");

	for (int j=0; j<maxnumpoints(); j++)
	{
		for (int i=0; i<nvars; i++)
		{
			if (i<BTC[i].n)
				fprintf(Fil, "%lf, %le,", BTC[i].t[j], BTC[i].C[j]);
			else
				fprintf(Fil, ", ,");
		}
		fprintf(Fil, "\n");
	}
	fclose(Fil);
}

void CBTCSet::writetofile(string outputfile)
{
	FILE *Fil;
	Fil = fopen(outputfile.c_str() , "w");

	for (int j=0; j<maxnumpoints(); j++)
	{
		for (int i=0; i<nvars; i++)
		{
			if (i<BTC[i].n)
				fprintf(Fil, "%lf, %le,", BTC[i].t[j], BTC[i].C[j]);
			else
				fprintf(Fil, ", ,");

		}
		fprintf(Fil, "\n");
	}

	fclose(Fil);

}

int CBTCSet::maxnumpoints()
{
	int m = 0;
	for (int i=0; i<nvars; i++)
		if (BTC[i].n>m) m = BTC[i].n;
	
	return m;
}

CBTCSet::CBTCSet(const CBTCSet &B)
{
	nvars = B.nvars;
	BTC.resize(nvars);

	for (int i=0; i<nvars; i++)
		BTC[i] = B.BTC[i]; 

}

CBTCSet& CBTCSet::operator = (const CBTCSet &B)
{
	nvars = B.nvars;
	BTC.resize(nvars);

	for (int i=0; i<nvars; i++)
		BTC[i] = B.BTC[i]; 
	
	return *this;

}

void CBTCSet::getfromfile(const string filename)
{
	ifstream BTCfile(filename);
	vector<string> s;
	int i=0;
	while (BTCfile.eof()==false)
	{
		s = getline(BTCfile);
		if (i==0)
		{
			nvars = s.size()-1;
			BTC.resize(nvars);
		}
		for (int j=0; j<nvars; j++)
		{	BTC[j].C.push_back(atof(s[j+1].c_str()));
			BTC[j].t.push_back(atof(s[0].c_str()));
			BTC[j].n++;
		}
		i++;
	}
	BTCfile.close();

}

CBTCSet::CBTCSet(const string filename)
{
	ifstream BTCfile(filename);
	vector<string> s;
	int i=0;
	if (BTCfile.good()) 
	while (BTCfile.eof()==false)
	{
		s = getline(BTCfile);
		if (i==0)
		{
			nvars = s.size()-1;
			BTC.resize(nvars);
		}
		for (int j=0; j<nvars; j++)
		{	if (s.size()-1==nvars)
			{	BTC[j].C.push_back(atof(s[j].c_str()));
				BTC[j].t.push_back(j);
				BTC[j].n++;
			}
		}
		i++;
	}
	BTCfile.close();

}

vector<string> getline(ifstream& file)
{
	string line;
	while (file.good())
	{
		getline(file, line);
		return split(line,',');
	}
}

vector<string> split(string s, char del)
{
	int lastdel=0;
	int j=0;
	vector<string> strings;
	char space = ' ';
	for (int i=0; i<s.size(); i++)
	{
		if (s[i]==del)
		{
			strings.push_back(s.substr(lastdel, i-lastdel));
			if (s[i+1] != space) lastdel = i+1; else lastdel = i+2;
			j++;
		}
	}
	strings.push_back(s.substr(lastdel, s.size()-lastdel));
	return strings;

}

vector<double> CBTCSet::getrandom()
{
	int a = int(GetRndUniF(0,BTC[0].n));
	vector<double> res(nvars);
	for (int i=0; i<nvars; i++)
		res[i] = BTC[i].C[a];

	return res;
}

vector<double> CBTCSet::percentile(double x)
{
	vector<double> v;
	for (int i=0; i<nvars; i++)
		v.push_back(BTC[i].percentile(x));

	return v;
}

vector<double> CBTCSet::percentile(double x, int limit)
{
	vector<double> v;
	for (int i=0; i<nvars; i++)
		v.push_back(BTC[i].percentile(x,limit));

	return v;
}
