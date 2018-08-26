// BTC.cpp: implementation of the CBTC class.
//
//////////////////////////////////////////////////////////////////////

#include "BTC.h"
#include "math.h"
#include "string.h"
#include "Header.h"
#include "stdafx.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CBTC::CBTC()
{
	
	n=0;
	t.resize(n);
	C.resize(n);
	for (int i=0; i<n; i++) {C[i]=0; t[i]=0;}
}

CBTC::CBTC(int n1)
{
	n=n1;
	t.resize(n);
	C.resize(n);
	for (int i=0; i<n; i++) {C[i]=0; t[i]=0;}
}

void CBTC::setnumpoints(int n1)
{
	
	n = n1;
	t.resize(n);
	C.resize(n);
	for (int i=0; i<n; i++) {C[i]=0; t[i]=0;}


}

CBTC::~CBTC()
{
	
}

CBTC::CBTC(const CBTC &CC)
{
	n=CC.n;
	t.resize(n);
	C.resize(n);
	for (int i=0; i<n; i++)
	{
		t[i] = CC.t[i];
		C[i] = CC.C[i];
	}
}

CBTC::CBTC(string Filename)
{
	FILE *FILEBTC;
	FILEBTC = fopen(Filename.c_str(), "r");
	int numpoints = 0;
	double tt, CC;
	fscanf(FILEBTC, "%i\n", &numpoints);
		
	n=numpoints;
	t.resize(n);
	C.resize(n);

	
	for (int i=0; i<numpoints; i++)
	{	
		fscanf(FILEBTC, "%le, %le\n", &t[i], &C[i]);
	}
	fclose(FILEBTC);
}

/*CBTC CBTC::operator = (const CBTC &CC)
{
	n=CC.n;
	t = new double[n];
	C = new double[n];
	for (int i=0; i<n; i++)
	{
		t[i] = CC.t[i];
		C[i] = CC.C[i];
	}
	
	return *this;
}*/

CBTC& CBTC::operator = (const CBTC &CC)
{
	n=CC.n;
	t.resize(n);
	C.resize(n);
	for (int i=0; i<n; i++)
	{
		t[i] = CC.t[i];
		C[i] = CC.C[i];
	}
	
	return *this;
}

CBTC CBTC::Log()
{
	CBTC BTC = CBTC(n);
	for (int i=0; i<n; i++)
	{
		BTC.t[i] = t[i];
		BTC.C[i] = log(C[i]);
	}
	return BTC;
}

CBTC CBTC::Log(double m)
{
	CBTC BTC(n);
	for (int i=0; i<n; i++)
	{
		BTC.t[i] = t[i];
		BTC.C[i] = log(max(C[i],m));
	}
	return BTC;
}

double CBTC::interpol(double x)
{
	double r=0;
	for (int i=0; i<n-1; i++)
	{
		if (t[i] <= x && t[i+1] >= x)
			r=(C[i+1]-C[i])/(t[i+1]-t[i])*(x-t[i]) + C[i];
	}
	if (x>t[n-1]) r=C[n-1];
	if (x<t[0]) r=C[0];

	return r;

}



double diff(CBTC BTC_p, CBTC BTC_d, int scale)
{
	double sum = 0;
	for (int i=0; i<BTC_d.n; i++)
	{
		if (BTC_d.C[i] > BTC_p.interpol(BTC_d.t[i]))
			sum += scale*pow(BTC_d.C[i] - BTC_p.interpol(BTC_d.t[i]),2)/sqrt(1.0+scale*scale);
		else
			sum += pow(BTC_d.C[i] - BTC_p.interpol(BTC_d.t[i]),2)/sqrt(1.0+scale*scale);
	}
	return sum;
}

double diff(CBTC &BTC_p, CBTC &BTC_d)
{
	double sum = 0;
	double sumvar1 = 0;
	double a;
	for (int i=0; i<BTC_d.n; i++)
	{	
		a = BTC_p.interpol(BTC_d.t[i]);
		sum += pow(BTC_d.C[i] - a,2)/BTC_d.n;
		sumvar1 += BTC_d.C[i]*BTC_d.C[i]/BTC_d.n;
	}

	return sum/sumvar1;
}

double diff_n_norm(CBTC &BTC_p, CBTC &BTC_d, double std)
{
	double sum = 0;
	double a;
	for (int i=0; i<BTC_d.n; i++)
	{	
		a = BTC_p.interpol(BTC_d.t[i]);
		sum += pow(BTC_d.C[i] - a,2)/(2*std*std) + log(std);
		
	}
	return sum;
}

double diff2(CBTC BTC_p, CBTC BTC_d)
{
	double sum = 0;
	double sumvar1 = 0;
	for (int i=0; i<BTC_d.n; i++)
	{	
		sum += pow(BTC_d.C[i] - BTC_p.interpol(BTC_d.t[i]),2);
		sumvar1 += BTC_d.C[i]*BTC_d.C[i];
	}

	return sum;
}

double CBTC::sum2()
{
	
	double sumvar1 = 0;
	for (int i=0; i<n; i++)
	{	
		sumvar1 += C[i]*C[i];
	}

	return sumvar1;
}


double R2(CBTC BTC_p, CBTC BTC_d)
{
	double sumcov = 0;
	double sumvar1 = 0;
	double sumvar2 = 0;
	double sum1 = 0;
	double sum2 = 0;
	for (int i=0; i<BTC_d.n; i++)
	{
		double x2 = BTC_p.interpol(BTC_d.t[i]);
		sumcov += BTC_d.C[i]*x2/BTC_d.n;
		sumvar1 += BTC_d.C[i]*BTC_d.C[i]/BTC_d.n;
		sumvar2 += x2*x2/BTC_d.n;	
		sum1 += BTC_d.C[i]/BTC_d.n;
		sum2 += x2/BTC_d.n;
	}

	return pow(sumcov-sum1*sum2,2)/(sumvar1-sum1*sum1)/(sumvar2-sum2*sum2);
}

double XYbar(CBTC BTC_p, CBTC BTC_d)
{
	double sumcov = 0;
	double sumvar1 = 0;
	double sumvar2 = 0;
	double sum1 = 0;
	double sum2 = 0;
	for (int i=0; i<BTC_d.n; i++)
	{
		double x2 = BTC_p.interpol(BTC_d.t[i]);
		sumcov += BTC_d.C[i]*x2/BTC_d.n;
		sumvar1 += BTC_d.C[i]*BTC_d.C[i]/BTC_d.n;
		sumvar2 += x2*x2/BTC_d.n;	
		sum1 += BTC_d.C[i]/BTC_d.n;
		sum2 += x2/BTC_d.n;
	}

	return sumcov;
}

double X2bar(CBTC BTC_p, CBTC BTC_d)
{
	double sumcov = 0;
	double sumvar1 = 0;
	double sumvar2 = 0;
	double sum1 = 0;
	double sum2 = 0;
	for (int i=0; i<BTC_d.n; i++)
	{
		double x2 = BTC_p.interpol(BTC_d.t[i]);
		sumcov += BTC_d.C[i]*x2/BTC_d.n;
		sumvar1 += BTC_d.C[i]*BTC_d.C[i]/BTC_d.n;
		sumvar2 += x2*x2/BTC_d.n;	
		sum1 += BTC_d.C[i]/BTC_d.n;
		sum2 += x2/BTC_d.n;
	}

	return sumvar1;
}

double Y2bar(CBTC BTC_p, CBTC BTC_d)
{
	double sumcov = 0;
	double sumvar1 = 0;
	double sumvar2 = 0;
	double sum1 = 0;
	double sum2 = 0;
	for (int i=0; i<BTC_d.n; i++)
	{
		double x2 = BTC_p.interpol(BTC_d.t[i]);
		sumcov += BTC_d.C[i]*x2/BTC_d.n;
		sumvar1 += BTC_d.C[i]*BTC_d.C[i]/BTC_d.n;
		sumvar2 += x2*x2/BTC_d.n;	
		sum1 += BTC_d.C[i]/BTC_d.n;
		sum2 += x2/BTC_d.n;
	}

	return sumvar2;
}

double Ybar(CBTC BTC_p, CBTC BTC_d)
{
	double sumcov = 0;
	double sumvar1 = 0;
	double sumvar2 = 0;
	double sum1 = 0;
	double sum2 = 0;
	for (int i=0; i<BTC_d.n; i++)
	{
		double x2 = BTC_p.interpol(BTC_d.t[i]);
		sumcov += BTC_d.C[i]*x2/BTC_d.n;
		sumvar1 += BTC_d.C[i]*BTC_d.C[i]/BTC_d.n;
		sumvar2 += x2*x2/BTC_d.n;	
		sum1 += BTC_d.C[i]/BTC_d.n;
		sum2 += x2/BTC_d.n;
	}

	return sum2;
}

double Xbar(CBTC BTC_p, CBTC BTC_d)
{
	double sumcov = 0;
	double sumvar1 = 0;
	double sumvar2 = 0;
	double sum1 = 0;
	double sum2 = 0;
	for (int i=0; i<BTC_d.n; i++)
	{
		double x2 = BTC_p.interpol(BTC_d.t[i]);
		sumcov += BTC_d.C[i]*x2/BTC_d.n;
		sumvar1 += BTC_d.C[i]*BTC_d.C[i]/BTC_d.n;
		sumvar2 += x2*x2/BTC_d.n;	
		sum1 += BTC_d.C[i]/BTC_d.n;
		sum2 += x2/BTC_d.n;
	}

	return sum1;
}


double diff(CBTC BTC_p, CBTC BTC_d, CBTC Q)
{
	double sum = 0;
	for (int i=0; i<BTC_d.n; i++)
	{
		sum += pow(BTC_d.C[i] - BTC_p.interpol(BTC_d.t[i]),2)*pow(Q.interpol(BTC_d.t[i]),2);
	}
	return sum;
}

void CBTC::readfile(string Filename)
{
	FILE *FILEBTC;
	FILEBTC = fopen(Filename.c_str(), "r");
	int numpoints = 0;
	double tt, CC;
	fscanf(FILEBTC, "%i\n", &numpoints);
	/*while (feof(FILEBTC)==false)
	{	
		fscanf(FILEBTC, "%lf, %lf", &tt, &CC);
		numpoints++;
	}
	//numpoints;
	fclose(FILEBTC);*/
	
	n=numpoints;
	t.resize(n);
	C.resize(n);

	//FILEBTC = fopen(Filename.c_str(), "r");
	//fscanf(FILEBTC, "%lf\n", &tt);
	for (int i=0; i<numpoints; i++)
	{	
		fscanf(FILEBTC, "%lf, %lf", &t[i], &C[i]);
	}
	fclose(FILEBTC);

}

/*void CBTC::readfile(CString Filename)
{
	FILE *FILEBTC;
	FILEBTC = fopen(Filename, "r");
	if (FILEBTC == NULL) 
		double e=1;
	int numpoints = 0;
	double tt, CC;
	while (feof(FILEBTC)==false)
	{	
		fscanf(FILEBTC, "%lf, %lf\n", &tt, &CC);
		numpoints++;
	}
	//numpoints--;
	fclose(FILEBTC);
	
	n=numpoints;
	t = new double[numpoints];
	C = new double[numpoints];

	FILEBTC = fopen(Filename, "r");
	for (int i=0; i<numpoints; i++)
	{	
		fscanf(FILEBTC, "%lf, %lf", &t[i], &C[i]);
	}
	fclose(FILEBTC);



}*/

void CBTC::writefile(string Filename)
{
	FILE *FILEBTC;
	FILEBTC = fopen(Filename.c_str(), "w");
	for (int i=0; i<n; i++)
		fprintf(FILEBTC, "%lf, %le\n", t[i], C[i]);
	
	fclose(FILEBTC);

}

/*void CBTC::writefile(CString Filename)
{
	FILE *FILEBTC;
	FILEBTC = fopen(Filename, "w");
	for (int i=0; i<n; i++)
		fprintf(FILEBTC, "%lf, %le\n", t[i], C[i]);
	
	fclose(FILEBTC);

}*/

/*double CBTC::GetS0(CBTC &M)
{
	double sumprod = 0;
	double sumsqr = 0;
	for (int i = 0; i<M.n; i++)
	{
		sumprod += M.C[i]*interpol(M.t[i]);
		sumsqr += interpol(M.t[i])*interpol(M.t[i]);
	}
	double S0 = sumprod/sumsqr;
	return S0;
}

double CBTC::GetS0(CBTC &M, CBTC &Q)
{
	double sumprod = 0;
	double sumsqr = 0;
	for (int i = 0; i<M.n; i++)
	{
		sumprod += M.C[i]*interpol(M.t[i])*pow(Q.interpol(M.t[i]),2);
		sumsqr += interpol(M.t[i])*interpol(M.t[i])*pow(Q.interpol(M.t[i]),2);
	}
	double S0 = sumprod/sumsqr;
	return S0;
}*/

CBTC operator*(double alpha, CBTC CBTC_T)
{
	CBTC S(CBTC_T.n);
	for (int i=0; i<CBTC_T.n; i++)
		S.C[i] = alpha*CBTC_T.C[i];

	return S;
}

CBTC operator*(CBTC CBTC_T, double alpha)
{
	CBTC S = CBTC_T;
	for (int i=0; i<CBTC_T.n; i++)
		S.C[i] = alpha*CBTC_T.C[i];

	return S;
}

/*double CBTC::EMC(CBTC &M)
{
	double sum = 0;
	double sumflow = 0;
	for (int i=0; i<n; i++)
	{	
		sum += C[i]*M.interpol(t[i]);
		sumflow += M.interpol(t[i]);
	}
	if (sumflow == 0.0)
		return 0;
	else
		return sum/sumflow;
}

double CBTC::Calculate_load(CBTC &M)
{
	double sum = 0;
	double sumflow = 0;
	for (int i=0; i<n; i++)
	{	
		sum += C[i]*M.interpol(t[i])*(t[2]-t[1]);
		
	}
	
	return sum;
}*/

double CBTC::maxC()
{
	double max = 0;
	for (int i=0; i<n; i++)
	{	if (C[i]>max)
			max = C[i];
	}
	return max;
}

double CBTC::percentile(double x)
{
	vector<double> X = QSort(C);
	int i = int(x*X.size());
	return X[i];

}

double CBTC::percentile(double x, int limit)
{
	vector<double> C1(C.size()-limit);
	for (int i=0; i<C1.size(); i++)
		C1[i] = C[i+limit];
	vector<double> X = QSort(C1);
	int ii = int(x*double(X.size()));
	return X[ii];

}

vector<double> QSort(vector<double> V)
{
	if (V.size()<=1) return V;
	int end = V.size();
	vector<double> less, greater;
	greater.push_back(V[end-1]);
	for (int i=0; i<end-1; i++)
		if (V[i]<V[end-1]) less.push_back(V[i]);
		else greater.push_back(V[i]);
		
	
	vector<double> res = QSort(less);
	if ((V==greater) & (less.size()==0)) return greater;
	vector<double> x2 = QSort(greater);
	
	res.insert(res.end(), x2.begin(), x2.end());
	
	return res;

}

