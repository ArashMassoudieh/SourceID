
#include "SourceID.h"
#include "Source.h"
#include "math.h"
#include "NormalDist.h"
#include "vector.h"
#include "matrix.h"
#include "DistributionNunif.h"
#include "NormalDist.h"
#include <string>
#include <omp.h>
#include <fstream>
#include "StringOP.h"

using namespace std;

CSourceID::CSourceID(void)
{
	
}

CSourceID::CSourceID(int num_Sources, int num_Constts, int num_Samples)
{
	
	numSources = num_Sources;
	numConstts = num_Constts;
	numSamples = num_Samples;
	numisotopes = 0;
	inc.resize(numSources);
	obsconsts.resize(numSamples);
	for (int i=0; i<numSamples; i++)
		obsconsts[i].SetnumConstts(numConstts);
	obsconst_mean.resize(numConstts);
	sourcefracs.resize(numSources);
	sourceconsts.resize(numSources);
	for (int i=0; i<numSources; i++)
		sourceconsts[i].SetnumConstts(numConstts);
}

CSourceID::CSourceID(int num_Sources, int num_Constts, int num_Samples,int num_iso)
{
	numSources = num_Sources;
	numConstts = num_Constts;
	numSamples = num_Samples;
	numisotopes = num_iso;
	inc.resize(numSources);
	obsconsts.resize(numSamples);
	
	for (int i=0; i<numSamples; i++)
		obsconsts[i].SetnumConstts(numConstts,numisotopes);
	obsconst_mean.resize(numConstts);
	sourcefracs.resize(numSources);
	sourceconsts.resize(numSources);
	for (int i=0; i<numSources; i++)
		sourceconsts[i].SetnumConstts(numConstts,numisotopes);

}

CSourceID& CSourceID::operator = (const CSourceID &C)
{
	numSources = C.numSources;
	numConstts = C.numConstts ;
	numSamples = C.numSamples;
	numisotopes = C.numisotopes;
	inc = C.inc;
	obsconsts = C.obsconsts;
	obsconst_mean = C.obsconst_mean;
	obsconst_std = C.obsconst_std;
	sourcefracs = C.sourcefracs;
	sourceconsts = C.sourceconsts;
	obs_iso_conts_mean = C.obs_iso_conts_mean;
	obs_iso_std = C.obs_iso_std;
	majorisotopes = C.majorisotopes;
	return *this;
}

CSourceID::CSourceID(const CSourceID &C)
{
	numSources = C.numSources;
	numConstts = C.numConstts ;
	numSamples = C.numSamples;
	inc = C.inc;
	obsconsts = C.obsconsts;
	obsconst_mean = C.obsconst_mean;
	obsconst_std = C.obsconst_std;
	sourcefracs = C.sourcefracs;
	sourceconsts = C.sourceconsts;
	obs_iso_conts_mean = C.obs_iso_conts_mean;
	obs_iso_std = C.obs_iso_std;
	majorisotopes = C.majorisotopes;
}
	
void CSourceID::getfromfile(string filename)
{
	FILE *fil;
//	fil = fopen(filename.c_str(), "r");
	ifstream file(filename);
	if (!file.good())
		return;

	string line;
//	while (file.eof() == false)
//	{
//		getline(file, line);
//		text.push_back(QString::fromStdString(line));
//	}
//	file.close();
//	return true;



	double ttt;
	for (int i=0; i<numConstts; i++)
	{
		string line;
		getline(file, line);  //read a line and split it by ","
		vector<string> s = split(line, ',');

		int c = 0;
		for (int j=0; j<numSamples; j++)
		{
			//fscanf(fil, "%le, ", &obsconsts[j].constts[i]);
			ttt = atof(s[c++].c_str());
			obsconsts[j].constts[i] = ttt;
		}

		for (int j=0; j<numSources; j++)
		{
			ttt = atof(s[c++].c_str());
//			fscanf(fil, "%le, ", &ttt);
			sourceconsts[j].constts_mean[i] = ttt;
			if (sourceconsts[j].constts_mean[i] == 0) sourceconsts[j].constts_mean[i] = 1E-8;
		}
//		fscanf(fil, "\n");
	}
	
	for (int i=0; i<numConstts; i++)
	{
//		vector<string> s = getline(file);  //read a line and split it by ","
		getline(file, line);  //read a line and split it by ","
		vector<string> s = split(line, ',');
		int c = 0;
		for (int j=0; j<numSamples; j++)
		{
			c++;
			//fscanf(fil, "%le, ", &ttt);
		}

		for (int j=0; j<numSources; j++)
		{
			ttt = atof(s[c++].c_str());
			//fscanf(fil, "%le, ", &ttt);
			sourceconsts[j].constts_var[i] = ttt;
		}
//		fscanf(fil, "\n");
	}
}


void CSourceID::getfromfile_isotope(string filename)
{
	FILE *fil;
	fil = fopen(filename.c_str(), "r");
	double ttt;
	for (int i=0; i<numisotopes; i++)
	{
		for (int j=0; j<numSamples; j++)
		{
			fscanf(fil, "%le, ", &obsconsts[j].iso[i]);
		}

		for (int j=0; j<numSources; j++)
		{
			fscanf(fil, "%le, ", &ttt);
			sourceconsts[j].iso_mean[i] = ttt;
			if (sourceconsts[j].constts_mean[i] == 0) sourceconsts[j].constts_mean[i] = 1E-8;
		}
		fscanf(fil, "\n");
	}
	
	for (int i=0; i<numisotopes; i++)
	{
		for (int j=0; j<numSamples; j++)
		{
			fscanf(fil, "%le, ", &ttt);
		}

		for (int j=0; j<numSources; j++)
		{
			fscanf(fil, "%le, ", &ttt);
			sourceconsts[j].iso_var[i] = ttt;
		}
		fscanf(fil, "\n");
	}
}


double CSourceID::calcerr()
{
	double err = 0;
	
	for (int k=0; k<numSamples; k++) 
	{
		for (int i=0; i<numConstts; i++)
		{
			double errcnst = 0;
			for (int j=0; j<numSources; j++)
				if (inc[j] == true) errcnst += sourcefracs[j]/sumfracs()*exp(sourceconsts[j].constts[i])/(1.0+exp(sourceconsts[j].constts[i]));
		
			err += pow(log(errcnst/(1.0-errcnst))-log(obsconsts[k].constts[i]/(1-obsconsts[k].constts[i])),2)/(2*obsconst_std*obsconst_std) + log(obsconst_std) + log(obsconsts[k].constts[i]) + log(1-obsconsts[k].constts[i]);
		}

		for (int i=0; i<numisotopes; i++)
		{	double num = 0;
			double denom = 0;
			for (int j=0; j<numSources; j++)
			{	num += sourcefracs[j]/sumfracs()*sourceconsts[j].iso[i]*exp(sourceconsts[j].constts[majorisotopes[i]])/(1+exp(sourceconsts[j].constts[majorisotopes[i]]));
				denom += sourcefracs[j]/sumfracs()*exp(sourceconsts[j].constts[majorisotopes[i]])/(1+exp(sourceconsts[j].constts[majorisotopes[i]]));
			}
			err += pow(log(num/denom)-log(obsconsts[k].iso[i]),2)/(2*obs_iso_std*obs_iso_std) + log(obs_iso_std) + log(obsconsts[k].iso[i]);

		}

	}
	return err;
}

CSourceID::~CSourceID(void)
{
	
}

double CSourceID::sumfracs()
{
	double sum = 0;
	for (int i=0; i<numSources; i++) sum+=sourcefracs[i];
	return sum;
}

//double CSourceID::maxcnst(int i)
//{
//	double maxcnsst=0;
//	for (int j=0; j<numSources; j++)
//		if (maxcnsst < sourceconsts[j].constts[i]) maxcnsst=sourceconsts[j].constts[i];
//	if (maxcnsst < obsconsts[i]) maxcnsst=obsconsts[i];
//	
//	return maxcnsst;
//
//}


double CSourceID::calcPostXYMCMC(int n, int interval, int nchains)
{
	double sumP=0;
	vector<CMatrix> Yl;
	vector<CMatrix> Yl_iso;
	vector<CVector> Xl;
	vector<CMatrix> Yt;
	vector<CMatrix> Yt_iso;
	vector<CVector> Xt;

	Yl.resize(nchains);
	Xl.resize(nchains);
	Yt.resize(nchains);
	Xt.resize(nchains);
	Yl_iso.resize(nchains);
	Yt_iso.resize(nchains);
	

	vector<double> obsconst_stdl;
	vector<double> obsconst_stdt;
	vector<double> obsiso_stdl;
	vector<double> obsiso_stdt;
	
	obsconst_stdl.resize(nchains);
	obsconst_stdt.resize(nchains);
	obsiso_stdl.resize(nchains);
	obsiso_stdt.resize(nchains);
	double unity = 1;
	Coutt.resize(nchains);Coutl.resize(nchains);
	Coutt_iso.resize(nchains);Coutl_iso.resize(nchains);
	CNormalDist ND;
	FILE *f;
	FILE *fnormal;
	FILE *fx_pop;
	FILE *fv_pop;
	FILE *fc_pop;
	FILE *fc_pop_iso;
	char str[3];
	FILE **fYoutput;
	FILE **fconstoutput;
	FILE **fY_iso_output;
	FILE **fisooutput;
	fYoutput = new FILE* [numConstts];
	fY_iso_output = new FILE* [numConstts];
	fconstoutput = new FILE* [elements.size()];
	fisooutput = new FILE* [numisotopes];
	for (int ii=0; ii<numConstts; ii++)
	{
		itoa (ii, str, 10 );
		string fYname = path + string("outputY_") + string(str) + string(".txt");
		fYoutput[ii] = fopen(fYname.c_str() , "w");
	}

	for (int ii=0; ii<numisotopes; ii++)
	{
		itoa (ii, str, 10 );
		string fYname = path + string("outputY_iso_") + string(str) + string(".txt");
		fY_iso_output[ii] = fopen(fYname.c_str() , "w");
	}

	for (int ii=0; ii<elements.size(); ii++)
	{
		itoa (elements[ii], str, 10 );
		string fYname = path + string("const_cont_") + string(str) + string(".txt");
		fconstoutput[ii] = fopen(fYname.c_str() , "w");
	}

	f = fopen((path + "output_MCMC.txt").c_str(), "w");
	fnormal = fopen((path+"output_normal.txt").c_str(), "w");
	fx_pop = fopen((path+"output_X.txt").c_str(), "w");
	fc_pop = fopen((path+"output_C.txt").c_str(), "w");
	fc_pop_iso = fopen((path+"output_iso_C.txt").c_str(), "w");
	fv_pop = fopen((path+"output_V.txt").c_str(), "w");
	
	vector<double> lnpt3, lnpt1, lnpt2, lnpl3, lnpl1, lnpl2, lnpl, lnpt;
	lnpt1.resize(nchains);lnpt2.resize(nchains);lnpt3.resize(nchains);lnpt.resize(nchains); lnpl.resize(nchains);
	lnpl1.resize(nchains);lnpl2.resize(nchains);lnpl3.resize(nchains);
	vector<double> crrctcoeff1;
	crrctcoeff1.resize(nchains);
	
	
	for (int i=0; i<nchains; i++)
	{	Yl[i] = createY(0.1);
		Yl_iso[i] = createY_iso(0.1);
		Xl[i] = createX(0.1);
		while (Xl[i].min()<=0)
			Xl[i] = createX(0.1);
		obsconst_stdl[i] = obsconst_std*(1 + ND.getnormalrand(0,0.05));
		obsiso_stdl[i] = obs_iso_std*(1 + ND.getnormalrand(0,0.05));

		if ((Xl[i][numSources-1]>0) && (Xl[i].min()>=0))
			crrctcoeff1[i] = 1;
		else
			crrctcoeff1[i] = 1e-200;

		lnpl1[i] = calcp1(Xl[i], Yl[i], obsconst_stdl[i]);// +calcp1_iso(Xl[i], Yl_iso[i], Yl[i], obsiso_stdl[i]);
		lnpl2[i] = calcp2(Yl[i]);// +calcp2_iso(Yl_iso[i]);
		lnpl3[i] = Log(Yl[i]).sum()+Log(1.0-Yl[i]).sum()+ Log(Xl[i]).sum();//+Log(Yl_iso[i]).sum()
	
		lnpl[i]=lnpl1[i]+lnpl2[i]+lnpl3[i]+log(crrctcoeff1[i]);
		Coutl[i] = Yl[i]*Xl[i];
		Coutl_iso[i] = Yl_iso[i]*Xl[i];
	}
	CVector count_stbl(nchains);
omp_set_num_threads(nchains);
cout<<"Number of Threads: "<<omp_get_num_threads()<<endl;
#pragma omp parallel for 
	for (int i=0; i<n; i++)
	{
		int k=i%nchains;
		Xt[k] = getX(Xl[k],purtfact/3);  
		CMatrix M = Log(Yl[k])-Log(1.0-Yl[k])+ND.getnormal(numConstts , numSources, 0.0, purtfact/10);
		Yt[k] = oneoneprod(Exp(M),1.0/(1.0+Exp(M)));
		Yt_iso[k] = CMatrix(numisotopes,numSources);
		for (int ii=0; ii<numSources; ii++)
			for (int jj=0; jj<numisotopes; jj++)
				Yt_iso[k][jj][ii] = Yl_iso[k][jj][ii]*exp(ND.getnormalrand(0.0, purtfact*sourceconsts[ii].iso_var[jj]));
		obsconst_stdt[k] = obsconst_stdl[k] + ND.getnormalrand(0,purtfact/3);
		obsiso_stdt[k] = obsiso_stdl[k] + ND.getnormalrand(0,purtfact*obs_iso_std/3);
		
		lnpt3[k] = Log(Yt[k]).sum()+Log(1.0-Yt[k]).sum()+Log(Xt[k]).sum();//+Log(Yt_iso[k]).sum()
		lnpt1[k] = calcp1(Xt[k], Yt[k], obsconst_stdt[k]);// +calcp1_iso(Xt[k], Yt_iso[k], Yt[k], obsiso_stdt[k]);
		lnpt2[k] = calcp2(Yt[k]);// +calcp2_iso(Yt_iso[k]);
		Coutt[k] = Yt[k]*Xt[k];		
		Coutt_iso[k] = Yt_iso[k]*Xt[k];	

		if ((Xt[k][numSources-1]>0) && (Xt[k].min()>=0) && (obsconst_stdl[k]>0))
			crrctcoeff1[k] = 1;
		else
			crrctcoeff1[k] = 1e-200;

		lnpt[k]=lnpt1[k]+lnpt2[k]+lnpt3[k]+log(crrctcoeff1[k]);
		double a1 = ND.unitrandom();
		double a2 = exp(lnpt[k]-lnpl[k]);
		if (a1<a2)
		{
			count_stbl[k] = 0;
			Xl[k] = Xt[k];
			Yl[k] = Yt[k];
			Yl_iso[k] = Yt_iso[k];
			lnpl[k] = lnpt[k];
			lnpl1[k] = lnpt1[k];
			lnpl2[k] = lnpt2[k];
			lnpl3[k] = lnpt3[k];
			obsconst_stdl[k] = obsconst_stdt[k];
			obsiso_stdl[k] = obsiso_stdt[k];
			Coutl[k] = Coutt[k];
			Coutl_iso[k] = Coutt_iso[k];
		}
		else
			count_stbl[k]++;
		
		if(i%10000 == 0)
			cout<<i<<"     "<<purtfact<<"     "<<count_stbl[k]<<endl;
		
		#pragma omp critical
		{
				
			if (int(i/interval)*interval==i)
			{	
				Coutl[k].writetofile(fc_pop);
				Coutl_iso[k].writetofile(fc_pop_iso);
				Xl[k].writetofile(fx_pop);
				fprintf(fv_pop,"%i, %i, %e, , %le, %e, %e, %e\n",i, k, obsconst_stdl[k], obsiso_stdl[k], lnpl[k], count_stbl[k], purtfact);
				for (int ii=0; ii<numConstts; ii++)
				{
					for (int jj=0;jj<numSources;jj++)
						fprintf(fYoutput[ii],"%e,",Yl[k][ii][jj]);
					fprintf(fYoutput[ii],"\n");
				}

				for (int ii=0; ii<numisotopes; ii++)
				{
					for (int jj=0;jj<numSources;jj++)
						fprintf(fY_iso_output[ii],"%e,",Yl_iso[k][ii][jj]);
					fprintf(fY_iso_output[ii],"\n");
				}

				for (int ii=0; ii<elements.size(); ii++)
				{
					CVector ConstCont = Xl[k]*Yl[k][ii]/(Xl[k]*Yl[k][ii]).sum(); 
					ConstCont.writetofile(fconstoutput[ii]);
				}
			}
				
		}
		
		if (((count_stbl.min()>20) || (count_stbl.max()>250)) && (purtfact>1e-6))
		{	purtfact*=0.75;
			count_stbl = 0;
		}
		
	}
	
	for (int ii=0; ii<elements.size(); ii++)
		fclose(fconstoutput[ii]);

	for (int jj=0;jj<numConstts;jj++)
		fclose(fYoutput[jj]);
	for (int jj=0;jj<numisotopes;jj++)
		fclose(fY_iso_output[jj]);
	fclose(fv_pop);
	fclose(f);
	fclose(fnormal);
	fclose(fx_pop);
	fclose(fc_pop); fclose(fc_pop_iso);
	return sumP;
}



CMatrix CSourceID::createY()
{
	CNormalDist NR;
	CMatrix h(numConstts, numSources);
	double a;
	for (int i=0; i<numSources; i++)
		for (int j=0; j<numConstts; j++)
		{	a = exp(sourceconsts[i].constts_mean[j]);
			h[j][i] = a/(1+a);
			if (h[j][i]<1e-18) j--;
		}

	return h;
}

CMatrix CSourceID::createY(double p)
{
	CNormalDist NR;
	CMatrix h(numConstts, numSources);
	double a;
	for (int i=0; i<numSources; i++)
		for (int j=0; j<numConstts; j++)
		{	a = exp(sourceconsts[i].constts_mean[j]+NR.getnormalrand(0,p));
			h[j][i] = a/(1+a);
			if (h[j][i]<1e-18) j--;
		}

	return h;
}

CMatrix CSourceID::createY_iso(double p)
{
	CNormalDist NR;
	CMatrix h(numisotopes, numSources);
	double a;
	for (int i=0; i<numSources; i++)
		for (int j=0; j<numisotopes; j++)
		{	h[j][i] = sourceconsts[i].iso_mean[j]*exp(NR.getnormalrand(0,p)*sourceconsts[i].iso_var[j]);
			if (h[j][i]<1e-18) j--;
		}

	return h;
}


CVector CSourceID::createX()
{
	CNormalDist NR;
	CVector V = CVector(sourcefracs);
	return V;
}

CVector CSourceID::createX(double p)
{
	CNormalDist NR;
	CVector V = CVector(sourcefracs)+NR.getnormal(numSources,0,p);
	double sum=0;
	for (int i=0; i<numSources-1; i++) sum+=V[i];
	V[numSources-1] = 1-sum;
	return V;
}

CVector CSourceID::getX(CVector Xl, double p)
{
	CNormalDist NR;
	CVector V1 = Exp(Log(Xl/(1.0-Xl))+NR.getnormal(numSources,0,p));
	CVector V = V1 / (1 + V1);
	double sum=0;
	for (int i=0; i<numSources-1; i++) sum+=V[i];
	V[numSources-1] = 1-sum;
	return V;
}

//CVector CSourceID::createX(double pp)
//{
//	CVector V = Exp(*X_m);
//
//	return V;
//}


double CSourceID::calcp1(CVector &X, CMatrix &Y, double obsstd)
{
	CVector Cout = Y*X/X.sum();
	double p = 0;
	for (int j=0; j<numSamples ; j++)
		for (int i=0; i<numConstts ; i++)
		{
			p += (-pow(log(obsconsts[j].constts[i]/(1-obsconsts[j].constts[i]))-log(Cout[i]/(1-Cout[i])),2)/(2*pow(obsstd,2))) - log(obsstd) - log(obsconsts[j].constts[i])-log(1-obsconsts[j].constts[i]);   
		}
	
	return p;

}

double CSourceID::calcp1_iso(CVector &X, CMatrix &Y, CMatrix &Y1, double obsstd_iso)
{
	
	CVector Cout(numisotopes);
	for (int i=0; i<numisotopes; i++)
	{
		double sumnumer = 0;
		double sumdeno = 0;
		for (int j=0; j<numSources; j++)
		{	sumnumer += Y[i][j]*Y1[majorisotopes[i]][j]*X[j];
			sumdeno += Y1[majorisotopes[i]][j]*X[j];
		}
		Cout[i] = sumnumer/sumdeno;
	}
	double p = 0;
	for (int j=0; j<numSamples ; j++)
		for (int i=0; i<numisotopes ; i++)
		{
			p += (-pow(log(obsconsts[j].iso[i])-log(Cout[i]),2)/(2*pow(obsstd_iso,2))) - log(obsstd_iso) - log(obsconsts[j].iso[i]);   
		}
	
	return p;

}
double CSourceID::calcp2(CMatrix &Y)
{
	
	double p = 0;
	for (int j=0; j<numSources ; j++)
		for (int i=0; i<numConstts ; i++)
		{
			p += -pow(sourceconsts[j].constts[i] -log(Y[i][j]/(1.0-Y[i][j])),2)/(2*pow(sourceconsts[j].constts_var[i],2)) - log(Y[i][j]) - log(1-Y[i][j]);   
		}
	
	return p;

}

double CSourceID::calcp2_iso(CMatrix &Y)
{
	double p = 0;
	for (int j=0; j<numSources ; j++)
		for (int i=0; i<numisotopes ; i++)
		{
			p += -pow(log(sourceconsts[j].iso[i]) -log(Y[i][j]),2)/(2*pow(sourceconsts[j].iso_var[i],2)) - log(Y[i][j]);   
		}
	
	return p;

}


void CSourceID::generate_sample()
{
	
	for (int i=0; i<numSources; i++)
	{	for (int j=0; j<numConstts; j++)
			//sourceconsts[i].constts[j] = normalrand(sourceconsts[i].constts_mean[j],sourceconsts[i].constts_var[j]);
			sourceconsts[i].constts[j] = sourceconsts[i].constts_mean[j];
		for (int j=0; j<numisotopes; j++)
			sourceconsts[i].iso[j] = sourceconsts[i].iso_mean[j];
	}
}	

double CSourceID::fluctuate()
{
	for (int j=0; j<numSources; j++)
		sourcefracs[j] = normalrand(sourcefracs[j], 0.1*log10(sourcefracs[j]+1E-16));
	return calcerr();	
}

double CSourceID::fluctuate(double x)
{
	for (int j=0; j<numSources; j++)
		sourcefracs[j] = normalrand(sourcefracs[j], x*0.1);
	return calcerr();	
}

double normalrand(double m, double std)
{
	double pi = atan(1.0)*4;
	double a1 = GetRndUniF(0,1);
	double a2 = GetRndUniF(0,1);

	double xx = sqrt(-2 * log(a1)) * cos(2 * pi * a2);
	
	return m*pow(10,xx*std);
}


//vector<string> getline(ifstream& file)
//{
//	string line;
//	while (file.good())
//	{
//		getline(file, line);
//		return split(line,',');
//	}
//}
//
//vector<string> split(string s, char del)
//{
//	int lastdel=0;
//	int j=0;
//	vector<string> strings;
//	for (int i=0; i<s.size(); i++)
//	{
//		if (s[i]==del)
//		{
//			strings.push_back(s.substr(lastdel, i-lastdel));
//			lastdel = i+1;
//			j++;
//		}
//	}
//	strings.push_back(trim(s.substr(lastdel, s.size()-lastdel)));
//	return strings;
//
//}
//
//string trim(string s)
//{
//	  return s.substr( s.find_first_not_of(' '), s.find_last_not_of(' ') + 1 );
//}

