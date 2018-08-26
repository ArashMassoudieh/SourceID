// Source_ID_MCMC.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "Matrix.h"
#include "Vector.h"
#include "NormalDist.h"
#include "Source.h"
#include "SourceID.h"
#include <iostream>
#include "GA.h"
#include "math.h"
#include <string>
#include <fstream>
#include "BTCSet.h"
#include "Header.h"
#include "StringOP.h"

using namespace std;


int _tmain(int argc, _TCHAR* argv[])
{
	string filename;
	cout<<"Input File Name: ";
	cin>>filename;
	
	ifstream infile(filename);
	if (!infile.good())
		return 0;
	int nrand, sampstart;
	int nparam, nconstts, npop, ngen, nsamples, npermutes, popmcmc, interval, nsample, n_chains, purt_fact, burnin, nisotopes;
	vector<int> majisotopes;
	bool forward;
	int forward1;
	string constfilename, isofilename;
	vector<int> elements;
	vector<double> percentiles;
	nsamples = 1;
	npermutes = 1;
	vector<string> ss;
	string path,GA_out_filename,Percent_out_filename, C_Percent_out_filename, iso_Percent_out_filename;
	double purturbfact, mutation, pcross, purt_red_rate;

	while (!infile.eof()) 
	{
		ss = getline(infile);
		if (ss[0]=="path") path=ss[1];
		if (ss[0]=="forward") forward = atoi(ss[1].c_str());
		if (ss[0]=="n_sources") nparam = atoi(ss[1].c_str());
		if (ss[0]=="n_generations") ngen = atoi(ss[1].c_str());
		if (ss[0]=="n_elems") nconstts = atoi(ss[1].c_str());
		if (ss[0]=="n_isotopes") nisotopes = atoi(ss[1].c_str());
		if (ss[0]=="maj_isotopes") 
			{for (int i=0; i<ss.size()-1; i++) majisotopes.push_back(atoi(ss[i+1].c_str()));} 
		if (ss[0]=="n_chains") n_chains = atoi(ss[1].c_str());
		if (ss[0]=="ga_pop") npop = atoi(ss[1].c_str());
		if (ss[0]=="n_sample") nsample = atoi(ss[1].c_str());
		if (ss[0]=="n_mcmcpop") {popmcmc = atoi(ss[1].c_str()); interval=atoi(ss[2].c_str());};
		if (ss[0]=="profile") constfilename = ss[1].c_str();
		if (ss[0]=="isotopes") isofilename = ss[1].c_str();
		if (ss[0]=="purturb") purturbfact = atof(ss[1].c_str());
		if (ss[0]=="mutation") mutation = atof(ss[1].c_str());
		if (ss[0]=="pcross") pcross = atof(ss[1].c_str());
		if (ss[0]=="purt_red_rate") purt_red_rate = atof(ss[1].c_str());
		if (ss[0]=="burn_in") burnin = atoi(ss[1].c_str());
		if (ss[0]=="elements") {elements.resize(ss.size()-1); for (int i=0; i<ss.size()-1; i++) elements[i] = atoi(ss[i+1].c_str());};
		if (ss[0]=="ga_output") GA_out_filename = ss[1];
		if (ss[0]=="percent_output") Percent_out_filename = ss[1];
		if (ss[0]=="percentiles") {percentiles.resize(ss.size()-1); for (int i=0; i<ss.size()-1; i++) percentiles[i] = atof(ss[i+1].c_str());}
		if (ss[0]=="c_percent_output") C_Percent_out_filename = ss[1];
		if (ss[0]=="iso_percent_output") iso_Percent_out_filename = ss[1];
	}



	CGA G(npop, nparam+2);
	G.N = 1;
	G.pcross = pcross;
	G.pmute = mutation;
	G.pertscale = purturbfact;
	G.shakescale = purturbfact;
	G.shakescalered = purt_red_rate;
	for (int i=0; i<nparam; i++)
		G.Setminmax(i,-3,3,3);
	G.Setminmax(nparam,0,6,4);
	G.Setminmax(nparam+1,0,6,4);
	G.initialize();
	G.RCGA = false;
	
	CSourceID CS(nparam,nconstts,nsample, nisotopes);
	CS.majorisotopes = majisotopes;
	if (forward == false) 
	{//	CS.getfromfile(constfilename);
		CS.getfromfile(path + constfilename);
		if (isofilename!="") CS.getfromfile_isotope(path + isofilename);
	}	
	else
//		CS.getfromfileX(constfilename);

	for (int i=0; i<CS.numSources; i++)
	{	CS.sourcefracs[i] = 1;
		CS.inc[i] = true;
	}
	
	*G.SID = CS;
	
	char outfilename [50];
	char ga_filename [50];
		
	//sprintf (ga_filename, (path + GA_out_filename).c_str());
	G.SID->generate_sample();
	FILE *foutput;
	/*sprintf (outfilename, (path + GA_out_filename).c_str());*/
	foutput = fopen((path + "output.txt").c_str(), "w");	

	if (forward == false) {
		G.optimize(ngen,path + GA_out_filename);
			
		double sumxx = 0;
		for (int i=0; i<nparam; i++)
			sumxx += max(pow(10, G.Ind[G.maxfitness()].x[i]),0.0001);

		for (int i=0; i<G.SID->numConstts; i++)
		{
			for (int j=0; j<G.SID->numSamples; j++)
				fprintf(foutput, "%le, ", G.SID->obsconsts[j].constts[i]);

			for (int j=0; j<G.SID->numSources; j++)
			{
				fprintf(foutput,"%le, ", G.SID->sourceconsts[j].constts[i]);
			}
			fprintf(foutput, "\n");
		}
		fprintf(foutput, "\n");

		for (int i=0; i<nparam; i++)
		{	
			G.SID->sourcefracs[i] = max(pow(10, G.Ind[G.maxfitness()].x[i]),0.0001)/sumxx;
			fprintf(foutput, "%le, ", pow(10, G.Ind[G.maxfitness()].x[i])/sumxx);
		}
		G.SID->obsconst_std = G.Ind[G.maxfitness()].x[nparam];
		fprintf(foutput, "\n");
		double reeval = G.SID->fluctuate(0);
		fprintf(foutput, "%le\n", reeval);
	}
	for (int i=0; i<G.SID->numConstts; i++)
	{
		double sumelem1 = 0;
		for (int j=0; j<G.SID->numSources; j++)
			sumelem1 += G.SID->sourcefracs[j]/G.SID->sumfracs()*G.SID->sourceconsts[j].constts[i];
			
		fprintf(foutput, "%i, %le, %le\n", i, sumelem1, G.SID->obsconst_mean[i]);
			
	}

	fclose(foutput);
	
		
	CSourceID SID(nparam,nconstts,nsample,nisotopes);
	SID.path = path;
	SID.getfromfile((path + constfilename).c_str());
	SID.getfromfile_isotope((path + isofilename).c_str());
	SID.elements = elements;
	SID = *G.SID;
	
	FILE *f;
	f = fopen((path+"output_MCMC.txt").c_str(), "w");
	char str[3];	
	SID.purtfact = purturbfact;
	for (int i=0; i<1; i++)
	{
		SID.iniXm = SID.sourcefracs;
		
		fprintf(f, "iteration %i\n",i);
		double sumP = SID.calcPostXYMCMC(popmcmc, interval,n_chains);
		fprintf(f, "sumP: %le\n", sumP);
		cout<<i<<endl;
	}
	fclose(f);
	
	CBTCSet MCMC_X(path + "Output_X.txt");
	CBTCSet MCMC_C(path + "Output_C.txt");
	CBTCSet MCMC_iso(path + "Output_iso_C.txt");
	vector<CBTCSet> CBTCSetConsts;
	for (int ii=0; ii<elements.size(); ii++)
	{	itoa (elements[ii], str, 10 );
		CBTCSetConsts.push_back(CBTCSet(path + "const_cont_" + string(str) + ".txt"));
	}
	

	
	vector<vector<vector<double>>> percentile_outs;
	vector<vector<double>> C_percentile_outs;
	C_percentile_outs.resize(percentiles.size());
	vector<vector<double>> C_iso_percentile_outs;
	C_iso_percentile_outs.resize(percentiles.size());

	percentile_outs.resize(elements.size()+1);
	for (int j=0; j<elements.size()+1; j++)
	{
		percentile_outs[j].resize(percentiles.size());
	}

	FILE *filepercentile;
	FILE *fileCpercentile;
	FILE *fileisopercentile;
	filepercentile = fopen((path + Percent_out_filename).c_str(),"w");
	fileCpercentile = fopen((path + C_Percent_out_filename).c_str(),"w");
	if (nisotopes>0) fileisopercentile = fopen((path + iso_Percent_out_filename).c_str(),"w");

	fprintf(filepercentile,"Total\n");
	for (int i=0; i<percentiles.size(); i++)
	{
		percentile_outs[0][i] = MCMC_X.percentile(percentiles[i],burnin);
		C_percentile_outs[i] = MCMC_C.percentile(percentiles[i],burnin);
		if (nisotopes>0) C_iso_percentile_outs[i] = MCMC_iso.percentile(percentiles[i],burnin);
		fprintf(filepercentile,"%f\n", percentiles[i]);
		CVector(percentile_outs[0][i]).writetofile(filepercentile);
		CVector(C_percentile_outs[i]).writetofile(fileCpercentile);
		if (nisotopes>0) CVector(C_iso_percentile_outs[i]).writetofile(fileisopercentile);
	}
	
	for (int j=0; j<elements.size(); j++)
	{	fprintf(filepercentile,"%i\n", elements[j]);
		for (int i=0; i<percentiles.size(); i++)
		{
			fprintf(filepercentile,"%f\n", percentiles[i]);
			percentile_outs[j+1][i] = CBTCSetConsts[j].percentile(percentiles[i],burnin);
			CVector(percentile_outs[j+1][i]).writetofile(filepercentile);
		}
	}


	fclose(filepercentile);

	return 0;


}
