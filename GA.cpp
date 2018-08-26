// GA.cpp: implementation of the CGA class.
//
//////////////////////////////////////////////////////////////////////

#include "GA.h"
#include "Individual.h"
#include "math.h"
#include <iostream>
#include "DistributionNUnif.h"
#include "Simplex.h"
#include "Matrix.h"
#include "SourceID.h"

using namespace std;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CGA::CGA()
{
	maxpop = 100;
	Ind.resize(maxpop);;
	Ind_old.resize(maxpop);
	fitdist = CDistribution(maxpop);
	N = 1;
	pcross = 1;
	cross_over_type = 1;
	MaxFitness = 0;
	SID = new CSourceID();
}

CGA::CGA(int n)
{
	maxpop = n;
	Ind.resize(maxpop);;
	Ind_old.resize(maxpop);
	fitdist = CDistribution(maxpop);
	N = 1;
	pcross = 1;
	cross_over_type = 1;
	for (int i = 0; i<15; i++)
	{		
		fixedinputs[i] = false;
	}
	MaxFitness = 0;
	//SID = new CSourceID();
}

CGA::CGA(int n, int nParam)
{
	maxpop = n;
	N = 1;
	pcross = 1;
	Ind.resize(maxpop);;
	Ind_old.resize(maxpop);
	for (int i=0; i<n; i++)
	{
		Ind[i] = CIndividual(nParam);
		Ind_old[i] = CIndividual(nParam);

	}
	fitdist = CDistribution(maxpop);
	cross_over_type = 1;
	for (int i = 0; i<15; i++)
	{		
		fixedinputs[i] = true;
	}
	MaxFitness = 0;
	SID = new CSourceID();
}

void CGA::setnparams(int n_params)
{
	Ind.resize(maxpop);;
	Ind_old.resize(maxpop);
	for (int i=0; i<maxpop; i++)
	{
		Ind[i] = CIndividual(n_params);
		Ind_old[i] = CIndividual(n_params);
	}
}


void CGA::setnumpop(int n)
{
	maxpop = n;
	CIndividual TempInd = Ind[0];
	
	int nParam = Ind[0].nParams;
	Ind.resize(maxpop);;
	Ind_old.resize(maxpop);
	for (int i=0; i<n; i++)
	{
		Ind[i] = CIndividual(nParam);
		Ind_old[i] = CIndividual(nParam);
		for (int j = 0; j<nParam; j++)
		{
			Ind[i].minrange[j] = TempInd.minrange[j];
			Ind[i].maxrange[j] = TempInd.maxrange[j];
			Ind[i].precision[j] = TempInd.precision[j];
			Ind_old[i].minrange[j] = TempInd.minrange[j];
			Ind_old[i].maxrange[j] = TempInd.maxrange[j];
			Ind_old[i].precision[j] = TempInd.precision[j];
		}

	}
	fitdist = CDistribution(maxpop);
}

CGA::CGA(const CGA &C)
{
	maxpop = C.maxpop;
	Ind.resize(maxpop);;
	Ind_old.resize(maxpop);
	for (int i=0; i<maxpop; i++)
	{
		Ind[i] = C.Ind[i];
		Ind_old[i] = C.Ind_old[i];
	}

	for (int i=0; i<15; i++)
	{
		fixedinputvale[i] = C.fixedinputvale[i];
		fixedinputs[i] = C.fixedinputs[i];
	}
	fitdist = C.fitdist;
	fitnesstype = C.fitnesstype;
	exponentcoeff = C.exponentcoeff;
	N = C.N;
	pcross = C.pcross;
	cross_over_type = C.cross_over_type;
	pmute = C.pmute;
	MaxFitness = C.MaxFitness;
	pertscale = C.pertscale;
	biasfact = C.biasfact;
	SID = C.SID;
}


CGA CGA::operator=(CGA &C)
{
	maxpop = C.maxpop;
	Ind.resize(maxpop);;
	Ind_old.resize(maxpop);
	for (int i=0; i<maxpop; i++)
	{
		Ind[i] = C.Ind[i];
		Ind_old[i] = C.Ind_old[i];
	}
	
	for (int i=0; i<15; i++)
	{
		fixedinputvale[i] = C.fixedinputvale[i];
		fixedinputs[i] = C.fixedinputs[i];
	}
	fitdist = C.fitdist;
	N = C.N;
	fitnesstype = C.fitnesstype;
	exponentcoeff = C.exponentcoeff;
	pcross = C.pcross;
	cross_over_type = C.cross_over_type;
	pmute = C.pmute;
	MaxFitness = C.MaxFitness;
	pertscale = C.pertscale;
	biasfact = C.biasfact;
	SID = C.SID;
	return *this;
	
}

CGA::~CGA()
{
	
}

void CGA::initialize()
{
	for (int i=0; i<maxpop; i++)
	{
		Ind[i].initialize();
	}


}

void CGA::Setminmax(int a, double minrange, double maxrange, int prec)
{
	for (int i=0; i<maxpop; i++)
	{	
		Ind[i].maxrange[a] = maxrange;
		Ind[i].minrange[a] = minrange;
		Ind[i].precision[a] = prec;
	}

}

void CGA::assignfitnesses()
{
	sumfitness = 0;
	int numinp = Ind[0].nParams;
	
	vector<double> inp;
	inp.resize(numinp);

	for (int k=0; k<maxpop; k++)
	
	{	int j=0;
		
		for (int i=0; i<numinp; i++)
		{
			{
				inp[i] = Ind[k].x[i];
			}
		}
		
		Ind[k].actual_fitness = 0;
		
		for (int i=0; i<numinp-2; i++)
		{	SID->sourcefracs[i] = pow(10,Ind[k].x[i]);
			SID->inc[i] = true;
		}

		SID->obsconst_std = Ind[k].x[numinp-2];
		SID->obs_iso_std = Ind[k].x[numinp-1];
		
		Ind[k].actual_fitness = SID->calcerr();
		
		
		Ind[k].fitness = pow(Ind[k].actual_fitness,N);
		//cout<<"Ind "<<k<<"  Fitness   "<<Ind[k].actual_fitness<<endl;
		sumfitness += Ind[k].fitness;
		//cout<<k<<"   "<<Ind[k].actual_fitness<<endl;
			
	}
	assignfitness_rank(N);
	
	
}


void CGA::crossover()
{
	
	for (int i=0; i<maxpop; i++)
		Ind_old[i] = Ind[i];
	int a = maxfitness();
	Ind[0] = Ind_old[a];
	Ind[1] = Ind_old[a];
	for (int i=2; i<maxpop; i+=2)
	{
		
		int j1 = fitdist.GetRand();
		int j2 = fitdist.GetRand();
		double x = GetRndUniF(0,1);
		if (x<pcross)
			if (cross_over_type == 1) 
				cross(Ind_old[j1], Ind_old[j2], Ind[i], Ind[i+1]);
			else
				cross2p(Ind_old[j1], Ind_old[j2], Ind[i], Ind[i+1]);
		else
		{	
			Ind[i] = Ind_old[j1];
			Ind[i+1] = Ind_old[j2];
		}
		
	}

}

void CGA::crossoverRC()
{
	
	for (int i=0; i<maxpop; i++)
		Ind_old[i] = Ind[i];
	int a = maxfitness();
	Ind[0] = Ind_old[a];
	Ind[1] = Ind_old[a];
	for (int i=2; i<maxpop; i+=2)
	{
		
		int j1 = fitdist.GetRand();
		int j2 = fitdist.GetRand();
		double x = GetRndUniF(0,1);
		if (x<pcross)
				cross_RC_L(Ind_old[j1], Ind_old[j2], Ind[i], Ind[i+1]);
		else
		{	
			Ind[i] = Ind_old[j1];
			Ind[i+1] = Ind_old[j2];
		}
		
	}

}

double CGA::avgfitness()
{
	double sum=0;
	for (int i=0; i<maxpop; i++)
		sum += Ind[i].fitness;
	return sum/maxpop;

}



int CGA::optimize(int nGens)
{

	FILE *FilOut;

	FilOut = fopen("GA_output.txt","w");

	initialize();
	for (int i=1; i<nGens; i++)
	{
		
		cout<<"Generation "<<i<<endl;
		assignfitnesses();
		crossover();
		mutate(pmute);
		int j = maxfitness();
		for (int k = 0; k<maxpop; k++)
			fprintf(FilOut, "%i, %lf\n", k,Ind[k].actual_fitness); 
		fprintf(FilOut, "Generation: %i,%lf, %lf, %lf", i, Ind[j].actual_fitness, Ind[j].x[0], Ind[j].x[1]);
		
		fprintf(FilOut, "\n");

	}
	assignfitnesses();
	return maxfitness();
	
	fclose(FilOut);
}


int CGA::optimize(int nGens, string RunFileName)
{

	FILE *FileOut;

	FileOut = fopen(RunFileName.c_str(),"w");
		
	double shakescaleini = shakescale;
	
	CVector X(Ind[0].nParams);
	//CSIMPLEX CSimp(Ind[0].nParams, X, 1.8);
	//*CSimp.SID = *SID;

	double ininumenhancements = numenhancements;
	numenhancements = 0;

	//initialize();
	
	CMatrix Fitness(nGens, 3);
	
	double l_MaxFitness;
	for (int i=0; i<nGens; i++)
	{
		l_MaxFitness = MaxFitness;
		assignfitnesses();
		fprintf(FileOut, "Generation: %i\n", i);
		for (int j=0; j<maxpop; j++)
		{
			fprintf(FileOut, "%i, ", j);
			
			for (int k=0; k<Ind[0].nParams; k++)
				fprintf(FileOut, "%le, ", Ind[j].x[k]);
			fprintf(FileOut, "%le, %le, %i", Ind[j].actual_fitness, Ind[j].fitness, Ind[j].rank);
			fprintf(FileOut, "\n");
		}
	
		int k = maxfitness();
		Fitness[i][0] = Ind[k].actual_fitness;

		if (i>10)
		{
			if ((Fitness[i][0]==Fitness[i-3][0]) && shakescale>pow(10.0,-Ind[0].precision[0]))
				shakescale *= shakescalered;
				
			if ((Fitness[i][0]>Fitness[i-1][0]) && (shakescale<shakescaleini))
				shakescale /= shakescalered;
				numenhancements = 0;
		}	
		if (i>50)
		{
			if ((Fitness[i][0]==Fitness[i-20][0]))
			{	numenhancements *= 1.05;
				if (numenhancements==0) numenhancements=ininumenhancements;
			}

			if ((Fitness[i][0]==Fitness[i-50][0]))
				numenhancements = ininumenhancements*10;
		}

		Fitness[i][1] = shakescale;
		Fitness[i][2] = pmute;

		if (i>20)
		{	if (shakescale == Fitness[i-20][1])
				shakescale = shakescaleini;
		}
		
		k = maxfitness();

		//enhance(num_enh, CSimp);

		k = maxfitness();
		MaxFitness = Ind[k].actual_fitness;
		Fitness[i][0] = Ind[k].actual_fitness;
		fillfitdist();
		if (RCGA == true)
			crossoverRC();
		else
			crossover();
		mutate(pmute);
		shake();
		
		cout<<i<<"    "<<Ind[k].actual_fitness<<"    "<<shakescale<<"     "<<numenhancements<<endl;
		
		if (numenhancements>0)
		{
			fprintf(FileOut, "Generation: %i (Enhanced)\n", i);
			for (int j=0; j<maxpop; j++)
			{
				fprintf(FileOut, "%i, ", j);
				
				for (int kk=0; kk<Ind[0].nParams; kk++)
					fprintf(FileOut, "%le, ", Ind[j].x[kk]);
				fprintf(FileOut, "%le, %le, %i", Ind[j].actual_fitness, Ind[j].fitness, Ind[j].rank);
				fprintf(FileOut, "\n");
			}
		}
	}
	
	assignfitnesses();
	fprintf(FileOut, "Final Enhancements\n");
	//enhance(numenhancements*100, CSimp);
	int k = maxfitness();
	cout<<Ind[k].actual_fitness<<endl;
	MaxFitness = Ind[k].actual_fitness;
	if (MaxFitness>=l_MaxFitness) pertscale/=2;
	double sumxx = 0;
	for (int kk=0; kk<Ind[0].nParams; kk++)
	{	sumxx += pow(10,Ind[k].x[kk]);
			
	}

	for (int kk=0; kk<Ind[0].nParams; kk++)
			fprintf(FileOut, "%le, ", pow(10,Ind[k].x[kk])/sumxx);
	fprintf(FileOut, "%i, %le, %le, %i, %le\n", k, Ind[k].actual_fitness, Ind[k].fitness, Ind[k].rank, pertscale);
		
	fclose(FileOut);
	return maxfitness();
	
}

void CGA::shake()
{
	for (int i=2; i<maxpop; i++)
		Ind[i].shake(shakescale);

}


void CGA::mutate(double mu)
{
	for (int i=2; i<maxpop; i++)
		Ind[i].mutate(mu);

}


int CGA::maxfitness()
{
	double max_fitness = +1E19;
	int i_max = 0;
	for (int i=0; i<maxpop; i++)
		if (max_fitness>Ind[i].actual_fitness)
		{	
			max_fitness = Ind[i].actual_fitness;
			i_max = i;
		}
	return i_max;

}

double CGA::variancefitness()
{
	double sum = 0;
	double a = avgfitness();
	for (int i=0; i<maxpop; i++)
		sum += (a - Ind[i].fitness)*(a - Ind[i].fitness);
	return sum;

}

double CGA::stdfitness()
{
	double sum = 0;
	double a = avg_inv_actual_fitness();
	for (int i=0; i<maxpop; i++)
		sum += (a - 1/Ind[i].actual_fitness)*(a - 1/Ind[i].actual_fitness);
	return sqrt(sum)/maxpop/a;

}

double CGA::avg_actual_fitness()
{
	double sum=0;
	for (int i=0; i<maxpop; i++)
		sum += Ind[i].actual_fitness;
	return sum/maxpop;

}

double CGA::avg_inv_actual_fitness()
{
	double sum=0;
	for (int i=0; i<maxpop; i++)
		sum += 1/Ind[i].actual_fitness;
	return sum/maxpop;

}



void CGA::assignrank()
{
	for (int i=0; i<maxpop; i++)
	{
		int r = 1;
		for (int j=0; j<maxpop; j++)
		{
			if (Ind[i].actual_fitness > Ind[j].actual_fitness) r++;
		}
		Ind[i].rank = r;
	}

}

void CGA::assignfitness_rank(double N)
{
	assignrank();
	double sum = 0;
	for (int i=0; i<maxpop; i++)
	{
		Ind[i].fitness = pow(1.0/static_cast<double>(Ind[i].rank),N);
		sum += Ind[i].fitness;
	}
	fitdist.s[0] = 0;
	fitdist.e[0] = Ind[0].fitness;
	for (int i=1; i<maxpop-1; i++)
	{
		fitdist.s[i] = fitdist.e[i-1];
		fitdist.e[i] = fitdist.s[i]+Ind[i].fitness;
	}
	fitdist.s[maxpop-1] = fitdist.e[maxpop-2];;
	fitdist.e[maxpop-1] = 1;

}

void CGA::enhance(CSIMPLEX &CSimp, CIndividual &InD)

{

      *CSimp.Ind = InD;
      if (numenhancements>0)
	  {
		  CSimp.nstepsimplex(numenhancements, 1.1);
		  for (int i=0; i<InD.nParams; i++)
				InD.x[i] = (CSimp.popt->vec[i]);
		  InD.actual_fitness = CSimp.fopt;
	  }
}

 

void CGA::enhance(int t_rank, CSIMPLEX &CSimp)

{

      for (int k=1; k<=t_rank; k++)
      {
            for (int j=0; j<maxpop; j++)
                  if (k == Ind[j].rank)
                  {
                        enhance(CSimp, Ind[j]);
                        j=maxpop;
                  }

      }

      assignfitness_rank(N);
     
}

void CGA::fillfitdist()
{
       double sum=0;
       for (int i=0; i<maxpop; i++)
       {
              sum+=Ind[i].fitness;
       }
 
       fitdist.s[0] = 0;
       fitdist.e[0] = Ind[0].fitness/sum;
       for (int i=1; i<maxpop-1; i++)
       {
              fitdist.e[i] = fitdist.s[i] + Ind[i].fitness/sum;
              fitdist.s[i] = fitdist.e[i];
       }
       fitdist.s[maxpop-1] = fitdist.e[maxpop-2];
       fitdist.e[maxpop-1] = 1;
 
}