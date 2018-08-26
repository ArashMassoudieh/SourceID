
#include "Source.h"

CSource::CSource(void)
{
}

CSource::CSource(int nn)
{
	n=nn;
	n_iso = 0;
	constts.resize(n);
	constts_mean.resize(n);
	constts_var.resize(n);
	
}

CSource::CSource(int nn, int nn_iso)
{
	n=nn;
	n_iso = nn_iso;
	constts.resize(n);
	constts_mean.resize(n);
	constts_var.resize(n);
	iso.resize(n_iso);
	iso_mean.resize(n_iso);
	iso_var.resize(n_iso);

}

CSource& CSource::operator = (const CSource &C)
{
	n = C.n;
	n_iso = C.n_iso;
	constts = C.constts;
	constts_mean = C.constts_mean;
	constts_var = C.constts_var;
	iso = C.iso;
	iso_mean = C.iso_mean;
	iso_var = C.iso_var;
	return *this;
}

CSource::CSource(const CSource &C)
{
	n = C.n;
	n_iso = C.n_iso;
	constts = C.constts;
	constts_mean = C.constts_mean;
	constts_var = C.constts_var;
	iso = C.iso;
	iso_mean = C.iso_mean;
	iso_var = C.iso_var;
}

void CSource::SetnumConstts(int nn)
{
	n=nn;
	n_iso = 0;
	constts.resize(n);
	constts_mean.resize(n);
	constts_var.resize(n);
	
}

void CSource::SetnumConstts(int nn, int nn_iso)
{
	n=nn;
	n_iso = nn_iso;
	constts.resize(n);
	constts_mean.resize(n);
	constts_var.resize(n);
	iso.resize(n_iso);
	iso_mean.resize(n_iso);
	iso_var.resize(n_iso);
	
}

CSource::~CSource(void)
{
	
}
