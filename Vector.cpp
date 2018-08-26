// Vector.cpp: implementation of the CVector class.
//
//////////////////////////////////////////////////////////////////////

#include "Vector.h"
#include "Math.h"
#include "Matrix.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CVector::CVector()
{
	num=0;
}

CVector::CVector(int n)
{
	num = n;
	vec.resize(num);
}

CVector::CVector(const vector<double> a, int n)
{
	num = n;
	vec = a;
}

CVector::CVector(const double x, int n)
{
	num = n;
	vec.resize(num);
	for (int i=0; i<num; i++) vec[i] = x;
}

CVector::CVector(const CVector &v)
{
	num = v.num;
	vec = v.vec;
}

CVector::CVector(const vector<double> &v)
{
	num = v.size();
	vec = v;
}

CVector::~CVector()
{
	vec.clear();
}

double& CVector::operator[](int i)
{
	double *p = 0;
	if ((i<num) & (i>-1))
		return vec[i];
	else
		return *p;
}

int CVector::range(int i)
{
	return i;
}

CVector& CVector::operator=(const CVector &v)
{
	num = v.num;
	vec = v.vec;
	return *this;
}

CVector& CVector::operator=(const vector<double> &v)
{
	num = v.size();
	vec = v;
	return *this;
}

CVector& CVector::operator=(const double &v)
{
	for (int i=0; i<num; ++i)
		vec[i] = v;
	return *this;
}

CVector& CVector::operator+() 
{return *this;}

void CVector::swap(int i, int j)
{	double tmp = vec[range(i)];
	vec[i] = vec[range(j)];
	vec[j] = tmp;
	
}

int CVector::getsize() {return num;}

CVector& CVector::operator*=(double x)
{
	for (int i=0; i<num; ++i)
		vec[i] *= x;
	return *this;

}

CVector& CVector::operator/=(double x)
{
	for (int i=0; i<num; ++i)
		vec[i] /= x;
	return *this;

}

CVector& CVector::operator+=(const CVector &v)
{
	for (int i=0; i<num; ++i)
		vec[i] += v.vec[i];
	return *this;
}

CVector& CVector::operator-=(const CVector &v)
{
	for (int i=0; i<num; ++i)
		vec[i] -= v.vec[i];
	return *this;
}

CVector operator+(CVector v1, CVector v2) 
{
	return v1 += v2;
}

CVector operator-(CVector v1, CVector v2) 
{
	return v1 -= v2;
}

double dotproduct(CVector v1, CVector v2) 
{
	double d;
	if (v1.num = v2.num) 
	{
	d = 0;
	for (int i=0; i<v1.num; ++i)
		d += v1.vec[i]*v2.vec[i];
	return d;
	}
}

CVector& CVector::operator*=(const CVector& v)
{
	for (int i=0; i<num; ++i)
		vec[i] *= v.vec[i];
	return *this;
}
	

CVector operator*(CVector v1, CVector v2) 
{
	return v1 *= v2;
}

double norm(CVector v)
{
	double sum = 0;
	for (int i=0; i<v.num; i++)
		sum += pow(v.vec[i],2);
	return sqrt(sum);
}

CVector operator*(double a, CVector v)
{
	return v*=a;

}

CVector operator/(CVector v, double a)
{
	return v*=(1/a);
}

CVector operator+(CVector v, double a)
{
	CVector v1(v.num);
	for (int i=0; i<v.num; i++)
		v1[i] = a + v[i];
	return v1;
}

CVector operator+(double a, CVector v)
{
	CVector v1(v.num);
	for (int i=0; i<v.num; i++)
		v1[i] = a + v[i];
	return v1;
}

CVector operator-(double a, CVector &v)
{
	CVector v1(v.num);
	for (int i=0; i<v.num; i++)
		v1[i] = a - v[i];
	return v1;

}


CVector operator/(CVector v1, CVector v2)
{
	CVector x(v1.getsize());
	for (int i = 0; i<v1.getsize(); ++i)
		x[i] = v1[i]/v2[i];
	return x;
}

CVector operator/(double a, CVector v2)
{
	CVector x(v2.getsize());
	for (int i = 0; i<v2.getsize(); ++i)
		x[i] = a/v2[i];
	return x;

}

bool CVector::operator==(double v)
{
	bool r=true;
	for (int i=0; i<num; ++i)
		if (vec[i] != v)
			r=false;
	return r;

}

double CVector::max()
{
	double a = -1E14;
	for (int i=0;i<num; i++)
	{
		if (vec[i]>a)
			a = vec[i];
	}
	return a;

}

double CVector::min()
{
	double a = 1E14;
	for (int i=0;i<num; i++)
	{
		if (vec[i]<a)
			a = vec[i];
	}
	return a;

}


double CVector::abs_max()
{
	double a = -1E14;
	for (int i=0;i<num; i++)
	{
		if (fabs(vec[i])>a)
			a = fabs(vec[i]);
	}
	return a;
}


double CVector::norm2()
{
	double a = 0;
	for (int i=0;i<num; i++)
	{
		a+=vec[i]*vec[i];
	}
	return pow(a,0.5);

}

double CVector::sum()
{
		double a = 0;
	for (int i=0;i<num; i++)
	{
		a+=vec[i];
	}
	return a;
}

CMatrix CVector::T()
{
	CMatrix K(1,num);
	for (int i=0; i<num; i++)
		K[0][i] = vec[i];
	return K;
}

CVector CVector::Log()
{
	CVector x(getsize());
	for (int i = 0; i<getsize(); i++)
		x[i] = log(vec[i]);
	return x;
}

CVector Log(CVector &V)
{
	return V.Log();

}

CVector CVector::Exp()
{
	CVector x(getsize());
	for (int i = 0; i<getsize(); i++)
		x[i] = exp(vec[i]);
	return x;
}

CVector Exp(CVector &V)
{
	return V.Exp();

}

void CVector::writetofile(FILE *f)
{
	for (int i=0; i<num; i++)
		fprintf(f, "%le, ", vec[i]);
	fprintf(f, "\n");
}

CMatrix CVector::diagmat()
{
	CMatrix A(num,num);
	for (int i=0; i<num; i++)
		A[i][i] = vec[i];

	return A;

}

CVector zeros(int i)
{
	CVector V(i);
	return V;
}