// Matrix.cpp: implementation of the CMatrix class.
//
//////////////////////////////////////////////////////////////////////

#include "Matrix.h"
#include "math.h"
#include <iostream>


using namespace std;
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CMatrix::CMatrix(int m, int n)
{
	numrows = m;
	numcols = n;
	mat.resize(numrows);
	
	for (int i = 0;i<numrows; ++i)
	{	
		mat[i].vec.resize(numcols);
		mat[i].num = numcols;
	}
	
}

CMatrix::CMatrix()
{
	numrows = 0;
	numcols = 0;
}

CMatrix::CMatrix(int m)
{
	numrows = m;
	numcols = m;
	mat.resize(numrows);
	
	for (int i = 0;i<numrows; ++i)
	{	
		mat[i].vec.resize(numcols);
		mat[i].num = numcols;
	}
}

CMatrix::CMatrix(const CMatrix &m)
{
	numrows = m.numrows;
	numcols = m.numcols;
	mat.resize(numrows);
	
	for (int i = 0;i<numrows; ++i)
	{	
		mat[i].vec.resize(numcols);
		mat[i].num = numcols;
	}
	for (int i=0; i<numrows; ++i)  mat[i]=m.mat[i];
}

CMatrix::CMatrix(const CVector &v)
{
	numrows = v.num;
	numcols = 1;
	mat.resize(numrows);
	for (int i = 0;i<numrows; ++i)
	{	
		mat[i].vec.resize(numcols);
		mat[i].num = numcols;
	}
	
	for (int i=0; i<numrows; ++i)  mat[i].vec[0] = v.vec[i];
}


CVector& CMatrix::operator[](int i)
{
	return mat[i];
}

CMatrix::~CMatrix()
{
	mat.clear();
}

int CMatrix::getnumrows() {return numrows;};
int CMatrix::getnumcols() {return numcols;};	

CMatrix& CMatrix::operator=(const CMatrix &m)
{
	
	numcols = m.numcols;
	numrows = m.numrows;
	mat.resize(numrows);
	for (int i = 0;i<numrows; ++i)
	{	
		mat[i].vec.resize(numcols);
		mat[i].num = numcols;
	}
	
	for (int i = 0; i<numrows; ++i)
		mat[i] = m.mat[i];
	return *this;
}

CMatrix& CMatrix::operator+=(const CMatrix &m)
{
	
	for (int i=0; i<numrows; i++)
		mat[i] += m.mat[i];
	return *this;
}

CMatrix& CMatrix::operator-=(const CMatrix &m)
{
	for (int i=0; i<numrows; i++)
		mat[i] -= m.mat[i];
	return *this;
}

void CMatrix::Print(FILE *FIL)
{
	for (int i=0; i<numrows; i++)
	{	for (int j=0; j<numcols; j++)
			fprintf(FIL, "%le ", mat[i].vec[j]);
		fprintf(FIL, "\n");
	}
	fclose(FIL);

}

CMatrix operator+(const CMatrix m1, const CMatrix m2)
{
	CMatrix mt = m1;
	mt += m2;
	return mt;
}

CMatrix operator-(const CMatrix m1, const CMatrix m2)
{
	CMatrix mt = m1;
	mt -= m2;
	return mt;
}

CMatrix mult(CMatrix &m1, CMatrix &m2)
{
	int nr = m1.numrows;
	int nc = m2.numcols;
	CMatrix mt(nr,nc);
	for (int i=0; i<nr; i++)
		for (int j=0; j<nc; j++)
			for (int k=0; k<m1.numcols; k++)
				mt[i][j] += m1[i][k]*m2[k][j];
	return mt;
}

CMatrix operator*(CMatrix m1, CMatrix m2)
{
	CMatrix a= mult(m1,m2);
	return a;
}


CVector mult(CMatrix &m1, CVector &v1)
{	
	int nr = m1.numrows;
	CVector vt(nr);
	for (int i=0; i<nr; i++)
		for (int k=0; k<v1.num ; k++)
				vt[i] += m1[i][k]*v1[k];
	return vt;
}

CMatrix operator*(double d, CMatrix m1)
{
	CMatrix TrM(m1.getnumrows(), m1.getnumcols());
	for (int i=0; i<m1.getnumrows(); i++)
		for (int j=0; j<m1.getnumcols(); j++)
			TrM[i][j] = d*m1[i][j];
	return TrM;

}

CMatrix operator+(double d, CMatrix m1)
{
	CMatrix TrM(m1.getnumrows(), m1.getnumcols());
	for (int i=0; i<m1.getnumrows(); i++)
		for (int j=0; j<m1.getnumcols(); j++)
			TrM[i][j] = d+m1[i][j];
	return TrM;

}

CMatrix operator-(double d, CMatrix m1)
{
	CMatrix TrM(m1.getnumrows(), m1.getnumcols());
	for (int i=0; i<m1.getnumrows(); i++)
		for (int j=0; j<m1.getnumcols(); j++)
			TrM[i][j] = d-m1[i][j];
	return TrM;

}

CMatrix operator+(CMatrix m1, double d)
{
	CMatrix TrM(m1.getnumrows(), m1.getnumcols());
	for (int i=0; i<m1.getnumrows(); i++)
		for (int j=0; j<m1.getnumcols(); j++)
			TrM[i][j] = d+m1[i][j];
	return TrM;

}

CMatrix operator-(CMatrix m1,double d)
{
	CMatrix TrM(m1.getnumrows(), m1.getnumcols());
	for (int i=0; i<m1.getnumrows(); i++)
		for (int j=0; j<m1.getnumcols(); j++)
			TrM[i][j] = m1[i][j]-d;
	return TrM;

}

CMatrix operator/(CMatrix m1,double d)
{
	CMatrix TrM(m1.getnumrows(), m1.getnumcols());
	for (int i=0; i<m1.getnumrows(); i++)
		for (int j=0; j<m1.getnumcols(); j++)
			TrM[i][j] = m1[i][j]/d;
	return TrM;
}

CMatrix operator/(double d, CMatrix m1)
{
	CMatrix TrM(m1.getnumrows(), m1.getnumcols());
	for (int i=0; i<m1.getnumrows(); i++)
		for (int j=0; j<m1.getnumcols(); j++)
			TrM[i][j] = d/m1[i][j];
	return TrM;
}


CVector operator*(CMatrix m, CVector v)
{
	return mult(m,v);
}

//CVector gauss(CMatrix, CVector)
//{
//}

void triangulate(CMatrix &m, CVector &v)
{
	int n=m.numrows;
	for (int i=0; i<n-1; i++)
	{	double diag = m[i][i];
		for (int j=i+1; j<n; j++)
		{	double p = m[j][i]/diag;
			m[j] -= p*m[i];
			v[j] -= p*v[i];
		}
	}
}

void backsubst(CMatrix& a , CVector& b, CVector& x)
{
	int n = a.numrows;
	for (int i = n-1; i>=0; i--)
	{	double diag = a[i][i];
		x[i] = (b[i] - dotproduct(x,a[i]))/diag;
	}
}

CVector gauss0(CMatrix M, CVector V)
{	
	int n = M.numrows;
	CVector b(n);
	
	triangulate(M, V);
	backsubst(M, V, b);
	return b;
}

CVector operator/(CVector V, CMatrix M)
{
	return gauss0(M,V);
}

CMatrix Log(CMatrix &M1)
{
	CMatrix TrM(M1.getnumrows(), M1.getnumcols());
	for (int i=0; i<M1.getnumrows(); i++)
		for (int j=0; j<M1.getnumcols(); j++)
			TrM[i][j] = log(M1[i][j]);
	return TrM;
}

CMatrix Exp(CMatrix &M1)
{
	CMatrix TrM(M1.getnumrows(), M1.getnumcols());
	for (int i=0; i<M1.getnumrows(); i++)
		for (int j=0; j<M1.getnumcols(); j++)
			TrM[i][j] = exp(M1[i][j]);
	return TrM;
}

CMatrix Sqrt(CMatrix &M1)
{
	CMatrix TrM(M1.getnumrows(), M1.getnumcols());
	for (int i=0; i<M1.getnumrows(); i++)
		for (int j=0; j<M1.getnumcols(); j++)
			TrM[i][j] = sqrt(M1[i][j]);
	return TrM;
}



CMatrix Invert(CMatrix M1)
{
	CMatrix InvM(M1.getnumcols(), M1.getnumcols());
	for (int i=0; i<M1.getnumcols(); i++)
	{
		CVector V(M1.getnumcols());
		V[i] = 1;
		InvM[i] = V/M1;
	}
	return Transpose(InvM);
}

CMatrix Cholesky_factor(CMatrix &M)
{
	int i;
	int j;
	int k;
	double s;
	CMatrix b(M.getnumcols(), M.getnumcols());
	int n = M.getnumcols();
	for ( j = 0; j < n; j++ )
	{	for ( i = 0; i < n; i++ )
		{
			b[i][j] = M[i][j];
		}
	}

	for ( j = 0; j < n; j++ )
	{	for ( k = 0; k <= j-1; k++ )
		{
			for ( i = 0; i <= k-1; i++ )
			{	b[k][j] = b[k][j] - b[i][k] * b[i][j];
			}
			b[k][j] = b[k][j] / b[k][k];
		}

		s = b[j][j];
		for ( i = 0; i <= j-1; i++ )
		{
			s = s - b[i][j] * b[i][j];
		}

		b[j][j] = sqrt ( s );
	}
//
//  Since the Cholesky factor is in R8GE format, zero out the lower triangle.
//
	for ( i = 0; i < n; i++ )
	{
		for ( j = 0; j < i; j++ )
		{
			b[i][j] = 0.0;
		}
	}

	return b;
}

CMatrix CMatrix::Cholesky_factor()
{
	CMatrix M = *this;
	int i;
	int j;
	int k;
	double s;
	CMatrix b(M.getnumcols(), M.getnumcols());
	int n = M.getnumcols();
	for ( j = 0; j < n; j++ )
	{	for ( i = 0; i < n; i++ )
		{
			b[i][j] = M[i][j];
		}
	}

	for ( j = 0; j < n; j++ )
	{	for ( k = 0; k <= j-1; k++ )
		{
			for ( i = 0; i <= k-1; i++ )
			{	b[k][j] = b[k][j] - b[i][k] * b[i][j];
			}
			b[k][j] = b[k][j] / b[k][k];
		}

		s = b[j][j];
		for ( i = 0; i <= j-1; i++ )
		{
			s = s - b[i][j] * b[i][j];
		}

		b[j][j] = sqrt ( s );
	}
//
//  Since the Cholesky factor is in R8GE format, zero out the lower triangle.
//
	for ( i = 0; i < n; i++ )
	{
		for ( j = 0; j < i; j++ )
		{
			b[i][j] = 0.0;
		}
	}

	return b;
}

CMatrix LU_decomposition(CMatrix &M)
{
	double AMAX,DUM, SUM;
	int  I,IMAX,J,K;
	int d;
	int n = M.getnumrows();
	CVector VV(n);
	CMatrix A = M;
	

	for  (I=0; I<n; I++)  {
    AMAX=0.0;
    for (J=0; J<n; J++)  
		if (fabs(A[I][J]) > AMAX)  AMAX=fabs(A[I][J]);

		if (AMAX<1E-200) return A;
		VV[I] = 1/AMAX;
	}
	for (J=0; J<n;J++)  {
		for (I=0; I<J; I++)  { 
			SUM = A[I][J];
			for (K=1; K<I; K++)  
				SUM = SUM - A[I][K]*A[K][J];
			A[I][J] = SUM;
			} // i loop 
		AMAX = 0.0;

		for (I=J; I<n; I++)  {
			SUM = A[I][J];
			for  (K=0; K<J; K++)  
				SUM = SUM - A[I][K]*A[K][J];
			A[I][J] = SUM;
			DUM = VV[I]*fabs(SUM);
			if (DUM >= AMAX) {
				IMAX = I;
				AMAX = DUM;
			}
		} // i loop
   
		if (J != IMAX)  {
			for (K=0; K<n; K++)  {
				DUM = A[IMAX][K];
				A[IMAX][K] = A[J][K];
				A[J][K] = DUM;
			} // k loop 
		d = -d;
		VV[IMAX] = VV[J];
    }
   
    if (fabs(A[J][J]) < 1E-200)   A[J][J] = 1E-200;

    if (J != n)  {
      DUM = 1.0 / A[J][J];
      for (I=J+1; I<n; I++)  
        A[I][J] *= DUM;
    } 
  } // j loop 

  return A;

}

double det(CMatrix &A)
{
	CMatrix D = LU_decomposition(A);
	double prod = 1;
	for (int i=0; i<A.getnumcols(); i++)
		prod *= A[i][i];

	return prod;

}

double CMatrix::det()
{
	CMatrix A = *this;
	CMatrix D = A.LU_decomposition();
	double prod = 1;
	for (int i=0; i<A.getnumcols(); i++)
		prod *= D[i][i];

	return prod;
}

CMatrix CMatrix::LU_decomposition()
{
    CMatrix A = *this;
	double AMAX,DUM, SUM;
	int  I,IMAX,J,K;
	int d;
	int n = A.getnumrows();
	CVector VV(n);
	
	for  (I=0; I<n; I++)  {
    AMAX=0.0;
    for (J=0; J<n; J++)  
		if (fabs(A[I][J]) > AMAX)  AMAX=fabs(A[I][J]);

		if (AMAX<1E-200) return A;
		VV[I] = 1/AMAX;
	}
	for (J=0; J<n;J++)  {
		for (I=0; I<J; I++)  { 
			SUM = A[I][J];
			for (K=1; K<I; K++)  
				SUM = SUM - A[I][K]*A[K][J];
			A[I][J] = SUM;
			} // i loop 
		AMAX = 0.0;

		for (I=J; I<n; I++)  {
			SUM = A[I][J];
			for  (K=0; K<J; K++)  
				SUM = SUM - A[I][K]*A[K][J];
			A[I][J] = SUM;
			DUM = VV[I]*fabs(SUM);
			if (DUM >= AMAX) {
				IMAX = I;
				AMAX = DUM;
			}
		} // i loop
   
		if (J != IMAX)  {
			for (K=0; K<n; K++)  {
				DUM = A[IMAX][K];
				A[IMAX][K] = A[J][K];
				A[J][K] = DUM;
			} // k loop 
		//d = -d;
		VV[IMAX] = VV[J];
    }
   
    if (fabs(A[J][J]) < 1E-200)   A[J][J] = 1E-200;

    if (J != n)  {
      DUM = 1.0 / A[J][J];
      for (I=J+1; I<n; I++)  
        A[I][J] *= DUM;
    } 
  } // j loop 

  return A;

}

CVector diag(CMatrix m)
{
	CVector v(m.getnumcols());
	for (int i=0; i<m.getnumcols(); ++i)
		v[i] = m[i][i];
	return v;
}

CMatrix operator*(CVector v, CMatrix m)
{
	CMatrix a = CMatrix(v)*m;
	return a;
}


CMatrix oneoneprod(CMatrix &m1, CMatrix &m2)
{
	CMatrix TrM(m1.getnumrows(), m1.getnumcols());
	for (int i=0; i<m1.getnumrows(); i++)
		for (int j=0; j<m1.getnumcols(); j++)
			TrM[i][j] = m1[i][j]*m2[i][j];
	return TrM;
}

void CMatrix::setval(double a)
{
	for (int i=0; i<numrows ; i++)
		for (int j=0; j<numcols ; j++)
			mat[i].vec[j] = a;


}

void CMatrix::setvaldiag(double a)
{
	for (int i=0; i<getnumrows(); i++)
		mat[i].vec[i] = a;

}

void CMatrix::writetofile(FILE *f)
{
	for (int i=0; i<numrows; i++)
	{	for (int j=0; j<numcols; j++)
			fprintf(f, "%le, ", mat[i].vec[j]);
		fprintf(f, "\n");
	}
}

CMatrix Transpose(CMatrix M1)
{
	CMatrix TrM(M1.getnumcols(), M1.getnumrows());
	for (int i=0; i<M1.getnumrows(); i++)
		for (int j=0; j<M1.getnumcols(); j++)
			TrM[i][j] = M1[j][i];
	return TrM;
}

double CMatrix::sum()
{
	double sum1 = 0;
	for (int i=0; i<getnumrows(); i++)
		for (int j=0; j<getnumcols(); j++)
			sum1+=mat[i].vec[j];

	return sum1;
}

double CMatrix::product()
{
	double sum1 = 1;
	for (int i=0; i<getnumrows(); i++)
		for (int j=0; j<getnumcols(); j++)
			sum1*=mat[i].vec[j];
	return sum1;
}


