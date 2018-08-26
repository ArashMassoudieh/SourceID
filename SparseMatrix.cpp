// SpareMatrix.cpp: implementation of the CSpareMatrix class.
//
//////////////////////////////////////////////////////////////////////

#include "Matrix.h"
#include "SparseMatrix.h"
#include "Vector.h"
#include "math.h"


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

int max(int i, int j)
{
	if (i>=j)
		return i;
	else
		return j;
}


CSpareMatrix::CSpareMatrix(int n_n, int n_nnz)
{
	nnz = n_nnz;
	n = n_n;
	value = new double[nnz];
	colsub = new int[n];
	rowind = new int[nnz];
	for (int i=0; i<nnz; i++)
	{
		rowind[i] = n+1;
		value[i] = 0;
	}
	for (int i=0; i<n; i++)
	{
		colsub[i] = 0;
	}
}

CSpareMatrix& CSpareMatrix::operator=(const CSpareMatrix& SM)
{
	nnz = SM.nnz;
	n = SM.n;
	value = new double[nnz];
	colsub = new int[n];
	rowind = new int[nnz];
	for (int i=0; i<nnz; i++)
	{
		rowind[i] = SM.rowind[i];
		value[i] = SM.value[i];
	}
	for (int i=0; i<n; i++)
	{
		colsub[i] = SM.colsub[i];
	}
	return *this;
}

CSpareMatrix::CSpareMatrix(CSpareMatrix& SM)
{
	nnz = SM.nnz;
	n = SM.n;
	value = new double[nnz];
	colsub = new int[n];
	rowind = new int[nnz];
	for (int i=0; i<nnz; i++)
	{
		rowind[i] = SM.rowind[i];
		value[i] = SM.value[i];
	}
	for (int i=0; i<n; i++)
	{
		colsub[i] = SM.colsub[i];
	}
}


CSpareMatrix::CSpareMatrix(int n_n)
{
	n = n_n;
	colsub = new int[n];
	for (int i=0; i<n; i++)
	{
		colsub[i] = 0;
	}
}

CSpareMatrix::CSpareMatrix()
{

}

CSpareMatrix::~CSpareMatrix()
{
	delete[] value;
	delete[] rowind;
	delete[] colsub;
}

double CSpareMatrix::getval(int ii,int ji)
{
	
	if (ii<0 || ii>n-1 || ji<0 || ji>n-1)
		return 0;
	else
	{
		int col2s;
		if (ji<n-1)
			col2s = colsub[ji+1];
		else
			col2s = nnz;

		for (int i=colsub[ji]; i<col2s; i++)
		{	
			if (ii==rowind[i])
				return value[i];
		}
		return 0;
	}
}

void CSpareMatrix::setval(int ii,int ji, double v)
{
	bool doesexist = false;
	int col2s;
	if (ji<n-1)
		col2s = colsub[ji+1];
	else
		col2s = nnz;
	for (int i=colsub[ji] ; i<col2s; i++)
		if (rowind[i] == ii)
		{	
			doesexist = true;
			value[i] = v;
		}
	if (doesexist == false)
	{
		if (ji<n-1)
			col2s = colsub[ji+1];
		else
			col2s = nnz;
	
		int shiftpoint=col2s;
		
		for (int j=colsub[ji]; j<col2s; j++)
			if (ii<rowind[j])
			{	
				shiftpoint = j;
				if (rowind[j] == n+1)
					j=col2s;
			}
		for (int j=nnz-2; j>shiftpoint-1; --j)
		{
			rowind[j+1] = rowind[j];
			value[j+1] = value[j];
		}
		rowind[shiftpoint] = ii;
		value[shiftpoint] = v;
		for (int i=ji+1; i<n; i++)
			colsub[i]++;

	}

}


CSpareMatrix mult(CSpareMatrix A, CSpareMatrix B)
{
	int nnnz = 0;
	int col2s;
	for (int i=0; i<A.n; i++)
		for (int j=0; j<A.n; j++)
		{	double sum = 0;
			if (j<B.n-1)
				col2s = B.colsub[j+1];
			else
				col2s = B.nnz;
			for (int is = B.colsub[j]; is<col2s ; is++)
				sum += B.value[is]*A.getval(i,A.rowind[is]);
			if (sum!=0)
				nnnz++;
		}
	
	CSpareMatrix V(A.n, nnnz);

	for (int i=0; i<A.n; i++)
		for (int j=0; j<A.n; j++)
		{	double sum = 0;
			if (j<B.n-1)
				col2s = B.colsub[j+1];
			else
				col2s = B.nnz;
			for (int is = B.colsub[j]; is<col2s ; is++)
				sum += B.value[is]*A.getval(i,A.rowind[is]);
			if (sum!=0)
				V.setval(i,j,sum);
		}
	return V;

}

int CSpareMatrix::colass(int i)
{
	for (int j=0; j<nnz; j++)
	{	
		if (colsub[j]>i)
			return j-1;
	}	
	return -1;
}


CVector mult(CSpareMatrix A, CVector B)
{
	CVector V(A.n);
	int col2s;
	for (int j=0; j<A.n; j++)
	{
		if (j<A.n-1)
			col2s = A.colsub[j+1];
		else
			col2s = A.nnz;
		for (int is = A.colsub[j]; is<col2s ; is++)	
			if (A.rowind[is]<A.n)
				V[j] += A.value[is]*B[A.rowind[is]];
	}

	return V;

}

CVector diag(CSpareMatrix A)
{
	CVector V(A.n);
	for (int i=0; i<A.n; i++)
		V[i] = A.getval(i,i);
	return V;

}

CSpareMatrix operator+(CSpareMatrix A1, CSpareMatrix A2)
{
	int nnnz = 0;
	for (int i=0; i<A1.n; i++)
		for (int j=0; j<A1.n; j++)
			if (A1.getval(i,j) + A2.getval(i,j) != 0)
				nnnz++;

	CSpareMatrix A(A1.n, nnnz);
	
	for (int i=0; i<A1.n; i++)
		for (int j=0; j<A1.n; j++)
			if (A1.getval(i,j) + A2.getval(i,j) != 0)
				A.setval(i,j,A1.getval(i,j)+A2.getval(i,j));

	return A;
}



CSpareMatrix operator*(CSpareMatrix A1, CSpareMatrix A2)
{
	return mult(A1, A2);

}


CVector operator*(CSpareMatrix A1, CVector V)
{
	return mult(A1, V);

}

CVector CSpareMatrix::getcol(int j)
{
	int col2s;
	CVector V(n);
	if (j<n-1)
		col2s = colsub[j+1];
	else
		col2s = nnz;
	for (int is = colsub[j]; is<col2s ; is++)
		V[rowind[is]] = value[is];

	return V;
	

}


CVector CSpareMatrix::getrow(int i)
{
	CVector V(n);
	for (int j = 0; j<n ; j++)
		V[j] = getval(i,j);

	return V;
}

CVector RichardsonIteration(CSpareMatrix A, CVector B, double tol)
{
	
	CVector r(A.n);
	CVector x(A.n);
	CVector z(A.n);

	double normb = norm(B);
	if (normb == 0) normb = 1;
	r = B-A*x;
	int i = 0;
	while ((norm(r)/normb>tol) && (i<10000))
	{
		i++;	
		z = r/diag(A);
		x+=z;
		r=B-A*x;
	}
	if (i>=10000)
	{	
		for (int kk=0; kk<B.num; kk++)
			B[kk] = -999;
		return B;
	}

	return x;
}

CSpareMatrix CSpareMatrix::Transpose()
{
	CSpareMatrix B(n, nnz);
	int col2s;
	for (int j=0; j<n; j++)
	{
		if (j<n-1)
			col2s = colsub[j+1];
		else
			col2s = nnz;
		for (int is = colsub[j]; is<col2s ; is++)	
			if (rowind[is]<n)
				B.setval(j, rowind[is],value[is]);
	}
	
	return B;
}

CVector operator/(CVector V, CSpareMatrix A)
{
	CSpareMatrix A1 = A.Transpose();
	return RichardsonIteration(A1, V, 1e-10);
}

CVector BiCG(CSpareMatrix A, CVector B, double tol)
{
	double resid;
	double ro_1, ro_2, alpha, beta;
	CVector z(A.n), ztilde(A.n), p(A.n), ptilde(A.n), q(A.n), qtilde(A.n);
	CVector x(A.n);

	double normb = norm(B);
	CVector r = B-A*x;
	CVector rtilde = r;

	if (normb == 0) 
		normb = 1;

	if ((resid=norm(r)/normb) <= tol)
		return x;
	
	int i=0;
	while ((resid>tol) && (i<2*A.n)) 
	{
		z = r/diag(A);
		ztilde = r/diag(A);
		ro_1 = dotproduct(z, rtilde);
		
		if (i==0)
		{
			p = z;
			ptilde = ztilde;
		}
		else
		{
			beta = ro_1/ro_2;
			p = z+beta*p;
			ptilde = ztilde+beta*ptilde;
		}

		q = A*p;
		qtilde = A.Transpose()*ptilde;
		alpha = ro_1/dotproduct(ptilde, q);
		x += alpha*p;
		r += -alpha*q;
		rtilde += -alpha*qtilde;
		ro_2 = ro_1;
		i++;
		resid = norm(r)/normb;
	}

	return x;
}
	
void CSpareMatrix::addtoval(int i,int j, double val)	
{
	if (val != 0.0)
		setval(i,j,getval(i,j)+val);
}

double CSpareMatrix::sumrows(int i)
{
	double sum = 0;
	for (int j = 0; j<n; j++)
		sum += getval(i,j);
	return sum;
	
}
	

double CSpareMatrix::sumcols(int i)
{
	double sum = 0;
	for (int j = 0; j<n; j++)
		sum += getval(j,i);
	return sum;
}

CMatrix inverse(CSpareMatrix A)
{
	CMatrix M(A.n, A.n);
	for (int i=0; i<A.n; i++)
	{	
		CVector V(A.n);
		V[i] = 1;
		M[i] = RichardsonIteration(A, V, 1e-12);
	}
	return M;
}

