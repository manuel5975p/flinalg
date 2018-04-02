#ifndef EIGEN
#define EIGEN
#include "QR.h"
#include "Matrix.h"
#include <stdexcept>
namespace core{
	/**
	Stores them into a diagonal Matrix at location D. If D points to a column vector
	of according length, the Eigenvalues will be written directly into *D.
	@brief Computes the Eigenvalues of a Matrix
	@param m The Matrix to be examined
	@param D A pointer to either a column vector or something else. In case of a column vector the eigenvalues will be written directly into *D.
	*/
	template<typename T>
	void eig(Matrix<T> m, Matrix<T>* D){
		
		if(m.m != m.n)
			throw std::invalid_argument("Matrix must be square for eig()");
		int rettype = 0;//0 for diagonal matrix, 1 for column vector
		if(D->n == 1 && D->m == m.m){
			rettype = 1;
		}
		Matrix<T> q(0),r(0);
		Matrix<T> mspace(m);
		for(int i = 0;i < 16;i++){
			QR(m,&q,&r);
			multInto(mspace,r,q);
			std::swap(mspace,m);
		}
		for(int i = 0;i < m.m;i++){
			for(int ih = i + 1;ih < m.n;ih++){
				m[i][ih] = 0;
			}
		}
		if(rettype == 1){
			for(int i = 0;i < m.m;i++)
				(*D)[i][0] = m[i][i];
		}
		else
		*D = m;
	}
	template<typename T>
	void eig(Matrix<T> m, Matrix<T>* V, Matrix<T>* D){
		if(m.m != m.n)
			throw std::invalid_argument("Matrix must be square for eig()");
	}
}
#endif