#ifndef EIGEN
#define EIGEN
#include "QR.h"
#include "Matrix.h"
#include "kernel.h"
#include <vector>
#include <stdexcept>
namespace core{
	/**
	Stores them into a diagonal Matrix at location D. If D points to a column vector
	of according length, the Eigenvalues will be written directly into *D.<br>
	This function does not find complex eigenvalues
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
			//rettype = 1;
		}
		Matrix<T> q(0),r(0);
		Matrix<T> mspace(m);
		for(int i = 0;i < 40;i++){
			QR(m,&q,&r);
			r.multInto(mspace,q);
			std::swap(mspace,m);
		}
		for(int i = 0;i < m.m;i++){
			for(int ih = i + 1;ih < m.n;ih++){
				m[i][ih] = 0;
				m[ih][i] = 0;
			}
		}
		if(rettype == 1){
			for(int i = 0;i < m.m;i++)
				(*D)[i][0] = m[i][i];
		}
		else
		*D = m;
	}
	/**
	In a future far away maybe this function will work and then a documentation will follow.
	@brief Hint: Do not use this function. It's bugged and it causes segmentation faults.
	*/
	template<typename T>
	void eig(Matrix<T> m, Matrix<T>* V, Matrix<T>* D){
		std::cout << "Eigstart" << std::endl;
		if(m.m != m.n)
			throw std::invalid_argument("Matrix must be square for eig()");
		*D = Matrix<T>(0);
		eig(m,D);
		std::vector<Matrix<T>> eigenvectors;
		Matrix<T> v (m.m,m.n);
		for(int i = 0;i < D->m;i++){
			T ev = (*D)[i][i];
			Matrix<T> ms(m);
			for(int x = 0;x < ms.m;x++){
				ms[x][x] -= ev;
			}
			eigenvectors.push_back(kernel(ms));
			unsigned int gi = 0;
			for(int i = 0;i < eigenvectors.size();i++){
				for(int x = 0;x < eigenvectors[i].n;x++){
					for(int ih = 0;ih < eigenvectors[i].m;ih++){
						v[ih][gi] = eigenvectors[i][ih][x];
					}
					gi++;
				}
			}
		}
		std::cout << v.m << std::endl;
		*V = std::move(v);
		std::cout << "Eigend" << std::endl;
	}
}
#endif