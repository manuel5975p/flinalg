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
		for(int i = 0;i < 100;i++){
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
	The eigenvalues are stored in a diagonal matrix at *D, and the eigenvectors are stored in the matrix *V.
	@brief Calculates eigenvalues and eigenvectors.
	@param m The matrix to be analyzed
	@param V Pointer to the eigenvector matrix
	@param D Pointer to the diagonal matrix containing the corresponding eigenvalues
	@throws std::invalid_argument if argument m is non square
	*/
	template<typename T>
	void eig(Matrix<T> m, Matrix<T>* V, Matrix<T>* D){
		if(m.m != m.n)
			throw std::invalid_argument("Matrix must be square for eig()");
		*D = Matrix<T>(0);
		eig(m,D);
		std::vector<Matrix<T>> eigenvectors;
		Matrix<T> v (m.m,m.n);
		int alg_multp = 1;
		for(int i = 0;i < D->m;i++){
			
			if(i < D->m - 1){
				if(std::abs((*D)[i][i] - (*D)[i + 1][i + 1]) < 0.0000001){
					continue;
				}
			}
			T ev = (*D)[i][i];
			Matrix<T> ms(m);
			for(int x = 0;x < ms.m;x++){
				ms[x][x] -= ev;
			}
			eigenvectors.push_back(std::move(kernel(ms)));
			unsigned int gi = 0;
			for(int i = 0;i < eigenvectors.size();i++){
				for(int x = 0;x < eigenvectors[i].n;x++){
					for(int ih = 0;ih < eigenvectors[i].m;ih++){
						v[ih][gi] = eigenvectors[i][ih][x];
					}
					gi++;
					/*if(gi >= v.n){
						std::cout << "Segv incoming" << std::endl;
					}*/
				}
			}
		}
		*V = std::move(v);
	}
}
#endif