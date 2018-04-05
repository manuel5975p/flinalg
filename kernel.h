#ifndef KERNEL
#define KERNEL
#include <cmath>
#include <vector>
#include <stdexcept>
#include "Matrix.h"
#include "LU.h"
#include "QR.h"
#define ZERO_THRESHOLD 0.000000001d
namespace core{
	/**
	It's basically mat.m - rank(m)
	@brief Computes the kernel dimension of mat
	@param mat The input matrix
	@return dim(ker(mat))
	*/
	template<typename T>
	unsigned int kernelDim(Matrix<T> mat){
		Matrix<T> l(0),u(0),p(0);
		LU(mat,&l,&u,&p);
		unsigned int a = 0;
		for(int i = 0;i < std::min(u.m,u.n);i++){
			if(std::abs(u[i][i]) <= ZERO_THRESHOLD){
				a++;
			}
		}
		return a;
	}
	/**
	@brief Computes a kernel basis of a matrix
	@return A kernel basis.
	@param m The input matrix
	*/
	template<typename T>
	Matrix<T> kernel(const Matrix<T>& m){
		if(m.m != m.n){
			throw std::invalid_argument("Non-square kernels not supported yet.");
		}
		Matrix<T> l(0),u(0),p(0);
		LU(m,&l,&u,&p);
		std::vector<Matrix<T>> kernelVectors;
		std::vector<unsigned int> emptyspots;
		for(int i = 0;i < std::min(u.m,u.n);i++){
			if(std::abs(u[i][i]) <= ZERO_THRESHOLD){
				emptyspots.push_back(i);
			}
		}
		for(int i = 0;i < emptyspots.size();i++){
			Matrix<T> solvec(u.m,1);
			for(int ih = 0;ih < solvec.m;ih++){
				solvec[ih][0] = (ih == emptyspots[i]);
			}
			for(int k = u.m - 1;k >= 0;k--){
				while(std::abs(u[k][k]) <= ZERO_THRESHOLD){
					k--;
					if(k < 0){goto end;}
				}
				T resid = 0;
				for(int r = k + 1;r < u.m;r++){
					resid += u[k][r] * solvec[r][0];
				}
				solvec[k][0] = -resid / u[k][k];
			}
			end:
			kernelVectors.push_back(std::move(solvec));
		}
		Matrix<T> ret(u.n,kernelVectors.size());
		for(int i = 0;i < kernelVectors.size();i++){
			for(int ih = 0;ih < kernelVectors[i].m;ih++){
				ret[ih][i] = kernelVectors[i][ih][0];
			}
		}
		return ret;
	}
}
#endif