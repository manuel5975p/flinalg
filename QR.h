#ifndef QRh
#define QRh
#include <exception>
#include <stdexcept>
#include <cmath>
#include "Matrix.h"
#include "ColumnPointer.h"
namespace core{
	template<typename T>
	T vectorNormSq(Matrix<T>& v){
		T sum = 0;
		if(v.n != 1){
			throw std::invalid_argument("Vectornorm called on matrix or dummy vector.");
		}
		for(int i = 0; i < v.m;i++){
			sum += v[i][0] * v[i][0];
		}
		return sum;
	}
	template<typename T>
	T vectorNorm(Matrix<T>& v){
		T sum = 0;
		if(v.n != 1){
			throw std::invalid_argument("Vectornorm called on matrix or dummy vector.");
		}
		for(int i = 0; i < v.m;i++){
			sum += v[i][0] * v[i][0];
		}
		return std::sqrt(sum);
	}
	template<typename T>
	T vectorNorm(ColumnPointer<T> cp, unsigned int length){
		T sum = 0;
		for(int i = 0; i < length;i++){
			sum += cp[i] * cp[i];
		}
		return std::sqrt(sum);
	}
	template<typename T>
	bool normCheck(Matrix<T>& c){
		for(int i = 0;i < c.m;i++){
			T sum = 0;
			for(int ih = 0;ih < c.n;ih++){
				sum += c[i][ih] * c[i][ih];
			}
			if(sum > 1.00001 || sum < 0.99999)
				return false;
		}
		for(int i = 0;i < c.n;i++){
			T sum = vectorNorm(c.getColumn(i),c.m);
			if(sum > 1.00001 || sum < 0.99999)
				return false;
		}
		return true;
	}
	template<typename T>
	T signum(T& t){
		return t >= 0 ? 1 : -1;
	}
	/**
	@brief Computes the QR decomposition of the Matrix m
	@param m The Matrix to be decomposed
	@param Q A pointer to the location where the orthogonal Matrix will be stored
	@param R A pointer to the location where the upper triangular Matrix will be stored
	*/
	template<typename T>
	void QR(Matrix<T> m, Matrix<T>* Q, Matrix<T>* R){
		Matrix<T> HHSpace(m.m,m.m);
		Matrix<T> q(m.m,m.m);
		Matrix<T> qspace(m.m,m.m);
		Matrix<T> mspace(m.m,m.n);
		for(int i = 0;i < m.m;i++){
			q[i][i] = 1;
		}
		for(int i = 0;i < std::min(m.n - 1,m.m - 1);i++){
			ColumnPointer<T> cp = m.getColumn(i);
			T norm = vectorNorm(cp + i,m.m - i);
			if(norm == 0)continue;
			Matrix<T> mirror(m.m,1);//Mirror is a column vector
			for(int x = i;x < m.m;x++){
				mirror[x][0] = cp[x];
			}
			mirror[i][0] += norm * signum(m[i][i]);
			norm = vectorNorm(mirror.getColumn(0) + i,m.m - i);
			for(int ii = 0;ii < mirror.m;ii++){
				mirror[ii][0] /= norm;
			}
			for(int ii = 0;ii < mspace.n;ii++){
				mspace[0][ii] = 0;
				for(int x = 0;x < mirror.m;x++){
					mspace[0][ii] += m[x][ii] * mirror[x][0] * 2;
				}
			}
			for(int ii = 1;ii < mspace.m;ii++){
				T mult = mirror[ii][0];
				if(mult == 0)
					for(int ih = 0;ih < mspace.n;ih++){
						mspace[ii][ih] = 0;
					}
				else
					for(int ih = 0;ih < mspace.n;ih++){
						mspace[ii][ih] = mspace[0][ih] * mult;
					}
			}
			for(int ih = 0;ih < mspace.n;ih++){
				mspace[0][ih] = mspace[0][ih] * mirror[0][0];
			}
			
			for(int ii = i;ii < q.m;ii++){
				for(int jj = i;jj < q.n;jj++){
					m[ii][jj] -= mspace[ii][jj];
				}
			}
			Matrix<T> mirspace = mirror;
			q.multInto(mirspace,mirror);
			for(int ii = 0;ii < mirspace.m;ii++){
				mirspace[ii][0] *= 2;
			}
			for(int ii = 0;ii < q.m;ii++){
				for(int jj = i;jj < q.n;jj++){//q is rectangular
					q[ii][jj] -= mirspace[ii][0] * mirror[jj][0];
				}
			}
		}
		*Q = q;
		*R = m;
	}
}
#endif