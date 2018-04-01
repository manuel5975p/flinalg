#ifndef EIGEN
#define EIGEN
#include "QR.h"
#include <stdexcept>
namespace core{
	template<typename T>
	void eig(Matrix<T> m, Matrix<T>* D){
		if(m.m != m.n)
			throw std::invalid_argument("Matrix must be square for eig()");
		
	}
	template<typename T>
	void eig(Matrix<T> m, Matrix<T>* V, Matrix<T>* D){
		if(m.m != m.n)
			throw std::invalid_argument("Matrix must be square for eig()");
	}
}
#endif