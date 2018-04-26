#ifndef SVD_H
#define SVD_H
#include "Matrix.h"
#include "eigen.h"
template<typename T>
void svd(Matrix<T>& m,Matrix<T>* U,Matrix<T>* S,Matrix<T>* V){
	Matrix<T> mt = m.transposed();
	Matrix<T> pre = m.mult(mt);
	Matrix<T> post = mt.mult(m);
	Matrix<T> u(0);
	Matrix<T> v(0);
}
#endif