#ifndef INVERSE_H
#define INVERSE_H
#include "Matrix.h"
namespace core{
	/**
	@brief Computes the inverse
	@param m The Matrix
	*/
	template<typename T>
	Matrix<T> inverse(Matrix<T> m){
		if(m.m != m.n)
			throw std::invalid_argument("Matrix must be square for inverse()");
		Matrix<T> I = identity<T>(m.m);
		for(int i = 0;i < m.n - 1;i++){
			ColumnPointer<T> cp = m.getColumn(i);
			int max_i = -1;
			T max_v = 0;
			for(int x = i;x < m.m;x++){
				if(((cp[x] > 0) ? cp[x] : -cp[x]) > max_v){
					max_v = ((cp[x] > 0) ? cp[x] : -cp[x]);
					max_i = x;//TODO: Eliminate copying into max_v
				}
			}
			if(max_v == 0){
				throw std::invalid_argument("Matrix must not be singular for inverse()");
			}
			I.permuteRows(i, max_i);
			m.permuteRows(i, max_i);
			for(int ih = i + 1;ih < m.m;ih++){
				T div = -(cp[ih] / cp[i]);
				m.addRow(i,ih,div,i);
				I.addRow(i,ih,div,i);
			}
		}
		//return I;
		for(int i = m.n - 1;i >= 0;i--){
			ColumnPointer<T> cp = m.getColumn(i);
			int max_i = -1;
			T max_v = 0;
			for(int x = i;x >= 0;x--){
				if(((cp[x] > 0) ? cp[x] : -cp[x]) > max_v){
					max_v = ((cp[x] > 0) ? cp[x] : -cp[x]);
					max_i = x;//TODO: Eliminate copying into max_v
				}
			}
			I.permuteRows(i, max_i);
			m.permuteRows(i, max_i);
			for(int ih = i - 1;ih >= 0;ih--){
				T div = -(cp[ih] / cp[i]);
				m.addRow(i,ih,div,0,i + 1);
				I.addRow(i,ih,div,0,i + 1);
			}
		}
		for(int i = 0;i < I.n;i++){
			T mult = 1/m[i][i];
			for(int x = 0;x < I.m;x++)
			I[i][x] *= mult;
		}
		return I;
	}
}
#endif