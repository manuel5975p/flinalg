#ifndef SOLVE
#define SOLVE
#include <exception>
#include <stdexcept>
#include "Matrix.h"
namespace core{
	/**
	@brief Solves a linear equation system
	@param m Coefficient Matrix
	@param b Right side of the equation system
	*/
	Matrix<T> solve(Matrix<T> m,Matrix<T> b){
		if(m.n != m.m){
			throw std::invalid_argument("Non-square Matrix!");
		}
		if(m.m != b.m){
			throw std::invalid_argument("Dimensions not matching!");
		}
		if(b.n != 1){
			throw std::invalid_argument("Simultaneous solving not implemented yet!");
		}
		for(int i = 0;i < m.n;i++){
			ColumnPointer<T> cp = m.getColumn(i);
			if(cp[i] == 0){
				int max_i = -1;
				T max_v;
				for(int x = i + 1;x < m.m;i++){
					if(((cp[x] > 0) ? cp[x] : -cp[x]) > max_v){
						max_v = ((cp[x] > 0) ? cp[x] : -cp[x]);
						max_i = x;
					}
				}
				if(max_i == -1){
					throw std::invalid_argument("Singular Matrix");
				}
				m.permuteRows(i, max_i);
				b.permuteRows(i, max_i);
			}
			for(int ih = i + 1;ih < m.m;ih++){
				T div = -(cp[ih] / cp[i]);
				m.addRow(i,ih,div);
				b.addRow(i,ih,div);
			}
		}
		Matrix<T> sol(b.m,b.n);
		sol[b.m - 1][0] = b[b.m - 1][0] / m[b.m - 1][b.m - 1];
		for(int i = m.m - 2;i >= 0;i--){
			T scalProd = 0;
			for(int ih = i + 1;ih < m.m;ih++){
				scalProd += sol[ih][0] * m[i][ih];
			}
			sol[i][0] = (b[i][0] - scalProd) / m[i][i];
		}
		return sol;
	}
}
#endif