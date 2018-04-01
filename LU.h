#ifndef LUh
#define LUh
#include "Matrix.h"
#include "Permutation.h"
#include "ColumnPointer.h"
/**
This namespace constains the implementations of QR, LU and singular value decomposition algorithms. Further it contains function templates for solving linear equation systems numerically
and calculating the eigensystem of matrices.
@brief The core namespace contains the core functions and classes such as the Matrix base class template and functions for basic decompositions.
*/
namespace core{
	std::ostream& operator<<(std::ostream& out,__float128 f){
		if(std::abs((double)(f)) < 0.00000001)return out << "0";
			std::cout << (double)(f);
		return out;
	}
	template<typename T>
	void LU(Matrix<T> m, Matrix<T>* L, Matrix<T>* U, Matrix<T>* P){
		int crank = 0;
		Matrix<T> l(m.m,m.m);
		Permutation p;
		for(int i = 0;i < m.n;i++){
			ColumnPointer<T> cp = m.getColumn(i);
			T a = cp[crank];
			if(true){		
				int max_i = -1;
				T max_v = 0;
				__builtin_prefetch(&(m[0][0]) + m.n * (crank + 1),1,0);
				for(int x = crank + 1;x < m.m;x++){
					if(((cp[x] > 0) ? cp[x] : -cp[x]) > max_v){
						max_v = ((cp[x] > 0) ? cp[x] : -cp[x]);
						max_i = x;//TODO: Eliminate copying into max_v
					}
				}
				if(max_v == 0){
					crank++;
					if(crank == m.m - 1)break;
					continue;
				}
				p.permute(crank, max_i);
				l.permuteRows(crank, max_i);
				m.permuteRows(crank, max_i);
			}
			for(int ih = crank + 1;ih < m.m;ih++){
				T div = -(cp[ih] / cp[crank]);
				m.addRow(crank,ih,div,i);
				l[ih][i] = -div;
			}
			crank++;
			if(crank == m.m - 1)break;
		}
		for(int i = 0; i < l.n;i++){
			l[i][i] = 1;
		}
		if(P != nullptr){
			*P = p.toMatrix<T>(m.m);
		}
		*L = std::move(l);
		*U = std::move(m);
	}
}
#endif