#ifndef PERM
#define PERM
#include<vector>
#include "Matrix.h"
namespace core{
	class Permutation{
		public:
		Permutation(){}
		void permute(unsigned int a,unsigned int b){
			perms.push_back((((unsigned long long)a) << 32) | b);
		}
		template<typename T> 
		void apply(Matrix<T>& o){
			for(int i = 0;i < perms.size();i++){
				o.permuteRows((unsigned int)(perms[i] >> 32),(unsigned int)perms[i]);
			}
		}
		template<typename T>
		Matrix<T> toMatrix(unsigned int n){
			Matrix<T> ret = identity<T>(n);
			apply(ret);
			return ret;
		}
		private:
		std::vector<unsigned long long> perms;
	};
}
#endif