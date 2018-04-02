#ifndef MATRIX 
#define MATRIX
#include <iostream>
#include <cassert>
#include <exception>
#include <string>
#include <stdexcept>
#include "ColumnPointer.h"
namespace core{
	template<typename T>
	class Matrix{
		public:
		unsigned int m, n;
		Matrix(int _m,int _n);
		Matrix(int dummy);
		Matrix(const Matrix& o);
		Matrix(Matrix&& o);
		~Matrix();
		/**
			Careful, this function costs quite some time for large matrices
			@brief Transposes the matrix
		*/
		void transpose();
		Matrix<T> transposed();
		std::string toMatlabString() const;
		std::string toWolframString() const;
		void multInto(Matrix& space, const Matrix& other) const;
		void multInto(Matrix& space, const Matrix& other, bool secondSymmetric) const;
		Matrix mult(const Matrix& other) const;
		void subtractInplace(const Matrix& other);
		Matrix mult(const Matrix& other, bool secondSymmetric) const;
		void multIntoFirstdiag(Matrix& space, const Matrix& other,unsigned int diag_ones) const;
		void multIntoSeconddiag(Matrix& space, const Matrix& other,unsigned int diag_ones) const;
		Matrix add(const Matrix& other);
		template<typename U>
		friend std::ostream& operator<<(std::ostream& out,const Matrix<U>& a);
		T* operator[](unsigned int pos);
		const T* operator[](unsigned int pos) const;
		Matrix<T> column(int c){
			Matrix<T> ret(m,1);
			for(int i = 0;i < m;i++){
				ret[i][0] = (*this)[i][c];
			}
			return ret;
		}
		Matrix<T>& operator=(Matrix<T>&& o);
		Matrix<T>& operator=(const Matrix<T>& o);
		ColumnPointer<T> getColumn(int c);
		void addRow(unsigned int from,unsigned int to);
		void swap(Matrix<T>& other){
			T* temp = data;
			data = other.data;
			other.data = temp;
			std::swap(m,other.m);
			std::swap(n,other.n);
		}
		void addRow(unsigned int from,unsigned int to,const T& factor);
		void addRow(unsigned int from,unsigned int to,const T& factor, unsigned int startpos);
		void permuteRows(unsigned int from,unsigned int to);
		private:
		T* data;
	};
	template<typename T>
	Matrix<T>::Matrix(int _m,int _n){
		m = _m;
		n = _n;
		data = new T[m*n];
		for(int i = 0;i < m*n;i++){
			data[i] = 0;
		}
	}
	template<typename T>
	Matrix<T>::Matrix(const Matrix<T>& o){
		if(o.data == nullptr){data = nullptr;m = 0;n = 0;}
		else{
			data = new T[o.n * o.m];
			for(int i = 0;i < o.m * o.n;i++){
				data[i] = o.data[i];
			}
		}
		m = o.m;
		n = o.n;
	}
	template<typename T>
	ColumnPointer<T> Matrix<T>::getColumn(int c){
		return ColumnPointer<T>(data + c, n);
	}
	template<typename T>
	Matrix<T>::Matrix(int dummy){
		data = nullptr;
		m = 0;
		n = 0;
	}
	template<typename T>
	Matrix<T>::Matrix(Matrix<T>&& o){
		data = o.data;
		o.data = nullptr;
		m = o.m;
		n = o.n;
	}
	template<typename T>
	T* Matrix<T>::operator[](unsigned int pos){
		return data + n * pos;
	}
	template<typename T>
	const T* Matrix<T>::operator[](unsigned int pos) const{
		return data + n * pos;
	}
	template<typename T>
	Matrix<T>::~Matrix(){
		if(data != nullptr)
		delete[] data;
		data = nullptr;
	}
	template<typename T>
	void Matrix<T>::transpose(){
		Matrix cp(n,m);
		for(int i = 0;i < m;i++){
			for(int ih = 0;ih < n;ih++){
				cp[i][ih] = this[0][ih][i];
			}
		}
		swap(cp);
	}
	template<typename T>
	Matrix<T> Matrix<T>::transposed(){
		Matrix cp(n,m);
		for(int i = 0;i < m;i++){
			for(int ih = 0;ih < n;ih++){
				cp[ih][i] = this[0][i][ih];
			}
		}
		return cp;
	}
	template<typename T>
	std::string Matrix<T>::toMatlabString() const{
		std::string ret;
		ret += "[";
		for(int i = 0;i < m;i++){
			for(int ih = 0;ih < n;ih++){
				ret += std::to_string((double)(*this)[i][ih]);
				if(ih < (n - 1)){
					ret += " ";
				}
			}
			ret += "; ";
		}
		ret += "]";
		return ret;
	}
	template<typename T>
	std::string Matrix<T>::toWolframString() const{
		std::string ret;
		ret += "{";
		for(int i = 0;i < m;i++){
			ret += "{";
			for(int ih = 0;ih < n;ih++){
				ret += std::to_string((double)(*this)[i][ih]);
				if(ih < (n - 1)){
					ret += ", ";
				}
			}
			ret += "}";
			if(i < (m - 1)){
				ret += ", ";
			}
		}
		ret += "}";
		return ret;
	}
	template<typename T>
	Matrix<T> Matrix<T>::mult(const Matrix<T>& other, bool secondSymmetric) const{
		assert(n == other.m);
		Matrix ret(m, other.n);
		for(int i = 0;i < ret.m;i++){
			for(int ih = 0;ih < ret.n;ih++){
				T sum = 0;
				for(int x = 0;x < n;x++){
					sum += data[i*n + x] * other.data[ih * n + x];
				}
				ret.data[i * ret.n + ih] = sum;
			}
		}
		return ret;
	}
	template<typename T>
	Matrix<T> Matrix<T>::mult(const Matrix<T>& other) const{
		assert(n == other.m);
		Matrix ret(m, other.n);
		for(int i = 0;i < ret.m;i++){
			for(int ih = 0;ih < ret.n;ih++){
				T sum = 0;
				for(int x = 0;x < n;x++){
					sum += data[i*n + x] * other.data[ih + x * other.n];
				}
				ret.data[i * ret.n + ih] = sum;
			}
		}
		return ret;
	}
	template<typename T>
	void Matrix<T>::multInto(Matrix<T>& space, const Matrix<T>& other,bool secondSymmetric) const{
		assert(n == other.m);
		assert(space.m == m && space.n == other.n);
		for(int i = 0;i < space.m;i++){
			for(int ih = 0;ih < space.n;ih++){
				T sum = 0;
				for(int x = 0;x < n;x++){
					sum += data[i * n + x] * other.data[ih * n + x];
				}
				space.data[i * space.n + ih] = sum;
			}
		}
	}
	template<typename T>
	void Matrix<T>::multInto(Matrix<T>& space, const Matrix<T>& other) const{
		assert(n == other.m);
		assert(space.m == m && space.n == other.n);
		for(int i = 0;i < space.m;i++){
			for(int ih = 0;ih < space.n;ih++){
				T sum = 0;
				for(int x = 0;x < n;x++){
					sum += data[i*n + x] * other.data[ih + x * other.n];
				}
				space.data[i * space.n + ih] = sum;
			}
		}
	}
	template<typename T>
	void Matrix<T>::multIntoSeconddiag(Matrix<T>& space, const Matrix<T>& other,unsigned int diag_ones) const{
		assert(n == other.m);
		assert(space.m == m && space.n == other.n);
		for(int i = 0;i < space.m;i++){
			for(int ih = 0;ih < diag_ones;ih++){
				space.data[i * space.n + ih] = data[i * space.n + ih];
			}
			for(int ih = diag_ones;ih < space.n;ih++){
				T sum = 0;
				for(int x = diag_ones;x < n;x++){
					sum += data[i*n + x] * other.data[ih + x * other.n];
				}
				space.data[i * space.n + ih] = sum;
			}
		}
	}
	template<typename T>
	void Matrix<T>::subtractInplace(const Matrix<T>& other){
		if(n != other.n || m != other.m)
			throw std::invalid_argument("Dimensions not matching. func: subtractInplace()");
		for(int i = 0;i < m*n;i++){
			data[i] -= other.data[i];
		}
	}
	template<typename T>
	void Matrix<T>::multIntoFirstdiag(Matrix<T>& space, const Matrix<T>& other,unsigned int diag_ones) const{
		assert(n == other.m);
		assert(space.m == m && space.n == other.n);
		for(int i = 0;i < diag_ones;i++){
			for(int ih = 0;ih < space.n;ih++){
				space.data[i * space.n + ih] = other.data[i * space.n + ih];
			}
		}
		for(int i = diag_ones;i < space.m;i++){
			for(int ih = 0;ih < space.n;ih++){
				T sum = 0;
				for(int x = 0;x < n;x++){
					sum += data[i*n + x] * other.data[ih + x * other.n];
				}
				space.data[i * space.n + ih] = sum;
			}
		}
	}
	template<typename T>
	Matrix<T> Matrix<T>::add(const Matrix<T>& other){
		if(n != other.n || m != other.m)
			throw std::invalid_argument("Dimensions not matching. func: add()");
		Matrix ret(m, other.n);
		for(int i = 0;i < ret.m;i++){
			for(int ih = 0;ih < ret.n;ih++){
				ret.data[i * ret.n + ih] = data[i * n + ih] + other.data[i * n + ih];
			}
		}
		return ret;
	}
	template<typename T>
	std::ostream& operator<<(std::ostream& out,const Matrix<T>& a){
		unsigned int space = 0;
		for(int i = 0;i < a.m;i++){
			for(int ih = 0;ih < a.n;ih++){
				unsigned int tspace = std::to_string((double)(a.data[i * a.n + ih])).size();
				space = std::max(space,tspace);
			}
		}
		for(int i = 0;i < a.m;i++){
			for(int ih = 0;ih < a.n;ih++){
				std::string frag = "";
				frag += std::to_string((double)(a.data[i * a.n + ih]));
				while(frag.size() < space){
					frag = " " + frag;
				}
				out << frag;
				if(ih < a.n - 1)
					out << " ";
			}
			out << std::endl;
		}
		return out;
	}
	template<typename T>
	void Matrix<T>::addRow(unsigned int from,unsigned int to){
		T* fr = data + from * n;
		T* t = data + to * n;
		for(int i = 0;i < n;i++){
			t[i] += fr[i];
		}
	}
	template<typename T>
	void Matrix<T>::addRow(unsigned int from, unsigned int to, const T& factor){
		T* fr = data + from * n;
		T* t = data + to * n;
		for(int i = 0;i < n;i++){
			t[i] += (fr[i] * factor);
		}
	}
	template<typename T>
	void Matrix<T>::addRow(unsigned int from, unsigned int to, const T& factor,unsigned int startpos){
		T* fr = data + from * n + startpos;
		T* t = data + to * n + startpos;
		for(int i = 0;i < n - startpos;i++){
			t[i] += (fr[i] * factor);
		}
	}
	template<typename T>
	void Matrix<T>::permuteRows(unsigned int from,unsigned int to){
		T* fr = data + from * n;
		T* t = data + to * n;
		for(int i = 0;i < n;i++){
			std::swap(fr[i],t[i]);
		}
	}
	template<typename T>
	Matrix<T>& Matrix<T>::operator=(Matrix<T>&& o){
		if(data != nullptr){
			delete[] data;
			data = nullptr;
		}
		data = o.data;
		o.data = nullptr;
		m = o.m;
		n = o.n;
	}
	template<typename T>
	Matrix<T>& Matrix<T>::operator=(const Matrix<T>& o){
		Matrix cpy(o);
		swap(cpy);
	}
	template<typename T>
	Matrix<T> identity(unsigned int n){
		Matrix<T> ret(n,n);
		for(int i = 0;i < n;i++){
			ret[i][i] = 1;
		}
		return ret;
	}
}
#endif
