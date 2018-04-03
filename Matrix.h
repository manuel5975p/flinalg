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
		/**
		Calls all the components assignment operators to with integer 0
		@brief Basic constructor for m x n matrix.
		@param m Number of rows
		@param n Number of columns
		*/
		Matrix(int _m,int _n);
		/**
		@brief Constructor for an empty matrix
		@param dummy Has absolutely no effect
		*/
		Matrix(int dummy);
		/**
		@brief Copy constructor
		*/
		Matrix(const Matrix& o);
		/**
		Assigns all components, deletes local data and assigns the other data to nullptr.
		@brief Move constructor
		*/
		Matrix(Matrix&& o);
		/**
		@brief Destructor
		*/
		~Matrix();
		/**
			Careful, this function costs quite some time for large matrices
			It is currently not more efficient than assigning to transposed()
			@see transposed()
			@brief Transposes the matrix
		*/
		void transpose();
		/**
		@brief Returns a transposed copy of this matrix
		*/
		Matrix<T> transposed();
		std::string toMatlabString() const;
		std::string toWolframString() const;
		/**
		@brief Computes the product into preallocated space
		@param space Reference to a matrix in which the result will be written
		@param other Right factor of the product
		*/
		void multInto(Matrix& space, const Matrix& other) const;
		/**
		This happens under the assumption that other is symmetric
		@brief Computes the product into preallocated space
		@param space Reference to a matrix in which the result will be written
		@param other Right factor of the product, must be symmetric
		@param secondSymmetric Dummy variable
		*/
		void multInto(Matrix& space, const Matrix& other, bool secondSymmetric) const;
		/**
		Allocates space for an ew matrix and returns the product
		@brief Computes the product
		@param other Right factor of the product
		*/
		Matrix mult(const Matrix& other) const;
		/**
		Basically this is operator-=(other)
		@brief Subtracts a matrix
		@param other The subtrahend
		*/
		void subtractInplace(const Matrix& other);
		Matrix mult(const Matrix& other, bool secondSymmetric) const;
		/**
		This happens under the assumption that *this is a block matrix of shape<br>
		I 0<br>
		0 M<br>
		with M being any kind of numbers
		@brief Computes the product into preallocated space
		@param space Reference to a matrix in which the result will be written
		@param other Right factor of the product
		@param diag_ones Count of diagonal ones in I
		*/
		void multIntoFirstdiag(Matrix& space, const Matrix& other,unsigned int diag_ones) const;
		/**
		This happens under the assumption that other is a block matrix of shape<br>
		I 0<br>
		0 M<br>
		with M being any kind of numbers
		@brief Computes the product into preallocated space
		@param space Reference to a matrix in which the result will be written
		@param other Right factor of the product
		@param diag_ones Count of diagonal ones in I
		*/
		void multIntoSeconddiag(Matrix& space, const Matrix& other,unsigned int diag_ones) const;
		/**
		Basically this is operator+=(other)
		@brief Adds a matrix
		@param other The addend
		*/
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
		/**
		@brief Swaps with another matrix
		*/
		void swap(Matrix<T>& other){
			T* temp = data;
			data = other.data;
			other.data = temp;
			std::swap(m,other.m);
			std::swap(n,other.n);
		}
		/**
		This function is mainly used in LU decomposition and solving by backsubstitution.
		@brief Adds a row multiplied with factor to another row
		@param from Source row
		@param to Destination row
		@param factor The multiplying factor
		*/
		void addRow(unsigned int from,unsigned int to,const T& factor);
		/**
		This function is mainly used in LU decomposition and solving by backsubstitution.
		@brief Adds a row multiplied with factor to another row
		@param from Source row
		@param to Destination row
		@param factor The multiplying factor
		@param startpos The first nonzero element of the source row
		*/
		void addRow(unsigned int from,unsigned int to,const T& factor, unsigned int startpos);
		/**
		@brief Permutes two rows
		*/
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
		if(pos >= m || pos < 0)
			throw std::out_of_range(std::string("Out of range ") + std::to_string(pos));
		return data + n * pos;
	}
	template<typename T>
	const T* Matrix<T>::operator[](unsigned int pos) const{
		if(pos >= m || pos < 0)
			throw std::out_of_range(std::string("Out of range ") + std::to_string(pos));
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
