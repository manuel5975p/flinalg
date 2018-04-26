#include <random>
#include "Matrix.h"
#include "Permutation.h"
#include "ColumnPointer.h"
#include "float128.h"
#include <iostream>
#include <cmath>
#include "LU.h"
#include <iomanip>
#include "QR.h"
#include "inverse.h"
#include "util.h"
#include "kernel.h"
#include "eigen.h"
#include <chrono>
#include <exception>
#include <stdexcept>

#define AE (unsigned char)142 
#define ae (unsigned char)132 
#define OE (unsigned char)153 
#define oe (unsigned char)148 
#define UE (unsigned char)154 
#define ue (unsigned char)129 
#define ss (unsigned char)225
std::mt19937_64 gen;
std::chrono::high_resolution_clock cloc;
using qu = __float128;
using type_f = float;
type_f s = 0.4564;
std::uniform_int_distribution<int> dis(-100,100);

/**
The test namespace
*/
namespace test{
	
	/**
	@brief This function is horrible
	@return The difference of the two matrices
	@throws std::invalid_argument If the dimensions are no matching
	*/
	type_f difference(const core::Matrix<type_f>& a,const core::Matrix<type_f>& b){
		if(a.n != b.n || a.m != b.m)
			throw std::invalid_argument("Dimensions not matching, func: difference()");
		type_f sum = 0;
		for(int i = 0;i < a.m;i++){
			for(int ih = 0;ih < a.n;ih++){
				type_f ex = (a[i][ih] - b[i][ih]);
				sum += ex < 0 ? -ex : ex;
			}
		}
		return sum;
	}
	/**
	@brief Computes a squared error between to matrices
	@return The square root of the average square error
	@throws std::invalid_argument If the dimensions are no matching
	*/
	type_f avgDist(const core::Matrix<type_f>& a,const core::Matrix<type_f>& b){
		if(a.n != b.n || a.m != b.m)
			throw std::invalid_argument("Dimensions not matching, func: difference()");
		type_f sum = 0;
		for(int i = 0;i < a.m;i++){
			for(int ih = 0;ih < a.n;ih++){
				type_f ex = (a[i][ih] - b[i][ih]);
				sum += ex * ex;
			}
		}
		return (std::sqrt(sum) / a.m) / a.n;
	}
	/**
	@brief Generates a Hilbert Matrix
	*/
	template<typename T>
	core::Matrix<T> hilbertMat(unsigned int si){
		core::Matrix<type_f> a(si,si);
		for(int i = 0;i < si;i++){
			for(int ih = 0;ih < si;ih++){
				a[i][ih] = ((type_f)1) / (i + ih + 1);
			}
		}
		return a;
	}
	#define diff 0
	template<typename T>
	core::Matrix<T> randomMat(unsigned int si){
		core::Matrix<type_f> a(si,si);
		for(int i = 0;i < si;i++){
			for(int ih = 0;ih < si;ih++){
				if(ih > si - diff){
					a[i][ih] = a[i][si - diff] * 2 + a[i][si - diff - 1] * (-3);
				}
				else{
					a[i][ih] = dis(gen);
					if(a[i][ih] < 50){
						//a[i][ih] /= 1000;
					}
					else if(a[i][ih] > 50){
						//a[i][ih] *= 1000;
					}
				}
			}
		}
		return a;
	}
}
int main(int argc, char** args){
	using namespace test;
	using namespace core;
	unsigned int size = 0;
	if(argc == 2){
		size = std::stoi(args[1]);
	}
	else
	std::cin >> size;
	Matrix<type_f> a = randomMat<type_f>(size);
	Matrix<type_f> inv(0);
	std::stopwatch sw;
	inv = inverse(a);
	unsigned long long elapsed = sw.elapsed();
	std::cout << elapsed / 1000 << "mikros" << std::endl;
	Matrix<type_f> I = identity<type_f>(size);
	std::cout << "Invertierfehler: " << avgDist(I,a.mult(inv)) << std::endl;
	/*
	Matrix<type_f> l(0), u(0), p(0),q(0),r(0);
	std::stopwatch sw;
	std::cout << sw.elapsed() << std::endl;
	QR(a,&q,&r);
	auto time = sw.elapsed();
	std::cout << (double)(time / 1000) / 1000 << " ms fuer QR" << std::endl;
	sw.reset();
	LU(a,&l,&u,&p);
	time = sw.elapsed();
	std::cout << (double)(time / 1000) / 1000 << " ms fuer LU" << std::endl;
	Matrix<type_f> prod1(0),prod2(0),prod3(0);
	prod1 = q.mult(r);
	prod2 = l.mult(u);
	prod3 = p.mult(a);
	std::cout << "Multing finished" << std::endl;
	std::cout << "Fehler QR: " << avgDist(prod1, a) << std::endl;
	std::cout << "Fehler LU: " << avgDist(prod2, prod3) << std::endl;
	Matrix<type_f> d(0);Matrix<type_f> v(3,3);
	Matrix<type_f> tr = a.transposed();
	Matrix<type_f> sq = a.mult(tr);
	eig(sq,&v,&d);
	std::cout << sq.toMatlabString() << std::endl;
	std::cout << "V:\n" << v << "D:\n" << d << std::endl;*/
	/*Matrix<type_f> eig1 = v.column(0);
	std::cout << eig1 << std::endl;
	std::cout << sq.mult(eig1) << std::endl;*/
	//std::cout << a.toMatlabString() << std::endl;
	//std::cout << kernelDim(a) << std::endl;
	/*std::cout << "Eigenvectors: " << std::endl << v << std::endl;*/
}