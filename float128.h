#ifndef LFLOAT
#define LFLOAT
#include <quadmath.h>
#include <string>
namespace std{
	/**
		@return The square root of arg
	*/
	inline __float128 sqrt(const __float128& arg){
		return sqrtq(arg);
	}
	/**
		@return The absolute value of arg
	*/
	inline __float128 abs(const __float128& arg){
		if(arg >= 0)
		return arg;
		return -arg;
	}
}
#endif