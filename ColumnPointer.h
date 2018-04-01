#ifndef CP
#define CP
namespace core{	
	template<typename T>
	class ColumnPointer{
		private:
		T* pointer;
		int rowlength;
		public:
		ColumnPointer(T* point,int rowlen);
		ColumnPointer<T> operator+(int plus);
		ColumnPointer<T> operator++();
		const T& operator[](int pos)const;
		T& operator[](int pos);
		T operator*() const;
		ColumnPointer<T>& operator+=(int plus);
	};
	
	template<typename T>
	ColumnPointer<T> ColumnPointer<T>::operator+(int plus){
		return ColumnPointer<T>(pointer + (plus * rowlength),rowlength); 
	}
	template<typename T>
	ColumnPointer<T>& ColumnPointer<T>::operator+=(int plus){
		pointer += (plus * rowlength,rowlength);
		return *this; 
	}
	template<typename T>
	ColumnPointer<T>::ColumnPointer(T* point,int rowlen){
		pointer = point;
		rowlength = rowlen;
	}
	template<typename T>
	ColumnPointer<T> ColumnPointer<T>::operator++(){
		ColumnPointer<T> ret = ColumnPointer<T>(pointer,rowlength);
		pointer += rowlength;
		return ret;
	}
	template<typename T>
	T ColumnPointer<T>::operator*() const{
		return *pointer;
	}
	template<typename T>
	const T& ColumnPointer<T>::operator[](int pos) const{
		return *(pointer + pos * rowlength);
	}
	template<typename T>
	T& ColumnPointer<T>::operator[](int pos) {
		return *(pointer + pos * rowlength);
	}
}
#endif