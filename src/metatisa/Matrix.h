#ifndef MAXTRIX_H
#define MAXTRIX_H

#include<iomanip>
#include"TypeDefBase.h"
#include<vector>
				
template<class Type_T>
class Matrix_T 
{
public:
	Matrix_T(unsigned rows, unsigned cols,Type_T init=Type_T(0)):
	  line(rows),colum(cols),data(rows*cols, init ){};
	  Matrix_T(){}
	inline Type_T& operator()(unsigned i,unsigned j)							//operator () 
	{
		assert(i<line&&j<colum);
		return data[i*colum+j];
	}
	inline const Type_T operator()(unsigned i,unsigned j)const					//operator() 
	{			
		assert(i<line&&j<colum);
		return data[i*colum+j];
	}
	Matrix_T<Type_T>& toBeAveraged();
	Matrix_T<Type_T>& toBeLoged();
	Matrix_T<Type_T>& toExp();
	inline int getLine()const {return line;}	
	inline int getColum()const {return colum;}
private:
	Matrix_T(std::vector<Type_T>& d,unsigned rows, unsigned cols)
		:line(rows),colum(cols),data(d){};
	std::vector<Type_T> data;
    int line, colum;
};

template<class Type_T>
Matrix_T<Type_T>& Matrix_T<Type_T>::toBeLoged()
{
	for(int col=0;col<colum;col++)
	{
		for(int lin=0;lin<line;lin++)
			if(data[ lin * colum + col ] != 0 )
				data[ lin * colum + col ] = log(data[ lin * colum + col ]);
	}
	return *this;
}

template<class Type_T>
Matrix_T<Type_T>& Matrix_T<Type_T>::toExp(){
	for(int col=0;col<colum;col++)
	{
		for(int lin=0;lin<line;lin++)
			if(data[ lin * colum + col ] != 0 )
				data[ lin * colum + col ] = exp(data[ lin * colum + col ]);
	}
	return *this;
}
template<class Type_T>
Matrix_T<Type_T>& Matrix_T<Type_T>::toBeAveraged()									
{
	for(int col=0;col<colum;++col)
	{
		Type_T sum=0;
		for(int lin=0;lin<line;lin++)
			sum+=data[ lin * colum + col ];
		for(int lin=0;lin<line;++lin)
			data[ lin * colum + col ]/=sum;
	}
//	toBeLogo();
	return *this;
}

#endif//difine MAXTRIX_H.
