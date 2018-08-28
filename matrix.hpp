#include <vector>
#include "math.h"
#include <vector>

#ifndef MATRIX
#define MATRIX

template <class T>
class matrix
{
public:
  std::vector<std::vector<T>> matData;
  int extent;

  matrix(size_t N)
  {
    extent=N;
    matData = std::vector<std::vector<T>>(std::vector<T>((size_t)N));
  }

  matrix(size_t N, T val)
  {
    extent=N;
    matData = std::vector<std::vector<T>>(std::vector<T>((size_t)N));
    for(int i=0;i<extent;i++)
 	for(int j=0;j<extent;j++)
	    i==j ? matData[i][j]=val : matData[i][j]=0;
  }

  // TODO: re-write using iterator functionality

  matrix(size_t N, std::vector<std::vector<T>> initial)
    : extent(N), matData(initial) {}
  
  std::vector<T> operator[](int i){
    return matData[i];}
  
  template <class vT>
  vT operator *(const vT vec)
  {
    vT retVec = vec;
    for(int i=0;i<extent;i++)
      {
	retVec[i]=0;
	for(int j=0;j<extent;j++)
	    retVec[i]+=matData[i][j]*vec[j];
      }
    return retVec;
  }

  matrix<T> operator *(matrix<T> mat)
  {
    matrix<T> retMat = mat;
    for(int i=0;i<extent;i++)
      {
	for(int j=0;j<extent;j++)
	  {
	    retMat[i][j]=0;
	    for(int k=0;k<extent;k++)
	      retMat[i][j]+=matData[i][k]*mat[k][j];
	  }
      }
    return retMat;
  }

  matrix<T> operator +(const matrix<T> mat)
  {
    matrix<T> retMat = mat;
    for(int i=0;i<extent;i++)
	for(int j=0;j<extent;j++)
	    retMat[i][j]= matData[i][j]+mat[i][j];
    return retMat;
  }


  matrix<T> operator -(const matrix<T> mat)
  {
    matrix<T> retMat = mat;
    for(int i=0;i<extent;i++)
	for(int j=0;j<extent;j++)
	    retMat[i][j]= matData[i][j]-mat[i][j];
    return retMat;
  }

};

#endif
