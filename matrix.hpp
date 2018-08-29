#include <vector>
#include "math.h"
#include <vector>
#include <iterator>

#ifndef MATRIX
#define MATRIX


/// Class representing a matrix supporting simple multiplication and
/// multiplication with vector types
template <class T>
class matrix
{
public:
  std::vector<std::vector<T>> matData;///< The storage for the NxN matrix
  int extent; ///< Size of each dimension of the square matrix

  /// initializes an empty matrix of size NxN
  /// \param N size of the matrix
  matrix(size_t N)
  {
    extent=N;
    matData = std::vector<std::vector<T>>(N,std::vector<T>(N));
  }

  /// initializes a diagonal matrix of size NxN with value T on the diagonal
  /// \param N extent of matrix
  /// \param val value on the diagonal
  matrix(size_t N, T val)
  {
    extent=N;
    matData = std::vector<std::vector<T>>(N,std::vector<T>(N));
    for(int i=0;i<extent;i++)
      for(int j=0;j<extent;j++)
	i==j ? matData[i][j]=val : matData[i][j]=0;
  }

  /// initalizes a diagonal matrix of size NxN with a full initial matrix
  /// \param N extent of matrix
  /// \param initial data to be copied to the internal matrix data
  matrix(size_t N, std::vector<std::vector<T>> initial)
    : extent(N), matData(initial) {}

  /// returns the iterator begin of the interior data type
  auto begin(){
    return matData.begin();}

  /// returns the iterator end of the interior data type
  auto end(){
    return matData.end();}
  
  /// extracts the vector at a particular row
  /// \param i row to extract
  std::vector<T> operator[](int i){
    return matData[i];}

  /// multiply a matrix by a vector, return same vector type
  /// \param vec vector to multiply
  /// \return resulting product vector
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

  /// multiply two matrices together
  /// \param mat matrix to multiply (right-left as expected)
  /// \return the product matrix
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

  /// add two matrices together
  /// \param mat summand
  /// \return sum
  matrix<T> operator +(const matrix<T> mat)
  {
    matrix<T> retMat = mat;
    for(int i=0;i<extent;i++)
	for(int j=0;j<extent;j++)
	    retMat[i][j]= matData[i][j]+mat[i][j];
    return retMat;
  }

  /// subtract two matrices
  /// \param mat differand
  /// \return difference
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
