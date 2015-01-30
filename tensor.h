#ifndef TENSORHEADERDEF                 //only if the macro TENSORHEADERDEF is not defined execute these lines of code

#define TENSORHEADERDEF                 //define the macro

//This file contains the header and implementation for the vector3d templated class. The implementation is included in 
//this file because the class is templated. If the implementation were not included, it would also have to be included in 
//any files that depend on the header. Leaving the implementation in streamlines the compiling process.

#include <math.h>    //sqrt()
#include <iostream>  //ostream
#include "vector3d.h"

using std::ostream;
using std::endl;

//Templated class for a 2D tensor holding 9 elements
template<class T>
class tensor {
  T data[9];
  //T xx,xy,xz,yx,yy,yz,zx,zy,zz;
 public:

  //constructor
 tensor( T a, T b, T c, T d, T e, T f, T g, T h, T i) : data{a,b,c,d,e,f,g,h,i} {}
 tensor() : data{0,0,0,0,0,0,0,0,0} {}
 tensor(T i) : data{i,0,0,0,i,0,0,0,i} {} 
 tensor(vector3d<T> v1, vector3d<T> v2, vector3d<T> v3) : data{v1.X(),v1.Y(),v1.Z(),v2.X(),v2.Y(),v2.Z(),v3.X(),v3.Y(),v3.Z()} {}

  //member functions
  //operator overloads
  tensor<T> operator + (const tensor&)const;
  tensor<T> operator - (const tensor&)const;
  tensor<T> operator * (const tensor&)const;
  tensor<T> operator * (const T&)const;
  tensor<T> operator / (const T&)const;
  template <class TT>
  friend tensor<TT> operator * (const TT&, const tensor<TT>&);
  template <class TT>
  friend tensor<TT> operator / (const TT&, const tensor<TT>&);
  template <class TT>
  friend ostream & operator<< (ostream &os, const tensor<TT>&);
  //assignment of data members
  void SetXX(const T& val){data[0] = val;}
  void SetXY(const T& val){data[1] = val;}
  void SetXZ(const T& val){data[2] = val;}
  void SetYX(const T& val){data[3] = val;}
  void SetYY(const T& val){data[4] = val;}
  void SetYZ(const T& val){data[5] = val;}
  void SetZX(const T& val){data[6] = val;}
  void SetZY(const T& val){data[7] = val;}
  void SetZZ(const T& val){data[8] = val;}
  //access of data members
  T XX()const{return data[0];}
  T XY()const{return data[1];}
  T XZ()const{return data[2];}
  T YX()const{return data[3];}
  T YY()const{return data[4];}
  T YZ()const{return data[5];}
  T ZX()const{return data[6];}
  T ZY()const{return data[7];}
  T ZZ()const{return data[8];}
  //math functions
  T Trace()const{return data[0]+data[4]+data[8];}
  tensor<T> Transpose()const;
  vector3d<T> MatMult(const vector3d<T>&)const;
  tensor<T> Identity()const;
  tensor<T> Zero()const;


  //destructor
  ~tensor() {}

};

//operator overload for addition - element wise addition
template <class T>
tensor<T> tensor<T>::operator+ (const tensor &v2)const{
  tensor<T> temp = (*this);
  temp.data[0] += v2.data[0];
  temp.data[1] += v2.data[1];
  temp.data[2] += v2.data[2];

  temp.data[3] += v2.data[3];
  temp.data[4] += v2.data[4];
  temp.data[5] += v2.data[5];

  temp.data[6] += v2.data[6];
  temp.data[7] += v2.data[7];
  temp.data[8] += v2.data[8];

  return temp;
}

//operator overload for subtraction - element wise subtraction
template <class T>
tensor<T> tensor<T>::operator- (const tensor &v2)const{
  tensor<T> temp = *this;
  temp.data[0] -= v2.data[0];
  temp.data[1] -= v2.data[1];
  temp.data[2] -= v2.data[2];

  temp.data[3] -= v2.data[3];
  temp.data[4] -= v2.data[4];
  temp.data[5] -= v2.data[5];

  temp.data[6] -= v2.data[6];
  temp.data[7] -= v2.data[7];
  temp.data[8] -= v2.data[8];

  return temp;
}

//operator overload for multiplication - matrix multiplication
template <class T>
tensor<T> tensor<T>::operator* (const tensor &v2)const{
  tensor<T> temp;

  temp.data[0] = (*this).data[0] * v2.data[0] + (*this).data[1] * v2.data[3] + (*this).data[2] * v2.data[6];
  temp.data[1] = (*this).data[0] * v2.data[1] + (*this).data[1] * v2.data[4] + (*this).data[2] * v2.data[7];
  temp.data[2] = (*this).data[0] * v2.data[2] + (*this).data[1] * v2.data[5] + (*this).data[2] * v2.data[8];

  temp.data[3] = (*this).data[3] * v2.data[0] + (*this).data[4] * v2.data[3] + (*this).data[5] * v2.data[6];
  temp.data[4] = (*this).data[3] * v2.data[1] + (*this).data[4] * v2.data[4] + (*this).data[5] * v2.data[7];
  temp.data[5] = (*this).data[3] * v2.data[2] + (*this).data[4] * v2.data[5] + (*this).data[5] * v2.data[8];

  temp.data[6] = (*this).data[6] * v2.data[0] + (*this).data[7] * v2.data[3] + (*this).data[6] * v2.data[6];
  temp.data[7] = (*this).data[6] * v2.data[1] + (*this).data[7] * v2.data[4] + (*this).data[6] * v2.data[7];
  temp.data[8] = (*this).data[6] * v2.data[2] + (*this).data[7] * v2.data[5] + (*this).data[6] * v2.data[8];

  return temp;
}


//operator overload for multiplication with a scalar - element wise multiplication
template <class T>
tensor<T> tensor<T>::operator* (const T &scalar)const{
  tensor<T> temp = *this;
  temp.data[0] *= scalar;
  temp.data[1] *= scalar;
  temp.data[2] *= scalar;

  temp.data[3] *= scalar;
  temp.data[4] *= scalar;
  temp.data[5] *= scalar;

  temp.data[6] *= scalar;
  temp.data[7] *= scalar;
  temp.data[8] *= scalar;

  return temp;
}

//operator overload for multiplication with a scalar - element wise multiplication
//this function is a friend function of the class so that double * tensor behaves as tensor * double
template <class TT>
tensor<TT> operator* (const TT &scalar, const tensor<TT> &v1){
  tensor<TT> temp;
  temp.data[0] = v1.data[0] * scalar;
  temp.data[1] = v1.data[1] * scalar;
  temp.data[2] = v1.data[2] * scalar;

  temp.data[3] = v1.data[3] * scalar;
  temp.data[4] = v1.data[4] * scalar;
  temp.data[5] = v1.data[5] * scalar;

  temp.data[6] = v1.data[6] * scalar;
  temp.data[7] = v1.data[7] * scalar;
  temp.data[8] = v1.data[8] * scalar;

  return temp;
}

//operator overload for division with a scalar - element wise division
template <class T>
tensor<T> tensor<T>::operator/ (const T &scalar)const{
  tensor<T> temp = *this;
  temp.data[0] /= scalar;
  temp.data[1] /= scalar;
  temp.data[2] /= scalar;

  temp.data[3] /= scalar;
  temp.data[4] /= scalar;
  temp.data[5] /= scalar;

  temp.data[6] /= scalar;
  temp.data[7] /= scalar;
  temp.data[8] /= scalar;

  return temp;
}

//operator overload for division with a scalar - element wise division
//this function is a friend function of the class so that double / tensor works
template <class TT>
tensor<TT> operator/ (const TT &scalar, const tensor<TT> &v1){
  tensor<TT> temp;
  temp.data[0] = scalar/v1.data[0];
  temp.data[1] = scalar/v1.data[1];
  temp.data[2] = scalar/v1.data[2];

  temp.data[3] = scalar/v1.data[3];
  temp.data[4] = scalar/v1.data[4];
  temp.data[5] = scalar/v1.data[5];

  temp.data[6] = scalar/v1.data[6];
  temp.data[7] = scalar/v1.data[7];
  temp.data[8] = scalar/v1.data[8];

  return temp;
}

//operator overload for << - allows use of cout, cerr, etc.
template <class TT>
ostream & operator<< (ostream &os, const tensor<TT> &v1){
  os << v1.data[0] << ", " << v1.data[1] << ", " << v1.data[2] << endl;
  os << v1.data[3] << ", " << v1.data[4] << ", " << v1.data[5] << endl;
  os << v1.data[6] << ", " << v1.data[7] << ", " << v1.data[8] << endl;
  return os;
}

//Function to return the transpose of the given tensor
template <class T>
tensor<T> tensor<T>::Transpose()const{

  tensor<T> temp;

  temp.data[0] = (*this).data[0];
  temp.data[1] = (*this).data[3];
  temp.data[2] = (*this).data[6];

  temp.data[3] = (*this).data[1];
  temp.data[4] = (*this).data[4];
  temp.data[5] = (*this).data[7];

  temp.data[6] = (*this).data[2];
  temp.data[7] = (*this).data[5];
  temp.data[8] = (*this).data[8];

  return temp;

}

//Function to return the matrix multiplication of a tensor and vector3d
template <class T>
vector3d<T> tensor<T>::MatMult(const vector3d<T> &vec)const{

  vector3d<T> temp;

  temp.SetX( (*this).data[0] * vec.X() + (*this).data[1] * vec.Y() + (*this).data[2] * vec.Z() );
  temp.SetY( (*this).data[3] * vec.X() + (*this).data[4] * vec.Y() + (*this).data[5] * vec.Z() );
  temp.SetZ( (*this).data[6] * vec.X() + (*this).data[7] * vec.Y() + (*this).data[8] * vec.Z() );

  return temp;

}

//Function to return the identity matrix
template <class T>
tensor<T> tensor<T>::Identity()const{

  tensor<T> temp;
  T one = 1;
  T zero = 0;

  temp.data[0] = one;
  temp.data[1] = zero;
  temp.data[2] = zero;

  temp.data[3] = zero;
  temp.data[4] = one;
  temp.data[5] = zero;

  temp.data[6] = zero;
  temp.data[7] = zero;
  temp.data[8] = one;

  return temp;

}

//Function to return a tensor with all zero entries
template <class T>
tensor<T> tensor<T>::Zero()const{

  tensor<T> temp;
  return temp;

}


#endif
