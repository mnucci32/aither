#ifndef VECTOR3DHEADERDEF                 //only if the macro VECTOR3DHEADERDEF is not defined execute these lines of code

#define VECTOR3DHEADERDEF                 //define the macro

//This file contains the header and implementation for the vector3d templated class. The implementation is included in 
//this file because the class is templated. If the implementation were not included, it would also have to be included in 
//any files that depend on the header. Leaving the implementation in streamlines the compiling process.

#include <math.h>    //sqrt()
#include <iostream>  //ostream

using std::ostream;

//Templated class for a vector holding 3 entries
template<class T>
class vector3d {
  T x,y,z;
 public:

  //constructor
  vector3d( T a, T b, T c) : x(a), y(b), z(c) {}
  vector3d() : x(0), y(0), z(0) {}

  //member functions
  //operator overloads
  vector3d<T> operator + (const vector3d&)const;
  vector3d<T> operator - (const vector3d&)const;
  vector3d<T> operator * (const T&)const;
  vector3d<T> operator / (const T&)const;
  template <class TT>
  friend vector3d<TT> operator * (const TT&, const vector3d<TT>&);
  template <class TT>
  friend vector3d<TT> operator / (const TT&, const vector3d<TT>&);
  template <class TT>
  friend ostream & operator<< (ostream &os, const vector3d<TT>&);
  //assignment of data members
  void SetX(const T& val){x = val;}
  void SetY(const T& val){y = val;}
  void SetZ(const T& val){z = val;}
  //access of data members
  T X()const{return x;}
  T Y()const{return y;}
  T Z()const{return z;}
  //math functions
  T DotProd( const vector3d&)const;
  vector3d<T> CrossProd(const vector3d&)const;
  inline T Mag()const;
  inline T MagSq()const;
  T SumElem()const;
  T Distance( const vector3d&)const;

  //destructor
  ~vector3d() {}

};

//operator overload for addition - element wise addition
template <class T>
vector3d<T> vector3d<T>::operator+ (const vector3d &v2)const{
  vector3d<T> temp = *this;
  temp.x += v2.x;
  temp.y += v2.y;
  temp.z += v2.z;
  return temp;
}

//operator overload for subtraction - element wise subtraction
template <class T>
vector3d<T> vector3d<T>::operator- (const vector3d &v2)const{
  vector3d<T> temp = *this;
  temp.x -= v2.x;
  temp.y -= v2.y;
  temp.z -= v2.z;
  return temp;
}

//operator overload for multiplication with a scalar - element wise multiplication
template <class T>
vector3d<T> vector3d<T>::operator* (const T &scalar)const{
  vector3d<T> temp = *this;
  temp.x *= scalar;
  temp.y *= scalar;
  temp.z *= scalar;
  return temp;
}

//operator overload for multiplication with a scalar - element wise multiplication
//this function is a friend function of the class so that double * vector3d behaves as vector3d * double
template <class TT>
vector3d<TT> operator* (const TT &scalar, const vector3d<TT> &v1){
  vector3d<TT> temp;
  temp.x = v1.x * scalar;
  temp.y = v1.y * scalar;
  temp.z = v1.z * scalar;
  return temp;
}

//operator overload for divisioin with a scalar - element wise division
template <class T>
vector3d<T> vector3d<T>::operator/ (const T &scalar)const{
  vector3d<T> temp = *this;
  temp.x /= scalar;
  temp.y /= scalar;
  temp.z /= scalar;
  return temp;
}

//operator overload for division with a scalar - element wise division
//this function is a friend function of the class so that double / vector3d works
template <class TT>
vector3d<TT> operator/ (const TT &scalar, const vector3d<TT> &v1){
  vector3d<TT> temp;
  temp.x = scalar/v1.x;
  temp.y = scalar/v1.y;
  temp.z = scalar/v1.z;
  return temp;
}

//operator overload for << - allows use of cout, cerr, etc.
template <class TT>
ostream & operator<< (ostream &os, const vector3d<TT> &v1){
  os << v1.x << ", " << v1.y << ", " << v1.z;
  return os;
}

//Function to calculate the dot product of two vectors
template <class T>
T vector3d<T>::DotProd( const vector3d &v2)const{

  return x*v2.x + y*v2.y + z*v2.z;

}


//Function to calculate the cross product of two vectors
template <class T>
vector3d<T> vector3d<T>::CrossProd(const vector3d &v2)const{

  vector3d<T> crossProd;

  crossProd.x = y*v2.z-z*v2.y;
  crossProd.y = -1.0*(x*v2.z-z*v2.x);
  crossProd.z = x*v2.y-y*v2.x;

  return crossProd;

}

//Function to calculate the magnitude of the vector
template <class T>
T vector3d<T>::Mag()const{
  return sqrt(x*x + y*y + z*z);
}

//Function to calculate the square of the magnitude of the vector
template <class T>
T vector3d<T>::MagSq()const{
  return x*x + y*y + z*z;
}

//Function to sum the elements in the vector
template <class T>
T vector3d<T>::SumElem()const{

  return x + y + z;

}

//Function to calculate the distance between two vector3ds
template <class T>
T vector3d<T>::Distance( const vector3d &v2)const{

  return sqrt( pow(x-v2.x , 2) + pow(y-v2.y , 2) + pow(z-v2.z , 2) );

}



#endif
