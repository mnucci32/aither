/*  An open source Navier-Stokes CFD solver.
    Copyright (C) 2015  Michael Nucci (michael.nucci@gmail.com)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>. */

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
  T data[3];
 public:

  //constructor
 vector3d( T a, T b, T c) : data{a, b, c} {}
 vector3d() : data{0, 0, 0} {}

  //member functions
  //operator overloads
  vector3d<T> operator + (const vector3d&)const;
  vector3d<T> operator - (const vector3d&)const;
  vector3d<T> operator * (const T&)const;
  vector3d<T> operator / (const T&)const;
  T operator[] (const int &a)const{return data[a];}
  T& operator[] (const int &a){return data[a];}
  bool operator == (const vector3d&)const;
  template <class TT>
  friend vector3d<TT> operator * (const TT&, const vector3d<TT>&);
  template <class TT>
  friend vector3d<TT> operator / (const TT&, const vector3d<TT>&);
  template <class TT>
  friend ostream & operator<< (ostream &os, const vector3d<TT>&);
  //assignment of data members
  void SetX(const T& val){data[0] = val;}
  void SetY(const T& val){data[1] = val;}
  void SetZ(const T& val){data[2] = val;}
  //access of data members
  T X()const{return data[0];}
  T Y()const{return data[1];}
  T Z()const{return data[2];}
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
  temp.data[0] += v2.data[0];
  temp.data[1] += v2.data[1];
  temp.data[2] += v2.data[2];
  return temp;
}

//operator overload for subtraction - element wise subtraction
template <class T>
vector3d<T> vector3d<T>::operator- (const vector3d &v2)const{
  vector3d<T> temp = *this;
  temp.data[0] -= v2.data[0];
  temp.data[1] -= v2.data[1];
  temp.data[2] -= v2.data[2];
  return temp;
}

//operator overload for multiplication with a scalar - element wise multiplication
template <class T>
vector3d<T> vector3d<T>::operator* (const T &scalar)const{
  vector3d<T> temp = *this;
  temp.data[0] *= scalar;
  temp.data[1] *= scalar;
  temp.data[2] *= scalar;
  return temp;
}

//operator overload for multiplication with a scalar - element wise multiplication
//this function is a friend function of the class so that double * vector3d behaves as vector3d * double
template <class TT>
vector3d<TT> operator* (const TT &scalar, const vector3d<TT> &v1){
  vector3d<TT> temp;
  temp.data[0] = v1.data[0] * scalar;
  temp.data[1] = v1.data[1] * scalar;
  temp.data[2] = v1.data[2] * scalar;
  return temp;
}

//operator overload for divisioin with a scalar - element wise division
template <class T>
vector3d<T> vector3d<T>::operator/ (const T &scalar)const{
  vector3d<T> temp = *this;
  temp.data[0] /= scalar;
  temp.data[1] /= scalar;
  temp.data[2] /= scalar;
  return temp;
}

//operator overload for division with a scalar - element wise division
//this function is a friend function of the class so that double / vector3d works
template <class TT>
vector3d<TT> operator/ (const TT &scalar, const vector3d<TT> &v1){
  vector3d<TT> temp;
  temp.data[0] = scalar/v1.data[0];
  temp.data[1] = scalar/v1.data[1];
  temp.data[2] = scalar/v1.data[2];
  return temp;
}

//operator overload for << - allows use of cout, cerr, etc.
template <class TT>
ostream & operator<< (ostream &os, const vector3d<TT> &v1){
  os << v1.data[0] << ", " << v1.data[1] << ", " << v1.data[2];
  return os;
}

//Function to calculate the dot product of two vectors
template <class T>
T vector3d<T>::DotProd( const vector3d &v2)const{

  return data[0]*v2.data[0] + data[1]*v2.data[1] + data[2]*v2.data[2];

}

//operator overload for comparison 
template <class T>
bool vector3d<T>::operator== (const vector3d &v2)const{
  bool test = false;
  vector3d<T> temp = *this;
  if ( temp.data[0] == v2.data[0] && temp.data[1] == v2.data[1] && temp.data[2] == v2.data[2] ){
    test = true;
  }
  return test;
}

//Function to calculate the cross product of two vectors
template <class T>
vector3d<T> vector3d<T>::CrossProd(const vector3d &v2)const{

  vector3d<T> crossProd;

  crossProd.data[0] = data[1]*v2.data[2]-data[2]*v2.data[1];
  crossProd.data[1] = -1.0*(data[0]*v2.data[2]-data[2]*v2.data[0]);
  crossProd.data[2] = data[0]*v2.data[1]-data[1]*v2.data[0];

  return crossProd;

}

//Function to calculate the magnitude of the vector
template <class T>
T vector3d<T>::Mag()const{
  return sqrt(data[0]*data[0] + data[1]*data[1] + data[2]*data[2]);
}

//Function to calculate the square of the magnitude of the vector
template <class T>
T vector3d<T>::MagSq()const{
  return data[0]*data[0] + data[1]*data[1] + data[2]*data[2];
}

//Function to sum the elements in the vector
template <class T>
T vector3d<T>::SumElem()const{

  return data[0] + data[1] + data[2];

}

//Function to calculate the distance between two vector3ds
template <class T>
T vector3d<T>::Distance( const vector3d &v2)const{

  return sqrt( pow(data[0]-v2.data[0] , 2) + pow(data[1]-v2.data[1] , 2) + pow(data[2]-v2.data[2] , 2) );

}



#endif
