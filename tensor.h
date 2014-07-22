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
  T xx,xy,xz,yx,yy,yz,zx,zy,zz;
 public:

  //constructor
  tensor( T a, T b, T c, T d, T e, T f, T g, T h, T i) : xx(a), xy(b), xz(c), yx(d), yy(e), yz(f), zx(g), zy(h), zz(i) {}
  tensor() : xx(0), xy(0), xz(0), yx(0), yy(0), yz(0), zx(0), zy(0), zz(0) {}
  tensor(T i) : xx(i), xy(0), xz(0), yx(0), yy(i), yz(0), zx(0), zy(0), zz(i) {} 
  tensor(vector3d<T> v1, vector3d<T> v2, vector3d<T> v3) : xx(v1.X()), yx(v1.Y()), zx(v1.Z()), xy(v2.X()), yy(v2.Y()), zy(v2.Z()), xz(v3.X()), yz(v3.Y()), zz(v3.Z()){}

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
  void SetXX(const T& val){xx = val;}
  void SetXY(const T& val){xy = val;}
  void SetXZ(const T& val){xz = val;}
  void SetYX(const T& val){yx = val;}
  void SetYY(const T& val){yy = val;}
  void SetYZ(const T& val){yz = val;}
  void SetZX(const T& val){zx = val;}
  void SetZY(const T& val){zy = val;}
  void SetZZ(const T& val){zz = val;}
  //access of data members
  T XX()const{return xx;}
  T XY()const{return xy;}
  T XZ()const{return xz;}
  T YX()const{return yx;}
  T YY()const{return yy;}
  T YZ()const{return yz;}
  T ZX()const{return zx;}
  T ZY()const{return zy;}
  T ZZ()const{return zz;}
  //math functions
  T Trace()const{return xx+yy+zz;}
  tensor<T> Transpose()const;
  vector3d<T> MatMult(const vector3d<T>&)const;
  tensor<T> Identity()const;


  //destructor
  ~tensor() {}

};

//constructor to place 3 vector3ds in a tensor
/* template <class T> */
/* tensor<T>::tensor(vector3d<T> v1, vector3d<T> v2, vector3d<T> v3){ */
/*   xx = v1.X(); */
/*   yx = v1.Y(); */
/*   zx = v1.Z(); */

/*   xy = v2.X(); */
/*   yy = v2.Y(); */
/*   zy = v2.Z(); */

/*   xz = v3.X(); */
/*   yz = v3.Y(); */
/*   zz = v3.Z(); */
/* } */

//operator overload for addition - element wise addition
template <class T>
tensor<T> tensor<T>::operator+ (const tensor &v2)const{
  tensor<T> temp = *this;
  temp.xx += v2.xx;
  temp.xy += v2.xy;
  temp.xz += v2.xz;

  temp.yx += v2.yx;
  temp.yy += v2.yy;
  temp.yz += v2.yz;

  temp.zx += v2.zx;
  temp.zy += v2.zy;
  temp.zz += v2.zz;

  return temp;
}

//operator overload for subtraction - element wise subtraction
template <class T>
tensor<T> tensor<T>::operator- (const tensor &v2)const{
  tensor<T> temp = *this;
  temp.xx -= v2.xx;
  temp.xy -= v2.xy;
  temp.xz -= v2.xz;

  temp.yx -= v2.yx;
  temp.yy -= v2.yy;
  temp.yz -= v2.yz;

  temp.zx -= v2.zx;
  temp.zy -= v2.zy;
  temp.zz -= v2.zz;

  return temp;
}

//operator overload for multiplication - matrix multiplication
template <class T>
tensor<T> tensor<T>::operator* (const tensor &v2)const{
  tensor<T> temp;

  temp.xx = (*this).xx * v2.xx + (*this).yx * v2.xy + (*this).zx * v2.xz;
  temp.xy = (*this).xy * v2.xx + (*this).yy * v2.xy + (*this).zy * v2.xz;
  temp.xz = (*this).xz * v2.xx + (*this).yz * v2.xy + (*this).zz * v2.xz;

  temp.yx = (*this).xx * v2.yx + (*this).yx * v2.yy + (*this).zx * v2.yz;
  temp.yy = (*this).xy * v2.yx + (*this).yy * v2.yy + (*this).zy * v2.yz;
  temp.yz = (*this).xz * v2.yx + (*this).yz * v2.yy + (*this).zz * v2.yz;

  temp.zx = (*this).xx * v2.zx + (*this).yx * v2.zy + (*this).zx * v2.zz;
  temp.zy = (*this).xy * v2.zx + (*this).yy * v2.zy + (*this).zy * v2.zz;
  temp.zz = (*this).xz * v2.zx + (*this).yz * v2.zy + (*this).zz * v2.zz;

  return temp;
}


//operator overload for multiplication with a scalar - element wise multiplication
template <class T>
tensor<T> tensor<T>::operator* (const T &scalar)const{
  tensor<T> temp = *this;
  temp.xx *= scalar;
  temp.xy *= scalar;
  temp.xz *= scalar;

  temp.yx *= scalar;
  temp.yy *= scalar;
  temp.yz *= scalar;

  temp.zx *= scalar;
  temp.zy *= scalar;
  temp.zz *= scalar;

  return temp;
}

//operator overload for multiplication with a scalar - element wise multiplication
//this function is a friend function of the class so that double * tensor behaves as tensor * double
template <class TT>
tensor<TT> operator* (const TT &scalar, const tensor<TT> &v1){
  tensor<TT> temp;
  temp.xx = v1.xx * scalar;
  temp.xy = v1.xy * scalar;
  temp.xz = v1.xz * scalar;

  temp.yx = v1.yx * scalar;
  temp.yy = v1.yy * scalar;
  temp.yz = v1.yz * scalar;

  temp.zx = v1.zx * scalar;
  temp.zy = v1.zy * scalar;
  temp.zz = v1.zz * scalar;

  return temp;
}

//operator overload for division with a scalar - element wise division
template <class T>
tensor<T> tensor<T>::operator/ (const T &scalar)const{
  tensor<T> temp = *this;
  temp.xx /= scalar;
  temp.xy /= scalar;
  temp.xz /= scalar;

  temp.yx /= scalar;
  temp.yy /= scalar;
  temp.yz /= scalar;

  temp.zx /= scalar;
  temp.zy /= scalar;
  temp.zz /= scalar;

  return temp;
}

//operator overload for division with a scalar - element wise division
//this function is a friend function of the class so that double / tensor works
template <class TT>
tensor<TT> operator/ (const TT &scalar, const tensor<TT> &v1){
  tensor<TT> temp;
  temp.xx = scalar/v1.xx;
  temp.xy = scalar/v1.xy;
  temp.xz = scalar/v1.xz;

  temp.yx = scalar/v1.yx;
  temp.yy = scalar/v1.yy;
  temp.yz = scalar/v1.yz;

  temp.zx = scalar/v1.zx;
  temp.zy = scalar/v1.zy;
  temp.zz = scalar/v1.zz;

  return temp;
}

//operator overload for << - allows use of cout, cerr, etc.
template <class TT>
ostream & operator<< (ostream &os, const tensor<TT> &v1){
  os << v1.xx << ", " << v1.yx << ", " << v1.zx << endl;
  os << v1.xy << ", " << v1.yy << ", " << v1.zy << endl;
  os << v1.xz << ", " << v1.yz << ", " << v1.zz << endl;
  return os;
}

//Function to return the transpose of the given tensor
template <class T>
tensor<T> tensor<T>::Transpose()const{

  tensor<T> temp;

  temp.xx = (*this).xx;
  temp.xy = (*this).yx;
  temp.xz = (*this).zx;

  temp.yx = (*this).xy;
  temp.yy = (*this).yy;
  temp.yz = (*this).zy;

  temp.zx = (*this).xz;
  temp.zy = (*this).yz;
  temp.zz = (*this).zz;

  return temp;

}

//Function to return the matrix multiplication of a tensor and vector3d
template <class T>
vector3d<T> tensor<T>::MatMult(const vector3d<T> &vec)const{

  vector3d<T> temp;

  temp.SetX( (*this).xx * vec.X() + (*this).yx * vec.Y() + (*this).zx * vec.Z() );
  temp.SetY( (*this).xy * vec.X() + (*this).yy * vec.Y() + (*this).zy * vec.Z() );
  temp.SetZ( (*this).xz * vec.X() + (*this).yz * vec.Y() + (*this).zz * vec.Z() );

  return temp;

}

//Function to return the identity matrix
template <class T>
tensor<T> tensor<T>::Identity()const{

  tensor<T> temp;
  T one = 1;
  T zero = 0;

  temp.xx = one;
  temp.xy = zero;
  temp.xz = zero;

  temp.yx = zero;
  temp.yy = one;
  temp.yz = zero;

  temp.zx = zero;
  temp.zy = zero;
  temp.zz = one;

  return temp;

}


#endif
