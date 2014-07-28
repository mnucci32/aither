#include <cstdlib>      //exit()
#include "matrix.h"
#include <iostream>     //cout

using std::cout;
using std::endl;
using std::cerr;
using std::copy;

//copy constructor
squareMatrix::squareMatrix( const squareMatrix &cp){
  (*this).size = cp.Size();
  (*this).data = new double[cp.Size()*cp.Size()];
  copy(&cp.data[0], &cp.data[0] + cp.Size()*cp.Size(), &(*this).data[0]);
}

//copy assignment operator
squareMatrix& squareMatrix::operator= (squareMatrix other){
  swap(*this, other);
  return *this;
}

//friend function to allow for swap functionality
void swap(squareMatrix &first, squareMatrix &second){
  std::swap(first.size, second.size);
  std::swap(first.data, second.data);
}

//member function to get the data from the matrix
double squareMatrix::Data( const int &r, const int &c )const{

  //test to see that row and column inputs are within bounds
  if ( (r > ((*this).size)) || (c > ((*this).size)) ){
    cerr << "ERROR: The requested data, does not lie within the matrix bounds. Check row and column inputs." << endl;
    exit(1);
  }
  return data[c + r * (*this).size];
}

//member function to set the data in the matrix
void squareMatrix::SetData( const int &r, const int &c, const double &val ){
  //test to see that row and column inputs are within bounds
  if ( (r > ((*this).size)) || (c > ((*this).size)) ){
    cerr << "ERROR: Cannot assign data to given location because it does not lie within the matrix bounds. Check row and column inputs." << endl;
    exit(1);
  }
  data[c + r * (*this).size] = val;
}


//operator overload for addition
squareMatrix squareMatrix::operator + (const squareMatrix& s2)const{
  squareMatrix s1 = *this;

  //check to see that matrix dimensions are the same
  if ( s1.size != s2.size ){
    cerr << "ERROR: Cannot add matricies, dimensions do not agree." << endl;
  }

  int cc = 0;
  int rr = 0;
  for( cc = 0; cc < s2.Size(); cc++ ){
    for( rr = 0; rr < s2.Size(); rr++ ){
      s1.SetData(rr,cc, s1.Data(rr,cc) + s2.Data(rr,cc));
    }
  }
  return s1;
}

//operator overload for addition with a scalar
squareMatrix squareMatrix::operator + (const double &scalar)const{
  squareMatrix s1 = *this;

  int cc = 0;
  int rr = 0;
  for( cc = 0; cc < s1.Size(); cc++ ){
    for( rr = 0; rr < s1.Size(); rr++ ){
      s1.SetData(rr,cc, s1.Data(rr,cc) + scalar);
    }
  }
  return s1;
}

//operator overload for addition with a scalar
squareMatrix operator+ (const double &scalar, const squareMatrix &s2){
  squareMatrix s1(s2.Size());

  int cc = 0;
  int rr = 0;
  for( cc = 0; cc < s2.Size(); cc++ ){
    for( rr = 0; rr < s2.Size(); rr++ ){
      s1.SetData(rr,cc, s2.Data(rr,cc) + scalar);
    }
  }
  return s1;
}


//operator overload for subtraction
squareMatrix squareMatrix::operator - (const squareMatrix& s2)const{
  squareMatrix s1 = *this;

  //check to see that matrix dimensions are the same
  if ( s1.size != s2.size ){
    cerr << "ERROR: Cannot add matricies, dimensions do not agree." << endl;
  }

  int cc = 0;
  int rr = 0;
  for( cc = 0; cc < s2.Size(); cc++ ){
    for( rr = 0; rr < s2.Size(); rr++ ){
      s1.SetData(rr,cc, s1.Data(rr,cc) - s2.Data(rr,cc));
    }
  }
  return s1;
}

//operator overload for subtraction with a scalar
squareMatrix squareMatrix::operator - (const double &scalar)const{
  squareMatrix s1 = *this;

  int cc = 0;
  int rr = 0;
  for( cc = 0; cc < s1.Size(); cc++ ){
    for( rr = 0; rr < s1.Size(); rr++ ){
      s1.SetData(rr,cc, s1.Data(rr,cc) - scalar);
    }
  }
  return s1;
}

//operator overload for subtraction with a scalar
squareMatrix operator- (const double &scalar, const squareMatrix &s2){
  squareMatrix s1(s2.Size());

  int cc = 0;
  int rr = 0;
  for( cc = 0; cc < s2.Size(); cc++ ){
    for( rr = 0; rr < s2.Size(); rr++ ){
      s1.SetData(rr,cc, scalar - s2.Data(rr,cc));
    }
  }
  return s1;
}

//operator overload for multiplication
squareMatrix squareMatrix::operator * (const squareMatrix& s2)const{
  squareMatrix s1 = *this;

  //check to see that matrix dimensions are the same
  if ( s1.size != s2.size ){
    cerr << "ERROR: Cannot multiply matricies, dimensions do not agree." << endl;
  }

  int cc = 0;
  int rr = 0;
  for( cc = 0; cc < s2.Size(); cc++ ){
    for( rr = 0; rr < s2.Size(); rr++ ){
      double newVal = 0.0;
      int ii = 0;
      for( ii = 0; ii < s2.Size(); ii++ ){
	newVal += ((*this).Data(rr,ii) * s2.Data(ii,cc));
      }
      s1.SetData(rr,cc, newVal);
    }
  }
  return s1;
}

//operator overload for multiplication with a scalar
squareMatrix squareMatrix::operator * (const double &scalar)const{
  squareMatrix s1 = *this;

  int cc = 0;
  int rr = 0;
  for( cc = 0; cc < s1.Size(); cc++ ){
    for( rr = 0; rr < s1.Size(); rr++ ){
      s1.SetData(rr,cc, s1.Data(rr,cc) * scalar);
    }
  }
  return s1;
}

//operator overload for multiplication with a scalar
squareMatrix operator* (const double &scalar, const squareMatrix &s2){
  squareMatrix s1(s2.Size());

  int cc = 0;
  int rr = 0;
  for( cc = 0; cc < s2.Size(); cc++ ){
    for( rr = 0; rr < s2.Size(); rr++ ){
      s1.SetData(rr,cc, s2.Data(rr,cc) * scalar);
    }
  }
  return s1;
}

//operator overload for division with a scalar
squareMatrix squareMatrix::operator / (const double &scalar)const{
  squareMatrix s1 = *this;

  int cc = 0;
  int rr = 0;
  for( cc = 0; cc < s1.Size(); cc++ ){
    for( rr = 0; rr < s1.Size(); rr++ ){
      s1.SetData(rr,cc, s1.Data(rr,cc) / scalar);
    }
  }
  return s1;
}

//operator overload for division with a scalar
squareMatrix operator/ (const double &scalar, const squareMatrix &s2){
  squareMatrix s1(s2.Size());

  int cc = 0;
  int rr = 0;
  for( cc = 0; cc < s2.Size(); cc++ ){
    for( rr = 0; rr < s2.Size(); rr++ ){
      s1.SetData(rr,cc, scalar / s2.Data(rr,cc));
    }
  }
  return s1;
}

//operation overload for << - allows use of cout, cerr, etc.
ostream & operator<< (ostream &os, const squareMatrix &m){

  int cc = 0;
  int rr = 0;
  for( rr = 0; rr < m.Size(); rr++ ){
    for( cc = 0; cc < m.Size(); cc++ ){
      cout << m.Data(rr,cc);
      if(cc != (m.Size()-1)){
	cout << ", ";
      }
      else{
	cout << endl;
      }
    }
  }
  return os;
}
