#include <cstdlib>      //exit()
#include "matrix.h"
#include <iostream>     //cout

using std::cout;
using std::endl;
using std::cerr;

//member function to get the data from the matrix
double squareMatrix::Data( const int &r, const int &c )const{

  //test to see that row and column inputs are within bounds
  if ( (r > ((*this).size)) || (c > ((*this).size)) ){
    cerr << "ERROR: The requested data, does not lie within the matrix bounds. Check row and column inputs." << endl;
    exit(1);
  }

  return data[(c-1) + (r-1) * (*this).size];

}

//member function to set the data in the matrix
void squareMatrix::SetData( const int &r, const int &c, const double &val ){
  //test to see that row and column inputs are within bounds
  if ( (r > ((*this).size)) || (c > ((*this).size)) ){
    cerr << "ERROR: Cannot assign data to given location because it does not lie within the matrix bounds. Check row and column inputs." << endl;
    exit(1);
  }

  data[(c-1) + (r-1) * (*this).size] = val;

}


//operator overload for addition
squareMatrix squareMatrix::operator + (const squareMatrix& s2)const{
  squareMatrix s1 = *this;

  //check to see that matrix dimensions are the same
  if ( s1.size != s2.size ){
    cerr << "ERROR: Cannot add matricies, dimensions do not agree." << endl;
  }

  int ii = 0;
  for ( ii = 0; ii < (s1.size * s1.size); ii++ ){
    s1.data[ii] += s2.data[ii];
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

  int ii = 0;
  for ( ii = 0; ii < (s1.size * s1.size); ii++ ){
    s1.data[ii] -= s2.data[ii];
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
      s1.SetData(rr,cc, s2.Data(rr,cc) - scalar);
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
	newVal += s1.Data(rr,ii) * s2.Data(ii,cc);
      }
      s1.SetData(rr,cc, newVal);
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
  for( cc = 0; cc < m.Size(); cc++ ){
    for( rr = 0; rr < m.Size(); rr++ ){
      cout << m.Data(rr,cc);
      if(rr != (m.Size()-1)){
	cout << ", ";
      }
      else{
	cout << endl;
      }
    }
  }

  return os;
}
