#include <cstdlib>      //exit()
#include "matrix.h"
#include <iostream>     //cout
#include <cmath>

using std::cout;
using std::endl;
using std::cerr;
using std::copy;
using std::swap_ranges;
using std::fabs;

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

//member function to swap rows of matrix
void squareMatrix::SwapRows(const int &r1, const int &r2){
  if(r1 != r2){
    for( int ii = 0; ii < size; ii++ ){
      int ind1 = ii + r1*size;
      int ind2 = ii + r2*size;
      std::swap((*this).data[ind1], (*this).data[ind2]);
    }
  }
}

//member function to invert matrix using Gauss-Jordan elimination
void squareMatrix::Inverse(){

  squareMatrix I((*this).size);
  int r = 0;
  int c = 0;

  //initialize identity matrix
  for( r = 0; r < I.Size(); r++ ){
    for( c = 0; c < I.Size(); c++ ){
      if(r == c){
	I.SetData(r,c,1.0);
      }
      else{
	I.SetData(r,c,0.0);
      }
    }
  }

  int cPivot = 0;
  for( cPivot = 0, r = 0; r < size; r++, cPivot++ ){

    //find pivot row
    int rPivot = (*this).FindMaxInCol(r,cPivot,size-1);

    //swap rows
    (*this).SwapRows(r,rPivot);
    I.SwapRows(r,rPivot);

    if(r != 0){  //if not on first row, need to get rid entries ahead of pivot
      for( int ii = 0; ii < cPivot; ii++ ){
	double factor = (*this).Data(r,ii) / (*this).Data(ii,ii);
	(*this).LinCombRow(ii, factor, r);
	I.LinCombRow(ii, factor, r);
      }
    }

    //normalize row by pivot
    if((*this).Data(r,cPivot) == 0.0){
      cerr << "ERROR: Singular matrix in Gauss-Jordan elimination!" << endl;
      exit(1);
    }
    double normFactor = 1.0/(*this).Data(r,cPivot);
    (*this).RowMultiply(r,cPivot,normFactor);  //only multiply entries from pivot and to the right
    I.RowMultiply(r,0,normFactor); //multiply all entries

  }

  //matrix is now upper triangular, work way back up to identity matrix
  cPivot = size - 2;  //start with second to last row
  for( cPivot = size - 2, r = size - 2; r >= 0; r--, cPivot-- ){
    for( int ii = size - 1; ii > cPivot; ii-- ){
      double factor = (*this).Data(r,ii);
      (*this).LinCombRow(ii, factor, r);
      I.LinCombRow(ii, factor, r);
    }

  }

  //set this matrix equal to its inverse
  (*this) = I;
}

//member function to add a linear combination of one row to another
void squareMatrix::LinCombRow(const int &r1, const double &factor, const int &r2){
  for( int ii = 0; ii < size; ii++ ){
    (*this).SetData(r2, ii, (*this).Data(r2,ii) - (*this).Data(r1,ii) * factor);
  }
}

//member function to multiply a row by a given factor
void squareMatrix::RowMultiply(const int &r, const int &c, const double &factor){
  for( int ii = c; ii < size; ii++ ){
    (*this).SetData(r, ii, (*this).Data(r,ii) * factor);
  }
}

//member function to find maximum absolute value in a given column and range within that column and return the corresponding row indice
int squareMatrix::FindMaxInCol(const int &c, const int &start, const int &end)const{
  double maxVal = 0.0;
  int maxRow = 0;
  for( int ii = start; ii < end+1; ii++ ){
    if( fabs((*this).Data(ii,c)) > maxVal ){
      maxVal = fabs((*this).Data(ii,c));
      maxRow = ii;
    }
  }
  return maxRow;
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
