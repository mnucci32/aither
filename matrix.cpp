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
  if ( (r >= ((*this).size)) || (c >= ((*this).size)) ){
    cerr << "ERROR: The requested data, does not lie within the matrix bounds. Check row and column inputs." << endl;
    exit(1);
  }
  return data[c + r * (*this).size];
}

//member function to set the data in the matrix
void squareMatrix::SetData( const int &r, const int &c, const double &val ){
  //test to see that row and column inputs are within bounds
  if ( (r >= ((*this).size)) || (c >= ((*this).size)) ){
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

  I.Identity();

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
      cerr << "ERROR: Singular matrix in Gauss-Jordan elimination! Matrix (mid inversion) is" << endl << *this << endl;
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
    cerr << "ERROR: Cannot subtract matricies, dimensions do not agree." << endl;
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


//member function to zero the matrix
void squareMatrix::Zero(){
  for(int cc = 0; cc < size; cc++){
    for(int rr = 0; rr < size; rr++){
      (*this).SetData(rr,cc,0.0);
    }
  }

}

//member function to set matrix to Identity
void squareMatrix::Identity(){
  for( int rr = 0; rr < (*this).Size(); rr++ ){
    for (int cc = 0; cc < (*this).Size(); cc++ ){
      if(rr == cc){
	(*this).SetData(rr,cc,1.0);
      }
      else{
	(*this).SetData(rr,cc,0.0);
      }
    }
  }

}

colMatrix squareMatrix::Multiply( const colMatrix &X )const{

  //Test to see that column matrix can be multiplied with square matrix
  if ( (*this).Size() != X.Size() ){
    cerr << "ERROR: Column matrix cannot be multiplied with square matrix. Sizes do not agree!" << endl;
    exit(1);
  }

  colMatrix B(X.Size());
  B.Zero();

  for (int rr = 0; rr < X.Size(); rr++ ){
    double tempData = 0.0;
    for ( int cc = 0; cc < X.Size(); cc++ ){
      tempData += (*this).Data(rr,cc) * X.Data(cc);
    }
    B.SetData(rr, tempData);
  }

  return B;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//function definitions for matrixDiagonal class

//copy constructor
matrixDiagonal::matrixDiagonal( const matrixDiagonal &cp){
  (*this).size = cp.Size();
  (*this).data = new squareMatrix[cp.Size()];
  copy(&cp.data[0], &cp.data[0] + cp.Size(), &(*this).data[0]);
}

//copy assignment operator
matrixDiagonal& matrixDiagonal::operator= (matrixDiagonal other){
  swap(*this, other);
  return *this;
}

//friend function to allow for swap functionality
void swap(matrixDiagonal &first, matrixDiagonal &second){
  std::swap(first.size, second.size);
  std::swap(first.data, second.data);
}

//member function to get the data from the matrix
squareMatrix matrixDiagonal::Data( const int &ind)const{

  //test to see that the index input is within bounds
  if ( ind >= (*this).size ){
    cerr << "ERROR: The requested data, does not lie within the matrix bounds. Check index input." << endl;
    exit(1);
  }
  return data[ind];
}

//member function to set the data in the matrix
void matrixDiagonal::SetData( const int &ind, const squareMatrix &val ){
  //test to see that the index input is within bounds
  if ( ind >= (*this).size ){
    cerr << "ERROR: The requested data, does not lie within the matrix bounds. Check index input." << endl;
    exit(1);
  }
  data[ind] = val;
}

//operation overload for << - allows use of cout, cerr, etc.
ostream & operator<< (ostream &os, const matrixDiagonal &m){

  int rr = 0;
  for( rr = 0; rr < m.Size(); rr++ ){
    cout << "In index " << rr << " matrix block is:" << endl;
    cout << m.Data(rr) << endl;
  }
  return os;
}


//member function to zero the matrix
void matrixDiagonal::Zero(const int &s){
  squareMatrix mZero(s);
  mZero.Zero();

  for(int cc = 0; cc < size; cc++){
    (*this).SetData(cc,mZero);
  }

}


//member function to delete the contents of the data structure and resize it
void matrixDiagonal::CleanResizeZero(const int &s, const int &m){
  squareMatrix mZero(m);
  mZero.Zero();

  delete [] (*this).data;
  (*this).data = new squareMatrix[s];
  (*this).size = s;

  for(int cc = 0; cc < size; cc++){
    (*this).SetData(cc,mZero);
  }


}

//member function to invert each matrix stored
void matrixDiagonal::Inverse(){
  for(int ii = 0; ii < (*this).Size(); ii++ ){
    squareMatrix temp = (*this).Data(ii);
    temp.Inverse();
    (*this).SetData(ii,temp);
  }
}

//functions ------------------------------------------------------------------------------------------------
//function to perform symmetric Gauss-Seidel relaxation to solver Ax=b
//when relax = 1.0, symmetric Gauss-Seidel is achieved. Values >1 result in symmetric successive over relaxation (SSOR)
//Values <1 result in under relaxation
double LUSGS( const colMatrix &Aii, const matrixDiagonal &Ail, const matrixDiagonal &Aiu, const matrixDiagonal &Ajl, const matrixDiagonal &Aju, const matrixDiagonal &Akl, const matrixDiagonal &Aku, vector<colMatrix> &x, const vector<colMatrix> &b, const vector<colMatrix> &solTimeMmN, const vector<colMatrix> &solDeltaNm1, const int &sweeps, const double &relax, const int &imax, const int &jmax, const double &theta){

  //Aii --> block matrix of the main diagonal
  //Ail --> block matrix of the lower i diagonal
  //Aiu --> block matrix of the upper i diagonal
  //Ajl --> block matrix of the lower j diagonal
  //Aju --> block matrix of the upper j diagonal
  //Akl --> block matrix of the lower k diagonal
  //Aku --> block matrix of the upper k diagonal
  //x   --> block vector of correction
  //b   --> block vector of residuals
  //sweeps --> number of symmetric sweeps to perform
  //relax  --> relaxation parameter >1 is overrelaxation, <1 is underrelaxation
  //imax   --> imax for block
  //jmax   --> jmax for block

  //initialize x to 0
  for (unsigned int ll = 0; ll < x.size(); ll++ ){
    x[ll].Zero();
  }

  //invert main diagonal
  // matrixDiagonal AiiInv = Aii;
  // AiiInv.Inverse();
  // squareMatrix AiiInv(x[0].Size());
  double AiiInv = 0.0;


  colMatrix newData(x[0].Size());
  colMatrix oldData(x[0].Size());

  colMatrix l2Resid(x[0].Size());

  double thetaInv = 1.0 / theta;

  for ( int kk = 0; kk < sweeps; kk++ ){
    //forward sweep
    for ( int ii = 0; ii < (int)x.size(); ii++ ){

      int il = GetDiagPosLowerI(ii);
      int iu = GetDiagPosUpperI(ii);
      int jl = GetDiagPosLowerJ(ii,imax);
      int ju = GetDiagPosUpperJ(ii,imax);
      int kl = GetDiagPosLowerK(ii,imax,jmax);
      int ku = GetDiagPosUpperK(ii,imax,jmax);

      newData.Zero();
      oldData.Zero();

      if ( il >=0 && il < (int)x.size() ){
	oldData = oldData + Ail.Data(ii).Multiply(x[il]);
      }
      if ( jl >=0 && jl < (int)x.size() ){
	oldData = oldData + Ajl.Data(ii).Multiply(x[jl]);
      }
      if ( kl >=0 && kl < (int)x.size() ){
	oldData = oldData + Akl.Data(ii).Multiply(x[kl]);
      }

      if ( iu >=0 && iu < (int)x.size() ){
	newData = newData + Aiu.Data(ii).Multiply(x[iu]);
      }
      if ( ju >=0 && ju < (int)x.size() ){
	newData = newData + Aju.Data(ii).Multiply(x[ju]);
      }
      if ( ku >=0 && ku < (int)x.size() ){
	newData = newData + Aku.Data(ii).Multiply(x[ku]);
      }

      AiiInv = 1.0 / Aii.Data(ii);
      // x[ii] = AiiInv * ( b[ii] - newData - oldData) ;


      x[ii] = (1.0 - relax) * x[ii] + relax * AiiInv * ( thetaInv * b[ii] + solDeltaNm1[ii] +
      	      solTimeMmN[ii] - newData - oldData) ;

    }

    //backward sweep
    for ( int ii = (int)x.size()-1; ii >= 0; ii-- ){

      int il = GetDiagPosLowerI(ii);
      int iu = GetDiagPosUpperI(ii);
      int jl = GetDiagPosLowerJ(ii,imax);
      int ju = GetDiagPosUpperJ(ii,imax);
      int kl = GetDiagPosLowerK(ii,imax,jmax);
      int ku = GetDiagPosUpperK(ii,imax,jmax);

      newData.Zero();
      oldData.Zero();

      if ( iu >=0 && iu < (int)x.size() ){
	oldData = oldData + Aiu.Data(ii).Multiply(x[iu]);
      }
      if ( ju >=0 && ju < (int)x.size() ){
	oldData = oldData + Aju.Data(ii).Multiply(x[ju]);
      }
      if ( ku >=0 && ku < (int)x.size() ){
	oldData = oldData + Aku.Data(ii).Multiply(x[ku]);
      }

      if ( il >=0 && il < (int)x.size() ){
	newData = newData + Ail.Data(ii).Multiply(x[il]);
      }
      if ( jl >=0 && jl < (int)x.size() ){
	newData = newData + Ajl.Data(ii).Multiply(x[jl]);
      }
      if ( kl >=0 && kl < (int)x.size() ){
	newData = newData + Akl.Data(ii).Multiply(x[kl]);
      }

      AiiInv = 1.0 / Aii.Data(ii);
      // x[ii] = AiiInv * ( b[ii] - newData - oldData) ;

      x[ii] = (1.0 - relax) * x[ii] + relax * AiiInv * ( thetaInv * b[ii] + solDeltaNm1[ii] +
              solTimeMmN[ii] - newData - oldData) ;

    }

    //calculate residual
    l2Resid.Zero();
    colMatrix resid(x[0].Size());

    for ( int ii = 0; ii < (int)x.size(); ii++ ){

      int il = ii-1;
      int iu = ii+1;
      int jl = ii-imax;
      int ju = ii+imax;
      int kl = ii-imax*jmax;
      int ku = ii+imax*jmax;

      // resid = b[ii] - Aii.Data(ii) * x[ii];

      resid = thetaInv * b[ii] + solDeltaNm1[ii] + solTimeMmN[ii] - Aii.Data(ii) * x[ii];

      if ( il >=0 && il < (int)x.size() ){
	resid = resid - Ail.Data(ii).Multiply(x[il]);
      }

      if ( iu >=0 && iu < (int)x.size() ){
	resid = resid - Aiu.Data(ii).Multiply(x[iu]);
      }

      if ( jl >=0 && jl < (int)x.size() ){
	resid = resid - Ajl.Data(ii).Multiply(x[jl]);
      }

      if ( ju >=0 && ju < (int)x.size() ){
	resid = resid - Aju.Data(ii).Multiply(x[ju]);
      }

      if ( kl >=0 && kl < (int)x.size() ){
	resid = resid - Akl.Data(ii).Multiply(x[kl]);
      }

      if ( ku >=0 && ku < (int)x.size() ){
	resid = resid - Aku.Data(ii).Multiply(x[ku]);
      }

      l2Resid = l2Resid + resid * resid;
    }


  } //loop for sweeps

  return l2Resid.Sum();

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////
//functions for colMatrix class
//
//copy constructor
colMatrix::colMatrix( const colMatrix &cp){
  (*this).size = cp.Size();
  (*this).data = new double[cp.Size()];
  copy(&cp.data[0], &cp.data[0] + cp.Size(), &(*this).data[0]);
}

//copy assignment operator
colMatrix& colMatrix::operator= (colMatrix other){
  swap(*this, other);
  return *this;
}

//friend function to allow for swap functionality
void swap(colMatrix &first, colMatrix &second){
  std::swap(first.size, second.size);
  std::swap(first.data, second.data);
}

//member function to get the data from the matrix
double colMatrix::Data( const int &r)const{

  //test to see that row and column inputs are within bounds
  if ( r >= (*this).size ){
    cerr << "ERROR: The requested data, does not lie within the column matrix bounds. Check row input." << endl;
    exit(1);
  }
  return data[r];
}

//member function to set the data in the matrix
void colMatrix::SetData( const int &r, const double &val ){
  //test to see that row and column inputs are within bounds
  if ( r >= (*this).size ){
    cerr << "ERROR: Cannot assign data to given location because it does not lie within the column matrix bounds. Check row input." << endl;
    exit(1);
  }
  data[r] = val;
}

//operator overload for addition
colMatrix colMatrix::operator + (const colMatrix& s2)const{
  colMatrix s1 = *this;

  //check to see that matrix dimensions are the same
  if ( s1.size != s2.size ){
    cerr << "ERROR: Cannot add column matricies, dimensions do not agree." << endl;
  }

  for( int rr = 0; rr < s2.Size(); rr++ ){
    s1.SetData(rr, s1.Data(rr) + s2.Data(rr));
  }
  return s1;
}

//operator overload for addition
colMatrix colMatrix::operator + (const vector<double>& v1)const{
  colMatrix s1 = *this;

  //check to see that matrix dimensions are the same
  if ( s1.size != (int)v1.size() ){
    cerr << "ERROR: Cannot add column matrix of size " << s1.size << " to vector class of size " << v1.size()<< ", dimensions do not agree." << endl;
  }

  for( int rr = 0; rr < s1.size; rr++ ){
    s1.SetData(rr, s1.Data(rr) + v1[rr]);
  }
  return s1;
}


//operator overload for addition with a scalar
colMatrix colMatrix::operator + (const double &scalar)const{
  colMatrix s1 = *this;

  for( int rr = 0; rr < s1.Size(); rr++ ){
    s1.SetData(rr, s1.Data(rr) + scalar);
  }
  return s1;
}

//operator overload for addition with a scalar
colMatrix operator+ (const double &scalar, const colMatrix &s2){
  colMatrix s1(s2.Size());

  for( int rr = 0; rr < s2.Size(); rr++ ){
    s1.SetData(rr, s2.Data(rr) + scalar);
  }
  return s1;
}


//operator overload for subtraction
colMatrix colMatrix::operator - (const colMatrix& s2)const{
  colMatrix s1 = *this;

  //check to see that matrix dimensions are the same
  if ( s1.size != s2.size ){
    cerr << "ERROR: Cannot subtract column matricies, dimensions do not agree." << endl;
  }

  for( int rr = 0; rr < s2.Size(); rr++ ){
    s1.SetData(rr, s1.Data(rr) - s2.Data(rr));
  }
  return s1;
}

//operator overload for addition
colMatrix colMatrix::operator - (const vector<double>& v1)const{
  colMatrix s1 = *this;

  //check to see that matrix dimensions are the same
  if ( s1.size != (int)v1.size() ){
    cerr << "ERROR: Cannot subtract column matrix of size " << s1.size << " to vector class of size " << v1.size()<< ", dimensions do not agree." << endl;
  }

  for( int rr = 0; rr < s1.size; rr++ ){
    s1.SetData(rr, s1.Data(rr) - v1[rr]);
  }
  return s1;
}

//operator overload for subtraction with a scalar
colMatrix colMatrix::operator - (const double &scalar)const{
  colMatrix s1 = *this;

  for( int rr = 0; rr < s1.Size(); rr++ ){
    s1.SetData(rr, s1.Data(rr) - scalar);
  }
  return s1;
}

//operator overload for subtraction with a scalar
colMatrix operator- (const double &scalar, const colMatrix &s2){
  colMatrix s1(s2.Size());

  for( int rr = 0; rr < s2.Size(); rr++ ){
    s1.SetData(rr, scalar - s2.Data(rr));
  }
  return s1;
}

//operator overload for elementwise multiplication
colMatrix colMatrix::operator * (const colMatrix& s2)const{
  colMatrix s1 = *this;

  //check to see that matrix dimensions are the same
  if ( s1.size != s2.size ){
    cerr << "ERROR: Cannot elementwise multiply column matricies, dimensions do not agree." << endl;
  }

  for( int rr = 0; rr < s2.Size(); rr++ ){
    s1.SetData(rr, s1.Data(rr) * s2.Data(rr));
  }
  return s1;
}

//operator overload for multiplication with a scalar
colMatrix colMatrix::operator * (const double &scalar)const{
  colMatrix s1 = *this;

  for( int rr = 0; rr < s1.Size(); rr++ ){
    s1.SetData(rr, s1.Data(rr) * scalar);
  }
  return s1;
}

//operator overload for multiplication with a scalar
colMatrix operator* (const double &scalar, const colMatrix &s2){
  colMatrix s1(s2.Size());

  for( int rr = 0; rr < s2.Size(); rr++ ){
    s1.SetData(rr, s2.Data(rr) * scalar);
  }
  return s1;
}

//operator overload for elementwise division
colMatrix colMatrix::operator / (const colMatrix& s2)const{
  colMatrix s1 = *this;

  //check to see that matrix dimensions are the same
  if ( s1.size != s2.size ){
    cerr << "ERROR: Cannot elementwise divide column matricies, dimensions do not agree." << endl;
  }

  for( int rr = 0; rr < s2.Size(); rr++ ){
    s1.SetData(rr, s1.Data(rr) / s2.Data(rr));
  }
  return s1;
}

//operator overload for division with a scalar
colMatrix colMatrix::operator / (const double &scalar)const{
  colMatrix s1 = *this;

  for( int rr = 0; rr < s1.Size(); rr++ ){
    s1.SetData(rr, s1.Data(rr) / scalar);
  }
  return s1;
}

//operator overload for division with a scalar
colMatrix operator/ (const double &scalar, const colMatrix &s2){
  colMatrix s1(s2.Size());

  for( int rr = 0; rr < s2.Size(); rr++ ){
    s1.SetData(rr, scalar / s2.Data(rr));
  }
  return s1;
}

//operation overload for << - allows use of cout, cerr, etc.
ostream & operator<< (ostream &os, const colMatrix &m){

  for( int rr = 0; rr < m.Size(); rr++ ){
    cout << m.Data(rr) << endl;
  }

  return os;
}


//member function to zero the matrix
void colMatrix::Zero(){
  for(int rr = 0; rr < size; rr++){
    (*this).SetData(rr,0.0);
  }
}

//member function to sum column matrix
double colMatrix::Sum(){
  double sum = 0.0;
  for( int ii = 0; ii < (*this).Size(); ii++ ){
    sum += (*this).Data(ii);
  }
  return sum;
}

//member function to delete the contents of the data structure and resize it
void colMatrix::CleanResizeZero(const int &s){

  delete [] (*this).data;
  (*this).data = new double[s];
  (*this).size = s;

  for(int cc = 0; cc < size; cc++){
    (*this).SetData(cc,0.0);
  }


}
