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

#ifndef MATRIXHEADERDEF             //only if the macro MATRIXHEADERDEF is not defined execute these lines of code
#define MATRIXHEADERDEF             //define the macro

#include <vector>  //vector
#include <string>  //string
#include <math.h>       //sqrt
#include <iostream>
#include "plot3d.hpp" //matrix location functions
#include "mpi.h" //parallelism
#include "macros.hpp"

using std::vector;
using std::string;
using std::ostream;

//class to store a column matrix
class colMatrix {
  int size;
  double *data;

 public:
  //constructor
  colMatrix( const int &a ): size(a) {data = new double[a];}
  colMatrix() : size(0), data(NULL){}

  //copy constructor
  colMatrix( const colMatrix &cp);

  //copy assignment operator
  colMatrix& operator= (colMatrix other);

  //move constructor
  colMatrix(colMatrix &&other) : colMatrix() {
    swap(*this,other);
    other.data=NULL;
  }

  //member functions
  double Data(const int &)const;
  void SetData(const int &, const double&);
  int Size()const{return size;}
  void Zero();
  double Sum();
  void CleanResizeZero(const int &);

  //operator overloads
  colMatrix operator + (const colMatrix&)const;
  colMatrix operator - (const colMatrix&)const;
  colMatrix operator * (const colMatrix&)const;
  colMatrix operator / (const colMatrix&)const;

  colMatrix operator + (const vector<double>&)const;
  colMatrix operator - (const vector<double>&)const;

  colMatrix operator + (const double&)const;
  colMatrix operator - (const double&)const;
  colMatrix operator * (const double&)const;
  colMatrix operator / (const double&)const;

  friend colMatrix operator + (const double&, const colMatrix&);
  friend colMatrix operator - (const double&, const colMatrix&);
  friend colMatrix operator * (const double&, const colMatrix&);
  friend colMatrix operator / (const double&, const colMatrix&);
  friend ostream & operator<< (ostream &os, const colMatrix&);

  friend void swap(colMatrix &first, colMatrix &second);

  //destructor
  ~colMatrix() {
    delete [] data;
    data = NULL;
  }

};

/*class to store an array of a fixed size equal to the number of variables being solved for. This is useful because a vector of these will be
contiguous in memory. */
class genArray {
  double data[NUMVARS];

 public:
  //constructor
  genArray() : data{0.0} {}
  genArray(const double& );
  genArray(const double &a, const double &b, const double &c, const double &d, const double &e ) : data{a, b, c, d, e} {}
  genArray(const double &a, const double &b, const double &c, const double &d, const double &e, const double &f, const double &g ) : data{a, b, c, d, e, f, g} {}

  //member functions
  void Zero();
  double Sum();

  //operator overloads
  genArray operator + (const genArray&)const;
  genArray operator - (const genArray&)const;
  genArray operator * (const genArray&)const;
  genArray operator / (const genArray&)const;

  genArray operator + (const vector<double>&)const;
  genArray operator - (const vector<double>&)const;

  genArray operator + (const double&)const;
  genArray operator - (const double&)const;
  genArray operator * (const double&)const;
  genArray operator / (const double&)const;

  double operator[] (const int &r)const{return data[r];}
  double& operator[] (const int &r){return data[r];}

  friend genArray operator + (const double&, const genArray&);
  friend genArray operator - (const double&, const genArray&);
  friend genArray operator * (const double&, const genArray&);
  friend genArray operator / (const double&, const genArray&);
  friend ostream & operator<< (ostream &os, const genArray&);

  void GlobalReduceMPI( const int&, const int& );

  //destructor
  ~genArray() {}

};

//class to store a square matrix
class squareMatrix {
  int size;
  double *data;

 public:
  //constructor
  squareMatrix( const int &a ): size(a) {data = new double[a*a];}
  squareMatrix() : size(0), data(NULL){}

  //copy constructor
  squareMatrix( const squareMatrix &cp);

  //copy assignment operator
  squareMatrix& operator= (squareMatrix other);
  //squareMatrix& operator= (const squareMatrix &other);

  //move constructor
  squareMatrix(squareMatrix &&other) : squareMatrix() {
    swap(*this,other);
    other.data=NULL;
  }

  //member functions
  double Data(const int &, const int &)const;
  void SetData(const int &, const int &, const double&);
  int Size()const{return size;}
  void SwapRows(const int &, const int &);
  void Inverse();
  int FindMaxInCol(const int &, const int &, const int &)const;
  void RowMultiply(const int &, const int &, const double &);
  void LinCombRow(const int &, const double &, const int &);
  void Zero();
  void Identity();
  colMatrix Multiply( const colMatrix & )const;

  //operator overloads
  squareMatrix operator + (const squareMatrix&)const;
  squareMatrix operator - (const squareMatrix&)const;
  squareMatrix operator * (const squareMatrix&)const;

  squareMatrix operator + (const double&)const;
  squareMatrix operator - (const double&)const;
  squareMatrix operator * (const double&)const;
  squareMatrix operator / (const double&)const;

  friend squareMatrix operator + (const double&, const squareMatrix&);
  friend squareMatrix operator - (const double&, const squareMatrix&);
  friend squareMatrix operator * (const double&, const squareMatrix&);
  friend squareMatrix operator / (const double&, const squareMatrix&);
  friend ostream & operator<< (ostream &os, const squareMatrix&);

  friend void swap(squareMatrix &first, squareMatrix &second);

  //destructor
  ~squareMatrix() {
    delete [] data;
    data = NULL;
  }

};

//Class to store the implicit flux jacobians for the entire mesh. Only values on populated diagonals are stored.
class matrixDiagonal {
  int size;
  squareMatrix *data;

 public:
  //constructor
  matrixDiagonal( const int &a ): size(a) {data = new squareMatrix[a];}
  matrixDiagonal() : size(0), data(NULL){}

  //copy constructor
  matrixDiagonal( const matrixDiagonal &cp);

  //copy assignment operator
  matrixDiagonal& operator= (matrixDiagonal other);

  //move constructor
  matrixDiagonal(matrixDiagonal &&other) : matrixDiagonal() {
    swap(*this,other);
    other.data=NULL;
  }

  //member functions
  squareMatrix Data(const int &)const;
  void SetData(const int &, const squareMatrix&);
  int Size()const{return size;}
  void Zero(const int &);
  void CleanResizeZero(const int &, const int &);
  void Inverse();

  //operator overloads
  friend ostream & operator<< (ostream &os, const matrixDiagonal&);

  friend void swap(matrixDiagonal &first, matrixDiagonal &second);

  //destructor
  ~matrixDiagonal() {
    delete [] data;
    data = NULL;
  }

};

class resid {
  double linf;
  int blk;
  int i;
  int j;
  int k;
  int eqn;

 public:
  //constructor
 resid() : linf(0.0), blk(0), i(0), j(0), k(0), eqn(0) {}
 resid( const double &a, const int &b, const int &c, const int &d, const int &e, const int &f) : linf(a), blk(b), i(c), j(d), k(e), eqn(f) {}

  //member functions
  double Linf()const{return linf;}
  int Block()const{return blk;}
  int ILoc()const{return i;}
  int JLoc()const{return j;}
  int KLoc()const{return k;}
  int Eqn()const{return eqn;}

  void UpdateMax(const double&, const int&, const int&, const int&, const int&, const int&);
  void GetAddressesMPI(MPI_Aint (&)[2])const;
  void GlobalReduceMPI( const int&, const MPI_Datatype&, const MPI_Op& );

  friend void MaxLinf( resid*, resid*, int*, MPI_Datatype*);

  void Zero(){
    linf = 0.0;
    blk = 0;
    i = 0;
    j = 0;
    k = 0;
    eqn = 0;
  }

  //destructor
  ~resid() {}

};


//function declarations



#endif
