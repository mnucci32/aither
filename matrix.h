#ifndef MATRIXHEADERDEF             //only if the macro MATRIXHEADERDEF is not defined execute these lines of code
#define MATRIXHEADERDEF             //define the macro

#include <vector>  //vector
#include <string>  //string
#include <math.h>       //sqrt
#include <iostream>

using std::vector;
using std::string;
using std::ostream;

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

//function declarations


#endif
