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
  const int size;
  double *data;

 public:
  //constructor
  squareMatrix( const int &a ): size(a) {data = new double[a];}

  //member functions
  double Data(const int &, const int &)const;
  void SetData(const int &, const int &, const double&);
  int Size()const{return size;}

  //operator overloads
  squareMatrix operator + (const squareMatrix&)const;
  squareMatrix operator - (const squareMatrix&)const;
  squareMatrix operator * (const squareMatrix&)const;

  friend squareMatrix operator + (const double&, const squareMatrix&);
  friend squareMatrix operator - (const double&, const squareMatrix&);
  friend squareMatrix operator * (const double&, const squareMatrix&);
  friend squareMatrix operator / (const double&, const squareMatrix&);
  friend ostream & operator<< (ostream &os, const squareMatrix&);

  //destructor
  ~squareMatrix() {delete [] data;}

};

//function declarations


#endif
