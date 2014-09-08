#ifndef BOUNDARYCONDITIONSHEADERDEF             //only if the macro BOUNDARYCONDITIONSHEADERDEF is not defined execute these lines of code
#define BOUNDARYCONDITIONSHEADERDEF             //define the macro

/* This header contains the class boundaryConditions.

This class stores the information needed to specify the boundary conditions for one block. */

#include <vector>  //vector
#include <string>  //string
#include <iostream> //ostream

using std::ostream;
using std::vector;
using std::string;
using std::cout;
using std::endl;

class boundaryConditions {

  int numSurfI;                           //number of i-surfaces to define boundary on in block
  int numSurfJ;                           //number of j-surfaces to define boundary on in block
  int numSurfK;                           //number of k-surfaces to define boundary on in block
  vector<string> bcTypes;                 //vector of boundary condition names for each surface
  vector<int> iMin;                       //vector of min i coordinate defining surface
  vector<int> iMax;                       //vector of max i coordinate defining surface
  vector<int> jMin;                       //vector of min j coordinate defining surface
  vector<int> jMax;                       //vector of max j coordinate defining surface
  vector<int> kMin;                       //vector of min k coordinate defining surface
  vector<int> kMax;                       //vector of max k coordinate defining surface
  vector<int> tag;                        //vector of boundary condition tags

 public:
  //constructor
  boundaryConditions();
  boundaryConditions( int, int, int );

  //member functions
  void SetNumSurfI( const int &a){numSurfI = a;}
  int NumSurfI()const{return numSurfI;}
  void SetNumSurfJ( const int &a){numSurfJ = a;}
  int NumSurfJ()const{return numSurfJ;}
  void SetNumSurfK( const int &a){numSurfK = a;}
  int NumSurfK()const{return numSurfK;}

  void SetBCTypes( const string &str, const int &a){bcTypes[a] = str;}
  string GetBCTypes( const int &a)const{return bcTypes[a];}
  void SetIMin( const int &i, const int &a){iMin[a] = i;}
  int GetIMin( const int &a)const{return iMin[a];}
  void SetJMin( const int &j, const int &a){jMin[a] = j;}
  int GetJMin( const int &a)const{return jMin[a];}
  void SetKMin( const int &k, const int &a){kMin[a] = k;}
  int GetKMin( const int &a)const{return kMin[a];}
  void SetIMax( const int &i, const int &a){iMax[a] = i;}
  int GetIMax( const int &a)const{return iMax[a];}
  void SetJMax( const int &j, const int &a){jMax[a] = j;}
  int GetJMax( const int &a)const{return jMax[a];}
  void SetKMax( const int &k, const int &a){kMax[a] = k;}
  int GetKMax( const int &a)const{return kMax[a];}
  void SetTag( const int &t, const int &a){tag[a] = t;}
  int GetTag( const int &a)const{return tag[a];}

  void ResizeVecs( const int &a);

  friend ostream & operator<< (ostream &os, const boundaryConditions&);

  string GetBCName(const int, const int, const int, const string&)const;

  //destructor
  ~boundaryConditions() {}

};

//function declarations


#endif
