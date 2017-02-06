/*  This file is part of aither.
    Copyright (C) 2015-17  Michael Nucci (michael.nucci@gmail.com)

    Aither is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Aither is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>. */

#ifndef INPUTSTATESHEADERDEF
#define INPUTSTATESHEADERDEF

/* This header file contains the bcState and icState classes and their children

These states are used to store data from the input file and apply boundary
conditions and initial conditions.
*/

#include <string>
#include <vector>
#include <fstream>
#include <memory>
#include "vector3d.hpp"

using std::string;
using std::vector;
using std::ifstream;
using std::unique_ptr;

// this is an abstract base class
class inputState {
  int tag_;

 public:
  // constructor
  inputState() : tag_(-1) {}
  explicit inputState(const int &t) : tag_(t) {}

  // move constructor and assignment operator
  inputState(inputState&&) noexcept = default;
  inputState& operator=(inputState&&) noexcept = default;

  // copy constructor and assignment operator
  inputState(const inputState&) = default;
  inputState& operator=(const inputState&) = default;

  // member functions
  const int Tag() const {return tag_;}
  void SetTag(const int &t) {tag_ = t;}

  virtual void Print(ostream &os) const = 0;  // abstract base class
  virtual const vector3d<double> Velocity() const {return {0, 0, 0};}
  virtual const double Density() const {return 0;}
  virtual const double Pressure() const {return 0;}
  virtual const double TurbulenceIntensity() const {return 0;}
  virtual const double EddyViscosityRatio() const {return 0;}
  virtual const vector3d<double> Direction() const {return {0, 0, 0};}
  virtual const double StagnationPressure() const {return 0;}
  virtual const double StagnationTemperature() const {return 0;}
  virtual const double Temperature() const {return 0;}
  virtual const double HeatFlux() const {return 0;}
  virtual bool IsIsothermal() const {return false;}
  virtual bool IsAdiabatic() const {return false;}
  virtual bool IsConstantHeatFlux() const {return false;}

  // destructor
  virtual ~inputState() noexcept {}
};


class icState : public inputState {
  vector3d<double> velocity_;
  double density_;
  double pressure_;
  double turbIntensity_ = 0.01;       // default values
  double eddyViscRatio_ = 10.0;
  bool specifiedTurbulence_ = false;

 public:
  // constructor
  icState() : velocity_{0.0, 0.0, 0.0}, density_(0.0), pressure_(0.0) {}
  explicit icState(string &str, const string name = "icState");

  // move constructor and assignment operator
  icState(icState&&) noexcept = default;
  icState& operator=(icState&&) noexcept = default;

  // copy constructor and assignment operator
  icState(const icState&) = default;
  icState& operator=(const icState&) = default;

  // Member functions
  const vector3d<double> Velocity() const override {return velocity_;}
  const double Density() const override {return density_;}
  const double Pressure() const override {return pressure_;}
  const double TurbulenceIntensity() const override {return turbIntensity_;}
  const double EddyViscosityRatio() const override {return eddyViscRatio_;}
  const bool SpecifiedTurbulence() const {return specifiedTurbulence_;}
  void SetSpecifiedTurbulence() {specifiedTurbulence_ = true;}
  void Print(ostream &os) const override;

  // Destructor
  virtual ~icState() noexcept {}
};

// data for characteristic bc is same is for initial conditions
class characteristic : public icState {
 public:
  // constructor
  explicit characteristic(string &str) : icState(str, "characteristic") {}

  // move constructor and assignment operator
  characteristic(characteristic&&) noexcept = default;
  characteristic& operator=(characteristic&&) noexcept = default;

  // copy constructor and assignment operator
  characteristic(const characteristic&) = default;
  characteristic& operator=(const characteristic&) = default;

  // Member functions
  void Print(ostream &os) const override;

  // Destructor
  ~characteristic() noexcept {}
};


class stagnationInlet : public inputState {
  vector3d<double> direction_;
  double p0_;
  double t0_;
  double turbIntensity_ = 0.01;  // default values
  double eddyViscRatio_ = 10.0;
  bool specifiedTurbulence_ = false;

 public:
  // constructor
  explicit stagnationInlet(string &str);

  // move constructor and assignment operator
  stagnationInlet(stagnationInlet&&) noexcept = default;
  stagnationInlet& operator=(stagnationInlet&&) noexcept = default;

  // copy constructor and assignment operator
  stagnationInlet(const stagnationInlet&) = default;
  stagnationInlet& operator=(const stagnationInlet&) = default;

  // Member functions
  const vector3d<double> Direction() const override {return direction_;}
  const double StagnationPressure() const override {return p0_;}
  const double StagnationTemperature() const override {return t0_;}
  const double TurbulenceIntensity() const override {return turbIntensity_;}
  const double EddyViscosityRatio() const override {return eddyViscRatio_;}
  const bool SpecifiedTurbulence() const {return specifiedTurbulence_;}
  void SetSpecifiedTurbulence() {specifiedTurbulence_ = true;}
  void Print(ostream &os) const override;

  // Destructor
  ~stagnationInlet() noexcept {}
};


class pressureOutlet : public inputState {
  double pressure_;

 public:
  // constructor
  explicit pressureOutlet(string &str, const string name = "pressureOutlet");

  // move constructor and assignment operator
  pressureOutlet(pressureOutlet&&) noexcept = default;
  pressureOutlet& operator=(pressureOutlet&&) noexcept = default;

  // copy constructor and assignment operator
  pressureOutlet(const pressureOutlet&) = default;
  pressureOutlet& operator=(const pressureOutlet&) = default;

  // Member functions
  const double Pressure() const override {return pressure_;}
  void Print(ostream &os) const override;

  // Destructor
  virtual ~pressureOutlet() noexcept {}
};


// data for supersonic inflow bc is same is for initial conditions
class supersonicInflow : public icState {
 public:
  // constructor
  explicit supersonicInflow(string &str) : icState(str, "supersonicInflow") {}

  // move constructor and assignment operator
  supersonicInflow(supersonicInflow&&) noexcept = default;
  supersonicInflow& operator=(supersonicInflow&&) noexcept = default;

  // copy constructor and assignment operator
  supersonicInflow(const supersonicInflow&) = default;
  supersonicInflow& operator=(const supersonicInflow&) = default;

  // Member functions
  void Print(ostream &os) const override;

  // Destructor
  ~supersonicInflow() noexcept {}
};


// data for subsonicOutflow bc is same is for pressureOutlet
class subsonicOutflow : public pressureOutlet {
 public:
  // constructor
  explicit subsonicOutflow(string &str) :
      pressureOutlet(str, "subsonicOutflow") {}

  // move constructor and assignment operator
  subsonicOutflow(subsonicOutflow&&) noexcept = default;
  subsonicOutflow& operator=(subsonicOutflow&&) noexcept = default;

  // copy constructor and assignment operator
  subsonicOutflow(const subsonicOutflow&) = default;
  subsonicOutflow& operator=(const subsonicOutflow&) = default;

  // Member functions
  void Print(ostream &os) const override;

  // Destructor
  ~subsonicOutflow() noexcept {}
};


class subsonicInflow : public inputState {
  vector3d<double> velocity_;
  double density_;
  double turbIntensity_ = 0.01;  // default values
  double eddyViscRatio_ = 10.0;
  bool specifiedTurbulence_ = false;

 public:
  // constructor
  explicit subsonicInflow(string &str);

  // move constructor and assignment operator
  subsonicInflow(subsonicInflow&&) noexcept = default;
  subsonicInflow& operator=(subsonicInflow&&) noexcept = default;

  // copy constructor and assignment operator
  subsonicInflow(const subsonicInflow&) = default;
  subsonicInflow& operator=(const subsonicInflow&) = default;

  // Member functions
  const vector3d<double> Velocity() const override {return velocity_;}
  const double Density() const override {return density_;}
  const double TurbulenceIntensity() const override {return turbIntensity_;}
  const double EddyViscosityRatio() const override {return eddyViscRatio_;}
  const bool SpecifiedTurbulence() const {return specifiedTurbulence_;}
  void SetSpecifiedTurbulence() {specifiedTurbulence_ = true;}
  void Print(ostream &os) const override;

  // Destructor
  ~subsonicInflow() noexcept {}
};


class viscousWall : public inputState {
  // default conditions for stationary adiabatic wall
  vector3d<double> velocity_ = {0.0, 0.0, 0.0};
  double temperature_ = 0.0;
  double heatFlux_ = 0.0;
  bool specifiedTemperature_ = false;
  bool specifiedHeatFlux_ = false;

 public:
  // constructor
  explicit viscousWall(string &str);

  // move constructor and assignment operator
  viscousWall(viscousWall&&) noexcept = default;
  viscousWall& operator=(viscousWall&&) noexcept = default;

  // copy constructor and assignment operator
  viscousWall(const viscousWall&) = default;
  viscousWall& operator=(const viscousWall&) = default;

  // Member functions
  const vector3d<double> Velocity() const override {return velocity_;}
  const double Temperature() const override {return temperature_;}
  const double HeatFlux() const override {return heatFlux_;}
  bool IsIsothermal() const override {
    return specifiedTemperature_ ? true : false;
  }
  bool IsAdiabatic() const override {
    return (!specifiedTemperature_ && heatFlux_ == 0.0) ? true : false;
  }
  bool IsConstantHeatFlux() const override {
    return (specifiedHeatFlux_ && heatFlux_ != 0.0) ? true : false;
  }
  void Print(ostream &os) const override;

  // Destructor
  virtual ~viscousWall() noexcept {}
};



// function declarations
ostream &operator<<(ostream &, const inputState &);
ostream &operator<<(ostream &, const icState &);
ostream &operator<<(ostream &, const characteristic &);
ostream &operator<<(ostream &, const stagnationInlet &);
ostream &operator<<(ostream &, const pressureOutlet &);
ostream &operator<<(ostream &, const supersonicInflow &);
ostream &operator<<(ostream &, const subsonicOutflow &);
ostream &operator<<(ostream &, const subsonicInflow &);
ostream &operator<<(ostream &, const viscousWall &);


vector<string> Tokenize(string, const string &, const unsigned int = 0);
string Trim(const string &, const string & = " \t\r\n");
vector3d<double> ReadVector(const string &);
vector<icState> ReadICList(ifstream &, string &);
vector<string> ReadStringList(ifstream &, string &);
vector<unique_ptr<inputState>> ReadBCList(ifstream &, string &);
string RemoveTrailing(const string &, const string &);
auto FindBCPosition(const string &, const vector<string> &, string &);
void AddBCToList(const string &, vector<unique_ptr<inputState>> &, string &);
void CheckICTags(const vector<icState> &, const int &);

#endif
