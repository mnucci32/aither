/*  This file is part of aither.
    Copyright (C) 2015-18  Michael Nucci (mnucci@pm.me)

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
#include <map>
#include <set>
#include <fstream>
#include <memory>
#include "vector3d.hpp"

using std::string;
using std::vector;
using std::map;
using std::set;
using std::ifstream;
using std::shared_ptr;

#define DEFAULT_TURB_INTENSITY 0.01
#define DEFAULT_EDDY_VISC_RATIO 0.01

// this is an abstract base class
class inputState {
  int tag_;
  bool nondimensional_ = false;

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
  const int StartTag() const {return this->Tag();}

  virtual const int EndTag() const {return this->Tag();}
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
  virtual bool IsTranslation() const {return false;}
  virtual bool IsRotation() const {return false;}
  virtual const double Rotation() const {return 0;}
  virtual const vector3d<double> Translation() const {return {0, 0, 0};}
  virtual const vector3d<double> Axis() const {return {0, 0, 0};}
  virtual const vector3d<double> Point() const {return {0, 0, 0};}
  virtual const double VonKarmen() const {return 0.0;}
  virtual const double WallConstant() const {return 0.0;}
  virtual bool IsWallLaw() const {return false;}
  virtual void Nondimensionalize(const double &rRef, const double &tRef,
                                 const double &lRef, const double &aRef) = 0;
  virtual double MassFraction(const string &species) const { return 0.0; };
  virtual map<string, double> MassFractions() const {
    return map<string, double>();
  }
  virtual int NumberSpecies() const { return 0; }
  virtual string File() const { return "undefined"; }
  virtual bool IsFromFile() const { return false; }
  virtual bool IsNonreflecting() const { return false; }
  virtual const double LengthScale() const { return 0.0; }
  bool IsNondimensional() const { return nondimensional_; }
  void SetNondimensional(const bool &nd) {nondimensional_ = nd;}
  virtual void Read(string &str) = 0;
  virtual string Name() const { return "inputState"; }

  // destructor
  virtual ~inputState() noexcept {}
};


class icState : public inputState {
  vector3d<double> velocity_ = {0.0, 0.0, 0.0};
  double density_ = 0.0;
  double pressure_ = 0.0;
  double turbIntensity_ = DEFAULT_TURB_INTENSITY;
  double eddyViscRatio_ = DEFAULT_EDDY_VISC_RATIO;
  map<string, double> massFractions_ = {{"air", 1.0}};
  string file_ = "undefined";
  bool specifiedTurbulence_ = false;
  bool specifiedMassFractions_ = false;
  bool specifiedFile_ = false;
  string name_ = "icState";

 public:
  // constructor
  icState() : inputState() {}

  // move constructor and assignment operator
  icState(icState&&) noexcept = default;
  icState& operator=(icState&&) noexcept = default;

  // copy constructor and assignment operator
  icState(const icState&) = default;
  icState& operator=(const icState&) = default;

  // Member functions
  void Read(string &str) override;
  string Name() const override { return name_; }
  void UpdateName(const string &n) { name_ = n; }
  const vector3d<double> Velocity() const override { return velocity_; }
  const double Density() const override {return density_;}
  const double Pressure() const override {return pressure_;}
  const double TurbulenceIntensity() const override {return turbIntensity_;}
  const double EddyViscosityRatio() const override {return eddyViscRatio_;}
  const bool SpecifiedTurbulence() const {return specifiedTurbulence_;}
  void SetSpecifiedTurbulence() {specifiedTurbulence_ = true;}
  const bool SpecifiedMassFractions() const {return specifiedMassFractions_;}
  void SetSpecifiedMassFractions() {specifiedMassFractions_ = true;}
  int NumberSpecies() const override {return massFractions_.size();}
  void Print(ostream &os) const override;
  void Nondimensionalize(const double &rRef, const double &tRef,
                         const double &lRef, const double &aRef) override;
  virtual void NondimensionalizeExtra(const double &rRef, const double &tRef,
                         const double &lRef, const double &aRef) {}
  double MassFraction(const string &species) const override {
    return massFractions_.find(species)->second;
  }
  map<string, double> MassFractions() const override { return massFractions_; }
  const bool SpecifiedFile() const { return specifiedFile_; }
  void SetSpecifiedFile() { specifiedFile_ = true; }
  string File() const override { return file_; }
  bool IsFromFile() const override { return specifiedFile_; }
  virtual void AssignExtraData(const string &s1, const string &s2) {}
  virtual void ExtraDataChecks() const {}
  virtual bool MatchesExtraData(const string &s1) const { return false; }

  // Destructor
  virtual ~icState() noexcept {}
};

// data for characteristic bc is same is for initial conditions
class characteristic : public icState {
 public:
  // constructor
  characteristic() : icState() { this->UpdateName("characteristic"); }

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

// data for inlet bc is same is for initial conditions
class inlet : public icState {
  bool nonreflecting_ = false;
  bool specifiedReflecting_ = false;
  int nonreflectingCount_ = 0;
  double lengthScale_ = 0.0;
  int lengthCount_ = 0;
  set<string> extraData_ = {"nonreflecting", "lengthScale"};

 public:
  // constructor
  inlet() : icState() { this->UpdateName("inlet"); }

  // move constructor and assignment operator
  inlet(inlet&&) noexcept = default;
  inlet& operator=(inlet&&) noexcept = default;

  // copy constructor and assignment operator
  inlet(const inlet&) = default;
  inlet& operator=(const inlet&) = default;

  // Member functions
  void Print(ostream &os) const override;
  bool IsNonreflecting() const override { return nonreflecting_; }
  bool SpecifiedReflecting() const { return specifiedReflecting_; }
  const double LengthScale() const override { return lengthScale_; }
  void NondimensionalizeExtra(const double &rRef, const double &tRef,
                              const double &lRef, const double &aRef) override;
  void AssignExtraData(const string &s1, const string &s2) override;
  void ExtraDataChecks() const override;
  bool MatchesExtraData(const string &s1) const override {
    return extraData_.find(s1) != extraData_.end();
  }

  // Destructor
  ~inlet() noexcept {}
};


class stagnationInlet : public inputState {
  vector3d<double> direction_ = {0.0, 0.0, 0.0};
  double p0_ = 0.0;
  double t0_ = 0.0;
  double turbIntensity_ = DEFAULT_TURB_INTENSITY;
  double eddyViscRatio_ = DEFAULT_EDDY_VISC_RATIO;
  map<string, double> massFractions_ = {{"air", 1.0}};
  bool specifiedTurbulence_ = false;
  bool specifiedMassFractions_ = false;
  string name_ = "stagnationInlet";

 public:
  // constructor
  stagnationInlet() : inputState() {}

  // move constructor and assignment operator
  stagnationInlet(stagnationInlet&&) noexcept = default;
  stagnationInlet& operator=(stagnationInlet&&) noexcept = default;

  // copy constructor and assignment operator
  stagnationInlet(const stagnationInlet&) = default;
  stagnationInlet& operator=(const stagnationInlet&) = default;

  // Member functions
  void Read(string &str) override;
  string Name() const override { return name_; }
  const vector3d<double> Direction() const override { return direction_; }
  const double StagnationPressure() const override {return p0_;}
  const double StagnationTemperature() const override {return t0_;}
  const double TurbulenceIntensity() const override {return turbIntensity_;}
  const double EddyViscosityRatio() const override {return eddyViscRatio_;}
  const bool SpecifiedTurbulence() const {return specifiedTurbulence_;}
  void SetSpecifiedTurbulence() {specifiedTurbulence_ = true;}
  const bool SpecifiedMassFractions() const {return specifiedMassFractions_;}
  void SetSpecifiedMassFractions() {specifiedMassFractions_ = true;}
  int NumberSpecies() const override {return massFractions_.size();}
  void Print(ostream &os) const override;
  void Nondimensionalize(const double &rRef, const double &tRef,
                         const double &lRef, const double &aRef) override;
  double MassFraction(const string &species) const override {
    return massFractions_.find(species)->second;
  }
  map<string, double> MassFractions() const override { return massFractions_; }

  // Destructor
  ~stagnationInlet() noexcept {}
};


class pressureOutlet : public inputState {
  double pressure_ = 0.0;
  bool nonreflecting_ = false;
  bool specifiedReflecting_ = false;
  double lengthScale_ = 0.0;
  string name_ = "pressureOutlet";

 public:
  // constructor
  pressureOutlet() : inputState() {}

  // move constructor and assignment operator
  pressureOutlet(pressureOutlet&&) noexcept = default;
  pressureOutlet& operator=(pressureOutlet&&) noexcept = default;

  // copy constructor and assignment operator
  pressureOutlet(const pressureOutlet&) = default;
  pressureOutlet& operator=(const pressureOutlet&) = default;

  // Member functions
  void Read(string &str) override;
  string Name() const override { return name_; }
  const double Pressure() const override {return pressure_;}
  bool IsNonreflecting() const override { return nonreflecting_; }
  bool SpecifiedReflecting() const { return specifiedReflecting_; }
  const double LengthScale() const override { return lengthScale_; }
  void Print(ostream &os) const override;
  void Nondimensionalize(const double &rRef, const double &tRef,
                         const double &lRef, const double &aRef) override;

  // Destructor
  virtual ~pressureOutlet() noexcept {}
};


// data for supersonic inflow bc is same is for initial conditions
class supersonicInflow : public icState {
 public:
  // constructor
  supersonicInflow() : icState() { this->UpdateName("supersonicInflow"); }

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


class viscousWall : public inputState {
  // default conditions for stationary adiabatic wall
  vector3d<double> velocity_ = {0.0, 0.0, 0.0};
  double temperature_ = 0.0;
  double heatFlux_ = 0.0;
  double vonKarmen_ = 0.41;
  double wallConstant_ = 5.5;
  string wallTreatment_ = "lowRe";
  bool specifiedTemperature_ = false;
  bool specifiedHeatFlux_ = false;
  string name_ = "viscousWall";

 public:
  // constructor
  viscousWall() : inputState() {}

  // move constructor and assignment operator
  viscousWall(viscousWall&&) noexcept = default;
  viscousWall& operator=(viscousWall&&) noexcept = default;

  // copy constructor and assignment operator
  viscousWall(const viscousWall&) = default;
  viscousWall& operator=(const viscousWall&) = default;

  // Member functions
  void Read(string &str) override;
  string Name() const override { return name_; }
  const vector3d<double> Velocity() const override { return velocity_; }
  const double Temperature() const override {return temperature_;}
  const double HeatFlux() const override {return heatFlux_;}
  const double VonKarmen() const override {return vonKarmen_;}
  const double WallConstant() const override {return wallConstant_;}
  bool IsWallLaw() const override {
    return (wallTreatment_ == "wallLaw") ? true : false;
  }
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
  void Nondimensionalize(const double &rRef, const double &tRef,
                         const double &lRef, const double &aRef) override;

  // Destructor
  virtual ~viscousWall() noexcept {}
};


class periodic : public inputState {
  vector3d<double> translation_ = {0.0, 0.0, 0.0};
  vector3d<double> axis_ = {0.0, 0.0, 0.0};
  vector3d<double> point_ = {0.0, 0.0, 0.0};
  double rotation_ = 0.0;
  int endTag_ = -1;
  string name_ = "periodic";

 public:
  // constructor
  periodic() : inputState() {}

  // move constructor and assignment operator
  periodic(periodic&&) noexcept = default;
  periodic& operator=(periodic&&) noexcept = default;

  // copy constructor and assignment operator
  periodic(const periodic&) = default;
  periodic& operator=(const periodic&) = default;

  // Member functions
  void Read(string &str) override;
  string Name() const override { return name_; }
  bool IsTranslation() const override {
    return translation_ != vector3d<double>(0.0, 0.0, 0.0);
  }
  bool IsRotation() const override {
    return axis_ != vector3d<double>(0.0, 0.0, 0.0);
  }
  const vector3d<double> Translation() const override {return translation_;}
  const vector3d<double> Axis() const override {return axis_;}
  const vector3d<double> Point() const override {return point_;}
  const double Rotation() const override {return rotation_;}
  const int EndTag() const override {return endTag_;}
  void Print(ostream &os) const override;
  void Nondimensionalize(const double &rRef, const double &tRef,
                         const double &lRef, const double &aRef) override;

  // Destructor
  virtual ~periodic() noexcept {}
};


// function declarations
ostream &operator<<(ostream &, const inputState &);
ostream &operator<<(ostream &, const icState &);
ostream &operator<<(ostream &, const characteristic &);
ostream &operator<<(ostream &, const inlet &);
ostream &operator<<(ostream &, const stagnationInlet &);
ostream &operator<<(ostream &, const pressureOutlet &);
ostream &operator<<(ostream &, const supersonicInflow &);
ostream &operator<<(ostream &, const viscousWall &);
ostream &operator<<(ostream &, const periodic &);


vector<string> Tokenize(string, const string &, const unsigned int = 0);
string Trim(const string &, const string & = " \t\r\n");
vector3d<double> ReadVector(const string &);
vector<double> ReadVectorXd(const string &);
map<string, double> ReadMassFractions(const string &);
vector<icState> ReadICList(ifstream &, string &);
vector<string> ReadStringList(ifstream &, string &);
vector<shared_ptr<inputState>> ReadBCList(ifstream &, string &);
string RemoveTrailing(const string &, const string &);
auto FindBCPosition(const string &, const vector<string> &, string &);
void AddBCToList(const string &, vector<shared_ptr<inputState>> &, string &);
void CheckICTags(const vector<icState> &, const int &);

#endif
