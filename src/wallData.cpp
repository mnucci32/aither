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

#include "wallData.hpp"
#include "vector3d.hpp"
#include "input.hpp"
#include "eos.hpp"
#include "primitive.hpp"

// wallVars functions
void wallVars::Pack(char *(&sendBuffer), const int &sendBufSize, int &position,
                    const MPI_Datatype &MPI_vec3d) const {
  // sendBuffer -- buffer to pack data into
  // sendBufSize -- size of buffer
  // position -- location within buffer
  // MPI_vec3d -- datatype for vector3d

  // pack wall shear stress
  MPI_Pack(&shearStress_, 1, MPI_vec3d, sendBuffer, sendBufSize, &position,
           MPI_COMM_WORLD);
  // pack wall scalars
  MPI_Pack(&heatFlux_, 1, MPI_DOUBLE, sendBuffer, sendBufSize, &position,
           MPI_COMM_WORLD);
  MPI_Pack(&yplus_, 1, MPI_DOUBLE, sendBuffer, sendBufSize, &position,
           MPI_COMM_WORLD);
  MPI_Pack(&temperature_, 1, MPI_DOUBLE, sendBuffer, sendBufSize, &position,
           MPI_COMM_WORLD);
  MPI_Pack(&turbEddyVisc_, 1, MPI_DOUBLE, sendBuffer, sendBufSize, &position,
           MPI_COMM_WORLD);
  MPI_Pack(&viscosity_, 1, MPI_DOUBLE, sendBuffer, sendBufSize, &position,
           MPI_COMM_WORLD);
  MPI_Pack(&density_, 1, MPI_DOUBLE, sendBuffer, sendBufSize, &position,
           MPI_COMM_WORLD);
  MPI_Pack(&frictionVelocity_, 1, MPI_DOUBLE, sendBuffer, sendBufSize,
           &position, MPI_COMM_WORLD);
  MPI_Pack(&tke_, 1, MPI_DOUBLE, sendBuffer, sendBufSize, &position,
           MPI_COMM_WORLD);
  MPI_Pack(&sdr_, 1, MPI_DOUBLE, sendBuffer, sendBufSize, &position,
           MPI_COMM_WORLD);
  MPI_Pack(&(*std::begin(mf_)), mf_.size(), MPI_DOUBLE, sendBuffer, sendBufSize,
           &position, MPI_COMM_WORLD);
}

void wallVars::PackSize(int &sendBufSize, const MPI_Datatype &MPI_vec3d) const {
  auto tempSize = 0;
  // add size for shear stress
  MPI_Pack_size(1, MPI_vec3d, MPI_COMM_WORLD, &tempSize);
  sendBufSize += tempSize;
  // 9 because 9 scalar varialbes
  MPI_Pack_size(9, MPI_DOUBLE, MPI_COMM_WORLD, &tempSize);
  sendBufSize += tempSize;
  // add size for mf
  MPI_Pack_size(mf_.size(), MPI_DOUBLE, MPI_COMM_WORLD, &tempSize);
  sendBufSize += tempSize;
}

void wallVars::Unpack(char *(&recvBuffer), const int &recvBufSize,
                      int &position, const MPI_Datatype &MPI_vec3d,
                      const int &numSpecies) {
  // recvBuffer -- buffer to unpack data from
  // recvBufSize -- size of buffer
  // position -- location within buffer

  // unpack shear stress
  MPI_Unpack(recvBuffer, recvBufSize, &position, &shearStress_, 1, MPI_vec3d,
             MPI_COMM_WORLD);
  // unpack wall scalars
  MPI_Unpack(recvBuffer, recvBufSize, &position, &heatFlux_, 1, MPI_DOUBLE,
             MPI_COMM_WORLD);
  MPI_Unpack(recvBuffer, recvBufSize, &position, &yplus_, 1, MPI_DOUBLE,
             MPI_COMM_WORLD);
  MPI_Unpack(recvBuffer, recvBufSize, &position, &temperature_, 1, MPI_DOUBLE,
             MPI_COMM_WORLD);
  MPI_Unpack(recvBuffer, recvBufSize, &position, &turbEddyVisc_, 1, MPI_DOUBLE,
             MPI_COMM_WORLD);
  MPI_Unpack(recvBuffer, recvBufSize, &position, &viscosity_, 1, MPI_DOUBLE,
             MPI_COMM_WORLD);
  MPI_Unpack(recvBuffer, recvBufSize, &position, &density_, 1, MPI_DOUBLE,
             MPI_COMM_WORLD);
  MPI_Unpack(recvBuffer, recvBufSize, &position, &frictionVelocity_, 1,
             MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Unpack(recvBuffer, recvBufSize, &position, &tke_, 1, MPI_DOUBLE,
             MPI_COMM_WORLD);
  MPI_Unpack(recvBuffer, recvBufSize, &position, &sdr_, 1, MPI_DOUBLE,
             MPI_COMM_WORLD);
  // unpack mass fractions
  mf_.clear();
  mf_.resize(numSpecies);
  MPI_Unpack(recvBuffer, recvBufSize, &position, &(*std::begin(mf_)),
             numSpecies, MPI_DOUBLE, MPI_COMM_WORLD);
}


// member functions
vector3d<double> wallData::WallShearStress(const int &ii, const int &jj,
                                           const int &kk) const {
  return (*this)(ii, jj, kk).shearStress_;
}

double wallData::WallHeatFlux(const int &ii, const int &jj,
                              const int &kk) const {
  return (*this)(ii, jj, kk).heatFlux_;
}

double wallData::Yplus(const int &ii, const int &jj, const int &kk) const {
  return (*this)(ii, jj, kk).yplus_;
}

double wallData::WallTke(const int &ii, const int &jj, const int &kk) const {
  return (*this)(ii, jj, kk).tke_;
}

double wallData::WallSdr(const int &ii, const int &jj, const int &kk) const {
  return (*this)(ii, jj, kk).sdr_;
}

double wallData::WallTemperature(const int &ii, const int &jj,
                                 const int &kk) const {
  return (*this)(ii, jj, kk).temperature_;
}

double wallData::WallEddyViscosity(const int &ii, const int &jj,
                                   const int &kk) const {
  return (*this)(ii, jj, kk).turbEddyVisc_;
}

double wallData::WallPressure(const int &ii, const int &jj, const int &kk,
                              const unique_ptr<eos> &eqnState) const {
  return eqnState->PressureRT(this->WallDensityVec(ii, jj, kk),
                              this->WallTemperature(ii, jj, kk));
}

double wallData::WallViscosity(const int &ii, const int &jj,
                               const int &kk) const {
  return (*this)(ii, jj, kk).viscosity_;
}

double wallData::WallDensity(const int &ii, const int &jj,
                             const int &kk) const {
  return (*this)(ii, jj, kk).density_;
}

vector<double> wallData::WallMassFractions(const int &ii, const int &jj,
                                           const int &kk) const {
  return (*this)(ii, jj, kk).mf_;
}

vector<double> wallData::WallDensityVec(const int &ii, const int &jj,
                                        const int &kk) const {
  const auto rho = this->WallDensity(ii, jj, kk);
  auto rhoVec = this->WallMassFractions(ii, jj, kk);
  std::for_each(rhoVec.begin(), rhoVec.end(),
                [&rho](double &val) { val *= rho; });
  return rhoVec;
}

double wallData::WallFrictionVelocity(const int &ii, const int &jj,
                                      const int &kk) const {
  return (*this)(ii, jj, kk).frictionVelocity_;
}

void wallData::PackWallData(char *(&sendBuffer), const int &sendBufSize,
                            int &position,
                            const MPI_Datatype &MPI_vec3d) const {
  // sendBuffer -- buffer to pack data into
  // sendBufSize -- size of buffer
  // position -- location within buffer

  // pack force counters
  MPI_Pack(&inviscidForce_, 1, MPI_DOUBLE, sendBuffer, sendBufSize, &position,
           MPI_COMM_WORLD);
  MPI_Pack(&viscousForce_, 1, MPI_DOUBLE, sendBuffer, sendBufSize, &position,
           MPI_COMM_WORLD);

  // pack number of species counter
  MPI_Pack(&numSpecies_, 1, MPI_INT, sendBuffer, sendBufSize, &position,
           MPI_COMM_WORLD);

  // pointer to bc data - remote processor can get data from input class

  // pack boundarySurface
  surf_.PackBoundarySurface(sendBuffer, sendBufSize, position);

  // pack wall variables
  for (auto &wv : data_) {
    wv.Pack(sendBuffer, sendBufSize, position, MPI_vec3d);
  }
}

void wallData::PackSize(int &sendBufSize, const MPI_Datatype &MPI_vec3d) const {
  auto tempSize = 0;
  // add sizes for force data
  MPI_Pack_size(2, MPI_DOUBLE, MPI_COMM_WORLD, &tempSize);
  sendBufSize += tempSize;
  // add sizes for number of species
  MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &tempSize);
  sendBufSize += tempSize;
  // 8 because iMin, iMax, jMin, jMax, kMin, kMax, tags, string sizes
  MPI_Pack_size(8, MPI_INT, MPI_COMM_WORLD, &tempSize);
  sendBufSize += tempSize;
  // add size for bc types (+1 for c_str end character)
  MPI_Pack_size(surf_.BCType().size() + 1, MPI_CHAR, MPI_COMM_WORLD, &tempSize);
  sendBufSize += tempSize;

  // add array of wallData
  for (auto &wv : data_) {
    wv.PackSize(sendBufSize, MPI_vec3d);
  }
}

void wallData::UnpackWallData(char *(&recvBuffer), const int &recvBufSize,
                              int &position, const MPI_Datatype &MPI_vec3d,
                              const input &inp) {
  // recvBuffer -- buffer to unpack data from
  // recvBufSize -- size of buffer
  // position -- location within buffer

  // unpack forces
  MPI_Unpack(recvBuffer, recvBufSize, &position, &inviscidForce_, 1, MPI_DOUBLE,
             MPI_COMM_WORLD);
  MPI_Unpack(recvBuffer, recvBufSize, &position, &viscousForce_, 1, MPI_DOUBLE,
             MPI_COMM_WORLD);

  // unpack number of species
  MPI_Unpack(recvBuffer, recvBufSize, &position, &numSpecies_, 1, MPI_INT,
             MPI_COMM_WORLD);

  // unpack BC surface
  surf_.UnpackBoundarySurface(recvBuffer, recvBufSize, position);
  
  // get bc data from tag
  bcData_ = inp.BCData(surf_.Tag());

  // unpack wall variables
  data_.ClearResize(surf_.NumI(), surf_.NumJ(), surf_.NumK(), 0,
                    wallVars(numSpecies_));
  for (auto &wv : data_) {
    wv.Unpack(recvBuffer, recvBufSize, position, MPI_vec3d, numSpecies_);
  }
}

// Split wallData at given direction and index
// The calling instance is modified to be the lower surface, and the upper
// surface is returned. If the calling instance is not split, the split
// boolean returns false and the returned wallData is garbage. The force
// variables are not split.
wallData wallData::Split(const string &dir, const int &ind, bool &split,
                         bool &low) {
  wallData upper = *this;
  // split if necessary
  auto upperSurf = surf_.Split(dir, ind, split, low);
  if (split) {  // surface split; upper and lower valid
    upper.surf_ = upperSurf;
    upper.data_ = data_.Slice(dir, {ind, data_.End(dir)});
    data_ = data_.Slice(dir, {data_.Start(dir), ind});
  } else if (!low) {  // not split; upper is valid
    *this = wallData();  // invalidate lower; upper already equals *this
  } else {               // if not split and lower is valid, invalidate upper
    upper = wallData();
  }
  return upper;
}

void wallData::Join(const wallData &upper, const string &dir, bool &joined) {
  // modify upper surface to see if it is possible to join
  auto upSurfMod = upper.surf_;
  upSurfMod.IncrementDirection(dir, surf_.Max(dir));

  // join if possible
  surf_.Join(upSurfMod, dir, joined);
  if (joined) {
    inviscidForce_ += upper.inviscidForce_;
    viscousForce_ += upper.viscousForce_;

    multiArray3d<wallVars> newVars(surf_.NumI(), surf_.NumJ(), surf_.NumK(), 0,
                                   1, wallVars(numSpecies_));
    newVars.Insert(dir, {data_.Start(dir), data_.PhysEnd(dir)},
                   data_.Slice(dir, {data_.Start(dir), data_.PhysEnd(dir)}));
    newVars.Insert(dir, {data_.PhysEnd(dir), newVars.End(dir)},
                   upper.data_.Slice(dir, {upper.data_.PhysStart(dir),
                                           upper.data_.End(dir)}));
    data_ = newVars;
  }
}

  void wallData::WallState(const int &ii, const int &jj, const int &kk,
                     const unique_ptr<eos> &eqnState, primitive &wState) const {
    auto rhoVec = this->WallDensityVec(ii, jj, kk);
    for (auto ii = 0; ii < wState.NumSpecies(); ++ii) {
      wState[ii] = rhoVec[ii];
    }
    wState[wState.MomentumXIndex()] = this->WallVelocity().X();
    wState[wState.MomentumYIndex()] = this->WallVelocity().Y();
    wState[wState.MomentumZIndex()] = this->WallVelocity().Z();
    wState[wState.EnergyIndex()] = this->WallPressure(ii, jj, kk, eqnState);
    for (auto ii = 0; ii < wState.NumTurbulence(); ++ii) {
      wState[wState.TurbulenceIndex() + ii] =
          (ii == 0) ? this->WallTke(ii, jj, kk) : this->WallSdr(ii, jj, kk);
    }
  }

  void wallData::Print(ostream &os) const {
    os << "Inviscid Force: " << inviscidForce_ << endl;
    os << "Viscous Force: " << viscousForce_ << endl;
    os << "Number of Species: " << numSpecies_ << endl;
    os << "BC Data: ";
    bcData_->Print(os);
    os << endl;
    os << "BC Surface: " << surf_ << endl;
    os << "Wall Data:" << endl;
    os << data_ << endl;
}

ostream &operator<<(ostream &os, const wallData &wd) {
  wd.Print(os);
  return os;
}

ostream &operator<<(ostream &os, const wallVars &wv) {
  os << wv.shearStress_ << "; " << wv.heatFlux_ << "; " << wv.yplus_ << "; "
     << wv.temperature_ << "; " << wv.turbEddyVisc_ << "; " << wv.viscosity_
     << "; " << wv.density_ << "; " << wv.frictionVelocity_ << "; " << wv.tke_
     << "; " << wv.sdr_;
  return os;
}