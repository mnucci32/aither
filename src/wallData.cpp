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
#include "primative.hpp"

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
  return eqnState->PressureRT(this->WallDensity(ii, jj, kk),
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

double wallData::WallFrictionVelocity(const int &ii, const int &jj,
                                      const int &kk) const {
  return (*this)(ii, jj, kk).frictionVelocity_;
}

void wallData::PackWallData(char *(&sendBuffer), const int &sendBufSize,
                            int &position,
                            const MPI_Datatype &MPI_wallData) const {
  // sendBuffer -- buffer to pack data into
  // sendBufSize -- size of buffer
  // position -- location within buffer

  // pack force counters
  MPI_Pack(&inviscidForce_, 1, MPI_DOUBLE, sendBuffer, sendBufSize, &position,
           MPI_COMM_WORLD);
  MPI_Pack(&viscousForce_, 1, MPI_DOUBLE, sendBuffer, sendBufSize, &position,
           MPI_COMM_WORLD);

  // pointer to bc data - remote processor can get data from input class

  // pack boundarySurface
  surf_.PackBoundarySurface(sendBuffer, sendBufSize, position);

  // pack wall variables
  MPI_Pack(&(*std::begin(data_)), data_.Size(), MPI_wallData, sendBuffer,
           sendBufSize, &position, MPI_COMM_WORLD);
}

void wallData::PackSize(int &sendBufSize,
                        const MPI_Datatype &MPI_wallData) const {
  auto tempSize = 0;
  // add sizes for force data
  MPI_Pack_size(2, MPI_DOUBLE, MPI_COMM_WORLD, &tempSize);
  sendBufSize += tempSize;
  // 8 because iMin, iMax, jMin, jMax, kMin, kMax, tags, string sizes
  MPI_Pack_size(8, MPI_INT, MPI_COMM_WORLD, &tempSize);
  sendBufSize += tempSize;
  // add size for bc types (+1 for c_str end character)
  MPI_Pack_size(surf_.BCType().size() + 1, MPI_CHAR, MPI_COMM_WORLD, &tempSize);
  sendBufSize += tempSize;

  // add array of wallData
  MPI_Pack_size(data_.Size(), MPI_wallData, MPI_COMM_WORLD, &tempSize);
  sendBufSize += tempSize;
}

void wallData::UnpackWallData(char *(&recvBuffer), const int &recvBufSize,
                              int &position, const MPI_Datatype &MPI_wallData,
                              const input &inp) {
  // recvBuffer -- buffer to unpack data from
  // recvBufSize -- size of buffer
  // position -- location within buffer

  // unpack forces
  MPI_Unpack(recvBuffer, recvBufSize, &position, &inviscidForce_, 1, MPI_DOUBLE,
             MPI_COMM_WORLD);
  MPI_Unpack(recvBuffer, recvBufSize, &position, &viscousForce_, 1, MPI_DOUBLE,
             MPI_COMM_WORLD);

  // unpack BC surface
  surf_.UnpackBoundarySurface(recvBuffer, recvBufSize, position);
  
  // get bc data from tag
  bcData_ = inp.BCData(surf_.Tag());

  // unpack wall variables
  data_.ClearResize(surf_.NumI(), surf_.NumJ(), surf_.NumK(), 0);
  MPI_Unpack(recvBuffer, recvBufSize, &position, &(*std::begin(data_)),
             data_.Size(), MPI_wallData, MPI_COMM_WORLD);
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

    multiArray3d<wallVars> newVars(surf_.NumI(), surf_.NumJ(), surf_.NumK(), 0);
    newVars.Insert(dir, {data_.Start(dir), data_.PhysEnd(dir)},
                   data_.Slice(dir, {data_.Start(dir), data_.PhysEnd(dir)}));
    newVars.Insert(dir, {data_.PhysEnd(dir), newVars.End(dir)},
                   upper.data_.Slice(dir, {upper.data_.PhysStart(dir),
                                           upper.data_.End(dir)}));
    data_ = newVars;
  }
}

  primative wallData::WallState(const int &ii, const int &jj, const int &kk,
                     const unique_ptr<eos> &eqnState) const {
    return primative(this->WallDensity(ii, jj, kk), this->WallVelocity(),
                    this->WallPressure(ii, jj, kk, eqnState),
                    this->WallTke(ii, jj, kk), this->WallSdr(ii, jj, kk));
  }

  void wallData::Print(ostream &os) const {
    os << "Inviscid Force: " << inviscidForce_ << endl;
    os << "Viscous Force: " << viscousForce_ << endl;
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