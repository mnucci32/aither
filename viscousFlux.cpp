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

#include <cmath>  // sqrt
#include <string>
#include "viscousFlux.hpp"
#include "eos.hpp"         // idealGas
#include "primVars.hpp"    // primVars
#include "turbulence.hpp"  // turbModel
#include "matrix.hpp"      // squareMatrix

using std::cout;
using std::endl;
using std::cerr;
using std::vector;
using std::string;
using std::max;

// constructor -- initialize flux from velocity gradient
/*
Viscous flux normal to face:
F = [ 0,
      taux,
      tauy,
      tauz,
      tau (dot) vel + K * tGrad (dot) area
      (mu + mut) * tkeGrad (dot) area
      (mu + mut) * omegaGrad (dot) area ]

In the above equation tau is the wall shear stress. Taux, tauy, and tauz are the
rows of the wall shear stress tensor i.e. taux = tauxx + tauxy + tauxz. K is the
thermal conductivity, tGrad is the temperature gradient, and area is the
normalized face area.

Wall shear stress:
tau = lambda * velGradTrace * area + mu * ( velGrad * area + velGrad' * area)

In the above equation lambda is the bulk viscosity, velGradTrace is the trace of
the velocity gradient, area is the normalized face area, mu is the dynamic
viscosity, and velGrad is the velocity gradient tensor.
*/
viscousFlux::viscousFlux(
    const tensor<double> &velGrad, const vector3d<double> &vel,
    const double &mu, const double &eddyVisc, const sutherland &suth,
    const idealGas &eqnState, const vector3d<double> &tGrad,
    const vector3d<double> &areaVec, const vector3d<double> &tkeGrad,
    const vector3d<double> &omegaGrad, const turbModel *turb,
    const primVars &state) {
  // velGrad -- velocity gradient tensor
  // vel -- velocity vector
  // mu -- dynamic viscosity
  // eddyVisc -- turbulent eddy viscosity
  // suth -- method to get viscosity as a function of temperature (Sutherland's
  // law)
  // eqnState -- equation of state
  // tGrad -- temperature gradient
  // areaVec -- area vector of face
  // tkeGrad -- tke gradient
  // omegaGrad -- omega gradient
  // turb -- turbulence model
  // state -- primative variables at face

  vector3d<double> normArea = areaVec / areaVec.Mag();  // normalize face area

  // get 2nd coefficient of viscosity assuming bulk viscosity is 0 (Stoke's)
  double lambda = suth.Lambda(mu + eddyVisc);

  double velGradTrace = velGrad.Trace();  // trace of velocity gradient
  // wall shear stress
  vector3d<double> tau =
      lambda * velGradTrace * normArea +
      (mu + eddyVisc) *
          (velGrad.MatMult(normArea) + velGrad.Transpose().MatMult(normArea));

  data_[0] = tau.X();
  data_[1] = tau.Y();
  data_[2] = tau.Z();
  data_[3] = tau.DotProd(vel) + (eqnState.Conductivity(mu) +
                                 eqnState.TurbConductivity(
                                     eddyVisc, turb->TurbPrandtlNumber())) *
                                 tGrad.DotProd(normArea);

  // turbulence viscous flux
  data_[4] = (mu + suth.NondimScaling() * turb->MolecDiff1Coeff() *
              turb->EddyViscNoLim(state)) * tkeGrad.DotProd(normArea);

  data_[5] = (mu + suth.NondimScaling() * turb->MolecDiff2Coeff() *
              turb->EddyViscNoLim(state)) * omegaGrad.DotProd(normArea);
}

// non-member functions
// ----------------------------------------------------------------------------
// operator overload for << - allows use of cout, cerr, etc.
ostream &operator<<(ostream &os, viscousFlux &flux) {
  os << "0.0, ";
  for (int ii = 0; ii < NUMVARS - 1; ii++) {
    os << flux.data_[ii];
    if (ii != NUMVARS - 2) {
      os << ", ";
    }
  }
  return os;
}

// member function for scalar multiplication
viscousFlux viscousFlux::operator*(const double &scalar) const {
  viscousFlux temp = *this;
  for (int ii = 0; ii < NUMVARS - 1; ii++) {
    temp.data_[ii] *= scalar;
  }
  return temp;
}

// friend function to allow multiplication from either direction
viscousFlux operator*(const double &scalar, const viscousFlux &flux) {
  viscousFlux temp;
  for (int ii = 0; ii < NUMVARS - 1; ii++) {
    temp.data_[ii] = flux.data_[ii] * scalar;
  }
  return temp;
}

// member function for scalar division
viscousFlux viscousFlux::operator/(const double &scalar) const {
  viscousFlux temp = *this;
  for (int ii = 0; ii < NUMVARS - 1; ii++) {
    temp.data_[ii] /= scalar;
  }
  return temp;
}

// friend function to allow division from either direction
viscousFlux operator/(const double &scalar, const viscousFlux &flux) {
  viscousFlux temp;
  for (int ii = 0; ii < NUMVARS - 1; ii++) {
    temp.data_[0] = scalar / flux.data_[ii];
  }
  return temp;
}

// function to calculate the thin shear layer flux jacobian -- NOT USED in LUSGS
// formulation
void CalcTSLFluxJac(const double &mu, const double &eddyVisc,
                    const idealGas &eqnState, const vector3d<double> &areaVec,
                    const primVars &left, const primVars &right,
                    const double &dist, squareMatrix &dFv_dUl,
                    squareMatrix &dFv_dUr, const sutherland &suth,
                    const double &prt) {
  // mu -- dynamic viscosity
  // eddyVisc -- turbulent eddy viscosity
  // eqnState -- equation of state
  // areaVec -- area vector of face
  // left -- left state (primative)
  // right -- right state (primative)
  // dist -- distance from centroid of left cell to centroid of right cell
  // dFv_dUl -- flux jacobian of viscous flux with respect to left state
  // dFV_dUr -- flux jacobian of viscous flux with respect to right state
  // suth -- method to get viscosity as a function of temperature

  // check to make sure square matrices are of correct size
  if (dFv_dUl.Size() != 5 || dFv_dUr.Size() != 5) {
    cerr << "ERROR: Error in viscousFlux.cpp:CalcTSLFluxJac. Problem with thin "
            "shear layer viscous jacobian calculation. The input jacobian "
            "matrices are not of the correct size!" << endl;
    exit(0);
  }

  // normalize area vector
  vector3d<double> normArea = areaVec / areaVec.Mag();

  // get velocity at face
  vector3d<double> vel = 0.5 * (right.Velocity() + left.Velocity());

  // calculate thin shear layer velocity gradients
  tensor<double> velGradTSL = CalcVelGradTSL(left, right, normArea, dist);

  // calculate bulk viscosity
  double lambda = suth.Lambda(mu + eddyVisc);

  // calculate shear stress at face
  double velGradTrace = velGradTSL.Trace();
  vector3d<double> tau =
      lambda * velGradTrace * normArea +
      (mu + eddyVisc) * (velGradTSL.MatMult(normArea) +
                         velGradTSL.Transpose().MatMult(normArea));

  // calculate coefficients (from Blazek)
  double theta = normArea.MagSq();
  double thetaX = (4.0 / 3.0) * normArea.X() * normArea.X() +
                  normArea.Y() * normArea.Y() + normArea.Z() * normArea.Z();
  double thetaY = normArea.X() * normArea.X() +
                  (4.0 / 3.0) * normArea.Y() * normArea.Y() +
                  normArea.Z() * normArea.Z();
  double thetaZ = normArea.X() * normArea.X() + normArea.Y() * normArea.Y() +
                  (4.0 / 3.0) * normArea.Z() * normArea.Z();

  double etaX = (1.0 / 3.0) * normArea.Y() * normArea.Z();
  double etaY = (1.0 / 3.0) * normArea.X() * normArea.Z();
  double etaZ = (1.0 / 3.0) * normArea.X() * normArea.Y();

  double piX = vel.X() * thetaX + vel.Y() * etaZ + vel.Z() * etaY;
  double piY = vel.X() * etaZ + vel.Y() * thetaY + vel.Z() * etaX;
  double piZ = vel.X() * etaY + vel.Y() * etaX + vel.Z() * thetaZ;

  double phiRhoL = -1.0 * (eqnState.Conductivity(mu) +
                           eqnState.TurbConductivity(eddyVisc, prt)) *
                   left.Temperature(eqnState) / ((mu + eddyVisc) * left.Rho());
  double phiRhoR = -1.0 * (eqnState.Conductivity(mu) +
                           eqnState.TurbConductivity(eddyVisc, prt)) *
                   right.Temperature(eqnState) /
                   ((mu + eddyVisc) * right.Rho());

  double phiPressL = (eqnState.Conductivity(mu) +
                      eqnState.TurbConductivity(eddyVisc, prt)) /
                     ((mu + eddyVisc) * left.Rho());
  double phiPressR = (eqnState.Conductivity(mu) +
                      eqnState.TurbConductivity(eddyVisc, prt)) /
                     ((mu + eddyVisc) * right.Rho());

  // calculate matrix - derivative of left primative vars wrt left conservative
  // vars
  squareMatrix dWl_dUl(5);
  dWl_dUl.Zero();

  // column 0
  dWl_dUl.SetData(0, 0, 1.0);
  dWl_dUl.SetData(1, 0, -1.0 * left.U() / left.Rho());
  dWl_dUl.SetData(2, 0, -1.0 * left.V() / left.Rho());
  dWl_dUl.SetData(3, 0, -1.0 * left.W() / left.Rho());
  dWl_dUl.SetData(4, 0,
                  0.5 * (eqnState.Gamma() - 1.0) * left.Velocity().MagSq());

  // column 1
  dWl_dUl.SetData(1, 1, 1.0 / left.Rho());
  dWl_dUl.SetData(4, 1, -1.0 * (eqnState.Gamma() - 1.0) * left.U());

  // column 2
  dWl_dUl.SetData(2, 2, 1.0 / left.Rho());
  dWl_dUl.SetData(4, 2, -1.0 * (eqnState.Gamma() - 1.0) * left.V());

  // column 3
  dWl_dUl.SetData(3, 3, 1.0 / left.Rho());
  dWl_dUl.SetData(4, 3, -1.0 * (eqnState.Gamma() - 1.0) * left.W());

  // column 4
  dWl_dUl.SetData(4, 4, eqnState.Gamma() - 1.0);

  //-------------------------------------------------------------------------
  // calculate matrix - derivative of right primative vars wrt right
  // conservative vars
  squareMatrix dWr_dUr(5);
  dWr_dUr.Zero();

  // column 0
  dWr_dUr.SetData(0, 0, 1.0);
  dWr_dUr.SetData(1, 0, -1.0 * right.U() / right.Rho());
  dWr_dUr.SetData(2, 0, -1.0 * right.V() / right.Rho());
  dWr_dUr.SetData(3, 0, -1.0 * right.W() / right.Rho());
  dWr_dUr.SetData(4, 0,
                  0.5 * (eqnState.Gamma() - 1.0) * right.Velocity().MagSq());

  // column 1
  dWr_dUr.SetData(1, 1, 1.0 / right.Rho());
  dWr_dUr.SetData(4, 1, -1.0 * (eqnState.Gamma() - 1.0) * right.U());

  // column 2
  dWr_dUr.SetData(2, 2, 1.0 / right.Rho());
  dWr_dUr.SetData(4, 2, -1.0 * (eqnState.Gamma() - 1.0) * right.V());

  // column 3
  dWr_dUr.SetData(3, 3, 1.0 / right.Rho());
  dWr_dUr.SetData(4, 3, -1.0 * (eqnState.Gamma() - 1.0) * right.W());

  // column 4
  dWr_dUr.SetData(4, 4, eqnState.Gamma() - 1.0);

  //------------------------------------------------------------------------
  // calculate matrix - derivative of viscous flux wrt left primative vars
  // column 0
  dFv_dUl.SetData(0, 0, 0.0);
  dFv_dUl.SetData(1, 0, 0.0);
  dFv_dUl.SetData(2, 0, 0.0);
  dFv_dUl.SetData(3, 0, 0.0);
  dFv_dUl.SetData(4, 0, phiRhoL * theta);

  // column 1
  dFv_dUl.SetData(0, 1, 0.0);
  dFv_dUl.SetData(1, 1, thetaX);
  dFv_dUl.SetData(2, 1, etaZ);
  dFv_dUl.SetData(3, 1, etaY);
  dFv_dUl.SetData(4, 1, -0.5 * (dist / (mu + eddyVisc)) * tau.X() + piX);

  // column 2
  dFv_dUl.SetData(0, 2, 0.0);
  dFv_dUl.SetData(1, 2, etaZ);
  dFv_dUl.SetData(2, 2, thetaY);
  dFv_dUl.SetData(3, 2, etaX);
  dFv_dUl.SetData(4, 2, -0.5 * (dist / (mu + eddyVisc)) * tau.Y() + piY);

  // column 3
  dFv_dUl.SetData(0, 3, 0.0);
  dFv_dUl.SetData(1, 3, etaY);
  dFv_dUl.SetData(2, 3, etaX);
  dFv_dUl.SetData(3, 3, thetaZ);
  dFv_dUl.SetData(4, 3, -0.5 * (dist / (mu + eddyVisc)) * tau.Z() + piZ);

  // column 4
  dFv_dUl.SetData(0, 4, 0.0);
  dFv_dUl.SetData(1, 4, 0.0);
  dFv_dUl.SetData(2, 4, 0.0);
  dFv_dUl.SetData(3, 4, 0.0);
  dFv_dUl.SetData(4, 4, phiPressL * theta);

  dFv_dUl = -1.0 * ((mu + eddyVisc) / dist) * dFv_dUl;

  //-----------------------------------------------------------------------
  // calculate matrix - derivative of viscous flux wrt right primative vars
  // column 0
  dFv_dUr.SetData(0, 0, 0.0);
  dFv_dUr.SetData(1, 0, 0.0);
  dFv_dUr.SetData(2, 0, 0.0);
  dFv_dUr.SetData(3, 0, 0.0);
  dFv_dUr.SetData(4, 0, phiRhoR * theta);

  // column 1
  dFv_dUr.SetData(0, 1, 0.0);
  dFv_dUr.SetData(1, 1, thetaX);
  dFv_dUr.SetData(2, 1, etaZ);
  dFv_dUr.SetData(3, 1, etaY);
  dFv_dUr.SetData(4, 1, 0.5 * (dist / (mu + eddyVisc)) * tau.X() + piX);

  // column 2
  dFv_dUr.SetData(0, 2, 0.0);
  dFv_dUr.SetData(1, 2, etaZ);
  dFv_dUr.SetData(2, 2, thetaY);
  dFv_dUr.SetData(3, 2, etaX);
  dFv_dUr.SetData(4, 2, 0.5 * (dist / (mu + eddyVisc)) * tau.Y() + piY);

  // column 3
  dFv_dUr.SetData(0, 3, 0.0);
  dFv_dUr.SetData(1, 3, etaY);
  dFv_dUr.SetData(2, 3, etaX);
  dFv_dUr.SetData(3, 3, thetaZ);
  dFv_dUr.SetData(4, 3, 0.5 * (dist / (mu + eddyVisc)) * tau.Z() + piZ);

  // column 4
  dFv_dUr.SetData(0, 4, 0.0);
  dFv_dUr.SetData(1, 4, 0.0);
  dFv_dUr.SetData(2, 4, 0.0);
  dFv_dUr.SetData(3, 4, 0.0);
  dFv_dUr.SetData(4, 4, phiPressR * theta);

  dFv_dUr = ((mu + eddyVisc) / dist) * dFv_dUr;

  // multiply by dW_dU to get flux jacobian derivative wrt conservative
  // variables
  dFv_dUl = dFv_dUl * dWl_dUl;
  dFv_dUr = dFv_dUr * dWr_dUr;

  // calculate spectral radius
  primVars faceState = 0.5 * (left + right);
  dFv_dUl.Identity();
  dFv_dUr.Identity();
  double specRad = (mu + eddyVisc) * eqnState.Gamma() /
                   (eqnState.Prandtl() * faceState.Rho() * dist);

  // add or subtract spectral radius to flux jacobian
  dFv_dUl = -1.0 * specRad * dFv_dUl;
  dFv_dUr = specRad * dFv_dUr;
}

// function to calculate the velocity gradients at a cell face using the Thin
// Shear Layer approximation
// NOT USED in LUSGS formulation
tensor<double> CalcVelGradTSL(const primVars &left, const primVars &right,
                              const vector3d<double> &areaVec,
                              const double &dist) {
  // left -- left state (primative)
  // right -- right state (primative)
  // areaVec -- area vector of face
  // dist -- distance between centroid of left cell and right cell

  // normalize area vector
  vector3d<double> normArea = areaVec / areaVec.Mag();

  // initialize velocity gradient tensor
  tensor<double> velGrad;

  // calculate velocity derivatives
  vector3d<double> velDeriv = (right.Velocity() - left.Velocity()) / dist;

  // populate velocity gradient tensor
  velGrad.SetXX(velDeriv.X() * normArea.X());
  velGrad.SetXY(velDeriv.Y() * normArea.X());
  velGrad.SetXZ(velDeriv.Z() * normArea.X());

  velGrad.SetYX(velDeriv.X() * normArea.Y());
  velGrad.SetYY(velDeriv.Y() * normArea.Y());
  velGrad.SetYZ(velDeriv.Z() * normArea.Y());

  velGrad.SetZX(velDeriv.X() * normArea.Z());
  velGrad.SetZY(velDeriv.Y() * normArea.Z());
  velGrad.SetZZ(velDeriv.Z() * normArea.Z());

  return velGrad;
}
