#ifndef WREMNANTS_CSVARIABLES_H
#define WREMNANTS_CSVARIABLES_H


#include "Math/GenVector/Boost.h"
#include "Math/Vector4D.h"
#include "Math/Vector3D.h"
#include "TLorentzVector.h"
#include <ROOT/RVec.hxx>
#include <iostream>

namespace wrem {

typedef ROOT::Math::PxPyPzEVector PxPyPzEVector;
typedef ROOT::Math::PxPyPzMVector PxPyPzMVector;
typedef ROOT::Math::PtEtaPhiMVector PtEtaPhiMVector;
typedef ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag> Vector;

template <typename T>
double dot(const T& vec1, const T& vec2) {
    return vec1.x()*vec2.x()+vec1.y()*vec2.y()+vec1.z()*vec2.z();
}

template <typename T>
T cross(const T& vec1, const T& vec2) {
    auto cross = ROOT::Math::Cross<double>({vec1.x(), vec1.y(), vec1.z()}, {vec2.x(), vec2.y(), vec2.z()});
    T res(cross[0], cross[1], cross[2]);
    return res;
}

Vector unitBoostedVector(ROOT::Math::Boost& boostOp, PxPyPzEVector& vec) {
    PxPyPzEVector boostvec = boostOp(vec);
    return Vector(boostvec.x(), boostvec.y(), boostvec.z()).Unit();
}

struct CSVars {

  double sintheta;
  double costheta;
  double sinphi;
  double cosphi;

};

CSVars csSineCosThetaPhi(const PtEtaPhiMVector& antilepton, const PtEtaPhiMVector& lepton) {
    PxPyPzEVector antilepton_v(antilepton);
    PxPyPzEVector dilepton = antilepton_v + PxPyPzEVector(lepton);
    const int zsign = std::copysign(1.0, dilepton.z());
    const double energy = 6500.;
    PxPyPzEVector proton1(0., 0., zsign*energy, energy);
    PxPyPzEVector proton2(0., 0., -1.*zsign*energy, energy);

    auto dilepCM = dilepton.BoostToCM();
    ROOT::Math::Boost dilepCMBoost(dilepCM);

    auto pro1boost = unitBoostedVector(dilepCMBoost, proton1);
    auto pro2boost = unitBoostedVector(dilepCMBoost, proton2);
    auto antilepton_boost = unitBoostedVector(dilepCMBoost, antilepton_v);
    auto csFrame = (pro1boost-pro2boost).Unit();
    auto csYaxis = cross(pro1boost, pro2boost).Unit();
    auto csXaxis = cross(csYaxis, csFrame).Unit();

    double costheta = dot(csFrame, antilepton_boost);
    auto csCross = cross(csFrame, antilepton_boost);
    double sintheta = csCross.R()/(csFrame.R()*antilepton_boost.R());

    double sinphi = dot(csYaxis, antilepton_boost)/sintheta;
    double cosphi = dot(csXaxis, antilepton_boost)/sintheta;

    CSVars angles = {sintheta, costheta, sinphi, cosphi};
    return angles;
}

}

#endif
