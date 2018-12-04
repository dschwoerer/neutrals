#ifndef __RADIATION_H__
#define __RADIATION_H__

#include <cmath>
#include <vector>

#include <bout_types.hxx>
#include <field3d.hxx>
#include <options.hxx>

class RadiatedPower {
public:
  const Field3D power(const Field3D &Te, const Field3D &Ne, const Field3D &Ni);

  virtual BoutReal power(BoutReal Te, BoutReal ne, BoutReal ni) {
    throw BoutException("power not implemented!");
  }

  virtual BoutReal ionisation(BoutReal Te) {
    throw BoutException("ionisation not implemented!");
  }

  virtual BoutReal recombination(BoutReal n, BoutReal Te) {
    throw BoutException("recombination not implemented!");
  }

  virtual BoutReal chargeExchange(BoutReal Te) {
    throw BoutException("chargeExchange not implemented!");
  }

  virtual BoutReal excitation(BoutReal Te) {
    throw BoutException("excitation not implemented!");
  }

private:
};

class InterpRadiatedPower : public RadiatedPower {
public:
  InterpRadiatedPower(const std::string &file);

  BoutReal power(BoutReal Te, BoutReal ne, BoutReal ni);

private:
  std::vector<BoutReal> te_array; // Te in eV
  std::vector<BoutReal> p_array;  // Radiative loss rate in Watts m^3
};

/// Rates supplied by Eva Havlicova
class HydrogenRadiatedPower : public RadiatedPower {
public:
  BoutReal power(BoutReal Te, BoutReal ne, BoutReal ni);

  // Collision rate coefficient <sigma*v> [m3/s]
  BoutReal ionisation(const BoutReal Te);

  //<sigma*v> [m3/s]
  BoutReal recombination(const BoutReal n, const BoutReal Te);

  // <sigma*v> [m3/s]
  BoutReal chargeExchange(const BoutReal Te);

  // <sigma*v> [m3/s]
  BoutReal excitation(const BoutReal Te);

private:
};

/*!
 * Hydrogen rates, fitted by Hannah Willett May 2015
 * University of York
 */
class UpdatedRadiatedPower : public RadiatedPower {
public:
  BoutReal power(BoutReal Te, BoutReal ne, BoutReal ni);

  // Ionisation rate coefficient <sigma*v> [m3/s]
  BoutReal ionisation(BoutReal T);

  // Recombination rate coefficient <sigma*v> [m3/s]
  BoutReal recombination(BoutReal n, BoutReal T);

  // Charge exchange rate coefficient <sigma*v> [m3/s]
  BoutReal chargeExchange(BoutReal Te);

  BoutReal excitation(BoutReal Te);

private:
};

/// Carbon in coronal equilibrium
/// From I.H.Hutchinson Nucl. Fusion 34 (10) 1337 - 1348 (1994)
class HutchinsonCarbonRadiation : public RadiatedPower {
  BoutReal power(BoutReal Te, BoutReal ne, BoutReal ni) {
    return ne * ni * 2e-31 * pow(Te / 10., 3) / (1. + pow(Te / 10., 4.5));
  }
};

/*
class PostJensen : public RadiatedPower {
public:
  BoutReal power(BoutReal Te, BoutReal ne, BoutReal ni) {
    if( (Te < Tmin) || (Te > Tmax) )
      return 0.0;

    return 0.0;
  }
protected:
  BoutReal Tmin, Tmax;
  BoutReal A[6];

  struct PJ_Data {
    const char* label; // Short name
    const char* name;  // Long name
    BoutReal Tmin, Tmax;
    BoutReal data[6];
  };

  static PJ_Data power_data[] = {
    {"C", "Carbon"},
    {0}
  };
};
*/

class TestingPower : public RadiatedPower {
public:
  TestingPower(Options *opt);
  BoutReal power(BoutReal Te, BoutReal ne, BoutReal ni);

  BoutReal ionisation(BoutReal Te);

  BoutReal recombination(BoutReal n, BoutReal Te);

  BoutReal chargeExchange(BoutReal Te);

  BoutReal excitation(BoutReal Te);

private:
  BoutReal _power;
  BoutReal _ionisation;
  BoutReal _recombination;
  BoutReal _chargeExchange;
  BoutReal _excitation;
};

#endif // __RADIATION_H__
