#pragma once

#include "radiation.hxx"
#include <field3d.hxx>

// Field3D crosssection_CX(const Field3D &v_cx);

class CrossSection {
public:
  Field3D cx_rate(const Field3D &T_e);
  BoutReal cx_rate(const BoutReal T_e);

  Field3D ionisation_rate(const Field3D &T_e);
  BoutReal ionisation_rate(const BoutReal T_e);

  Field3D recombination_rate(const Field3D &n, const Field3D &T_e);
  Field3D recombination_rate(const Field3D &n, const BoutReal T_e);
  Field3D recombination_rate(const BoutReal n, const Field3D &T_e);
  BoutReal recombination_rate(const BoutReal n, const BoutReal T_e);

  // constructor;
  CrossSection(RadiatedPower *atom);

private:
  // CrossSection(const CrossSection &); // disable copying

  RadiatedPower &atom;
};
