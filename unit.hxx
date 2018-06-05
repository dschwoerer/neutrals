#pragma once

#include <math.h>

#include <bout/constants.hxx>

class Unit {
public:
  virtual BoutReal getDensity() const = 0;
  virtual BoutReal getSpeed() const = 0;
  virtual BoutReal getLength() const = 0;
  virtual BoutReal getTime() const = 0;
  virtual BoutReal getEnergy() const = 0;
  virtual BoutReal getTemperature() const = 0;
};

class SIUnit : public Unit {
public:
  BoutReal getDensity() const { return 1; };
  BoutReal getSpeed() const { return 1; };
  BoutReal getLength() const { return 1; };
  BoutReal getTime() const { return 1; };
  BoutReal getEnergy() const { return 1; };
  BoutReal getTemperature() const { return 1; };
};

class BohmUnit : public Unit {
public:
  BohmUnit(BoutReal B, BoutReal T, BoutReal n, int m_i)
      : B(B), T(T), n(n), m_i(m_i){
    printf("\tBohmUnit: %g %g %g %d\n",B,T,n,m_i);
  };
  BoutReal getDensity() const { return n; };
  BoutReal getSpeed() const { return sqrt(T * SI::qe / (m_i * SI::Mp)); };
  BoutReal getLength() const { return getSpeed() * getTime(); };
  BoutReal getTime() const { return (m_i * SI::Mp) / (B * SI::qe); };
  BoutReal getEnergy() const { return T * SI::qe; };
  BoutReal getTemperature() const { return T * SI::qe; };

private:
  BohmUnit() = delete;
  BoutReal B;
  BoutReal T;
  BoutReal n;
  BoutReal m_i;
};
