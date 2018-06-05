#pragma once

#include "radiation.hxx"
#include "radiation_factory.hxx"
#include <field3d.hxx>
#include <output.hxx>

// Field3D crosssection_CX(const Field3D &v_cx);

class CrossSection {
public:
  virtual Field3D cx_rate(const Field3D &T_e);
  //BoutReal cx_rate(const BoutReal T_e);

  virtual Field3D ionisation_rate(const Field3D &T_e);
  //BoutReal ionisation_rate(const BoutReal T_e);

  virtual Field3D recombination_rate(const Field3D &n, const Field3D &T_e);
  // Field3D recombination_rate(const Field3D &n, const BoutReal T_e);
  // Field3D recombination_rate(const BoutReal n, const Field3D &T_e);
  // BoutReal recombination_rate(const BoutReal n, const BoutReal T_e);

  // constructor;
  CrossSection(RadiatedPower *atom);

  ~CrossSection() {
    delete atom;
  }
private:
  // CrossSection(const CrossSection &); // disable copying

  RadiatedPower * atom;
};


class TestCrossSection: public CrossSection{
public:
  virtual Field3D cx_rate(const Field3D &T_e) override {
    return mycx;
  }
  virtual Field3D ionisation_rate(const Field3D &T_e) override {
    return myion;
  }
  virtual Field3D recombination_rate(const Field3D &n, const Field3D &T_e) override {
    return myrec;
  }

  Field3D mycx, myion, myrec;

  TestCrossSection(RadiatedPower *atom=nullptr):CrossSection(atom){
  };
};


class CrossSectionFactory {
public:
  static CrossSection * create(Options * opt){
    std::string type;
    opt->get("CrossSection",type,"radiation");
    if (type == "radiation") {
      return new CrossSection(RadiatedPowerFactory::create(opt));
    } else if (type == "testing") {
      return new TestCrossSection(RadiatedPowerFactory::create(opt));
    } else {
      throw BoutException("Unknown CrossSection requested for CrossSectionFactory::create - '%s'",type.c_str());
    }
  }
};
