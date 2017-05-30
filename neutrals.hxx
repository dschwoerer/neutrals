/* License: GPLv3+
 *
 * Implementation of an interface to the neutral routines to be used
 * in BOUT++.
 * 
 * Should allow to switch between neutral models without changing the
 * code significantly.
 */

#pragma once
#include <bout/solver.hxx>
#include <bout/constants.hxx>

class Unit{
public:
  virtual BoutReal getDensity() const =0;
  virtual BoutReal getSpeed() const =0;
  virtual BoutReal getLength() const =0;
  virtual BoutReal getTime() const =0;
  virtual BoutReal getEnergy() const =0;
  virtual BoutReal getTemperature() const =0;
};

class SIUnit: public Unit{
public:
  BoutReal getDensity()const {return 1;};
  BoutReal getSpeed()const {return 1;};
  BoutReal getLength()const {return 1;};
  BoutReal getTime()const {return 1;};
  BoutReal getEnergy() const {return 1;};
  BoutReal getTemperature() const {return 1;};
};

class BohmUnit: public Unit{
public:
  BohmUnit(BoutReal B_,BoutReal T_,BoutReal n_,int m_i_):
    B(B_),T(T_),n(n_),m_i(m_i_){};
  BoutReal getDensity() const {return n;};
  BoutReal getSpeed() const {return sqrt(T*SI::qe/(m_i*SI::Mp));};
  BoutReal getLength() const {return getSpeed()*getTime();};
  BoutReal getTime() const {return (m_i*SI::Mp)/(B*SI::qe);};
  BoutReal getEnergy() const {return T*SI::qe;};
  BoutReal getTemperature() const {return T*SI::qe;};
private:
  BohmUnit()=delete;
  BoutReal B;
  BoutReal T;
  BoutReal n;
  BoutReal m_i;
};

class Neutrals{
public:
  Neutrals(Solver * solver_, Mesh * mesh_);
  virtual ~Neutrals(){};
  /// what fields are needed may depend on the neutral model used
  virtual void setPlasmaDensity(const Field3D &n);
  virtual void setElectronTemperature(const Field3D &Te);
  virtual void setIonTemperature(const Field3D &Ti);
  virtual void setElectronVelocity(const Field3D &U);
  virtual void setIonVelocity(const Field3D &U);
  virtual void setPotential(const Field3D &phi);
  // Set the unit system
  virtual void setUnit(const Unit &unit);
  /// dump infos
  virtual void dumpRates(Datafile & dump);
  /// Update rates. If needed, this also evolves the neutrals by
  /// e.g. setting the time derivative of the needed neutral
  /// variables.
  virtual void update()=0;
  /// return the information about neutrals, needed by the plasma code
  /// maybe rather provide getter for momenta & velocities and
  /// densities for electrons and ions?
  /// Might be easier for the user and also provide a more general
  /// interface ...
  /// Derived classes might miss to override all routines, which is
  /// most likely needed, if any is overwritten ...
  virtual const Field3D & getCXRate() const;
  virtual const Field3D & getRecombinationRate() const;
  virtual const Field3D & getIonisationRate() const;
  
  // virtual Field3D getCXOverN() const {
  //   ASSERT2(gamma_CX != nullptr);
  //   if (gamma_CX_over_n == nullptr){
  //     return *gamma_CX/ *n;
  //   } else {
  //     return * gamma_CX_over_n;
  //   }
  // }
  // virtual Field3D getRecOverN() const {
  //   ASSERT2(gamma_rec != nullptr);
  //   if (gamma_rec_over_n == nullptr){
  //     return *gamma_rec/ *n;
  //   } else {
  //     return * gamma_rec_over_n;
  //   }
  // }
  // virtual Field3D getIonOverN() const {
  //   ASSERT2(gamma_ion != nullptr);
  //   if (gamma_ion_over_n == nullptr){
  //     return *gamma_ion/ *n;
  //   } else {
  //     return * gamma_ion_over_n;
  //   }
  // }
  
  /// Source terms (if these terms are actually sinks, the will have a
  /// negative sign)
  //virtual Field3D getIonMomentumSource() const;
  //virtual Field3D getElectronMomentumSource() const;
  virtual Field3D getIonVelocitySource() const;
  virtual Field3D getElectronVelocitySource() const;
  virtual Field3D getDensitySource() const;
  //virtual Field3D getElectronEnergySource() const;
  virtual Field3D getVorticitySource() const;
  virtual Field3D getElectronTemperatureSource() const;
  
protected:
  const Field3D * n; ///< density
  const Field3D * Te; ///< Electron temperature
  const Field3D * Ti; ///< Ion Temperature
  const Field3D * Ui; ///< Ion velocity
  const Field3D * Ue; ///< Electron velocity
  const Field3D * phi; ///< Electrostatic potential
  Field3D gamma_CX; ///< charge exchange rate
  Field3D gamma_ion; ///< ionisation rate
  Field3D gamma_rec; ///< recombination rate
  const Unit * unit;
  BoutReal mu;
  Mesh * mesh;
  Solver * solver;
};


class NeutralsFactory{
public:
  static   std::unique_ptr<Neutrals> create(Solver * solver, Mesh * mesh, Options * option);
  static   std::unique_ptr<Neutrals> create(Solver * solver, Mesh * mesh, std::string option_section);
private:
  NeutralsFactory();
};
