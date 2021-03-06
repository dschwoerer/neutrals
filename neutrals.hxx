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

class Neutrals{
public:
  virtual ~Neutrals(){};
  /// what fields are needed may depend on the neutral model used
  virtual void setPlasmaDensity(const Field3D &n);
  virtual void setElectronTemperature(const Field3D &Te);
  virtual void setIonTemperature(const Field3D &Ti);
  virtual void setElectronVelocity(const Field3D &U);
  virtual void setIonVelocity(const Field3D &U);
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
  
  virtual Field3D getCXOverN() const {
    ASSERT2(gamma_CX != nullptr);
    if (gamma_CX_over_n == nullptr){
      return *gamma_CX/ *n;
    } else {
      return * gamma_CX_over_n;
    }
  }
  virtual Field3D getRecOverN() const {
    ASSERT2(gamma_rec != nullptr);
    if (gamma_rec_over_n == nullptr){
      return *gamma_rec/ *n;
    } else {
      return * gamma_rec_over_n;
    }
  }
  virtual Field3D getIonOverN() const {
    ASSERT2(gamma_ion != nullptr);
    if (gamma_ion_over_n == nullptr){
      return *gamma_ion/ *n;
    } else {
      return * gamma_ion_over_n;
    }
  }
  
  /// Source terms (if these terms are actually sinks, the will have a
  /// negative sign)
  //virtual Field3D getIonMomentumSource() const;
  //virtual Field3D getElectronMomentumSource() const;
  virtual Field3D getIonVelocitySource() const;
  virtual Field3D getElectronVelocitySource() const;
  virtual Field3D getDensitySource() const;
  //virtual Field3D getElectronEnergySource() const;
  virtual Field3D getElectronTemperatureSource() const;
protected:
  const Field3D * n; ///< density
  const Field3D * Te; ///< Electron temperature
  const Field3D * Ti; ///< Ion Temperature
  const Field3D * Ui; ///< Ion velocity
  const Field3D * Ue; ///< Electron velocity
  Field3D * gamma_CX; ///< charge exchange rate
  Field3D * gamma_ion; ///< ionisation rate
  Field3D * gamma_rec; ///< recombination rate
  /// fields before multipling with density (if applicable to the model)
  Field3D * gamma_CX_over_n;
  Field3D * gamma_ion_over_n;
  Field3D * gamma_rec_over_n;
  BoutReal mu;
};


class NeutralsFactory{
public:
  static   std::unique_ptr<Neutrals> create(Solver * solver, Options * option);
  static   std::unique_ptr<Neutrals> create(Solver * solver, std::string option_section);
private:
  NeutralsFactory();
};
