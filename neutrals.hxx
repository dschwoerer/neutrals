/* License: GPLv3+
 *
 * Implementation of an interface to the neutral routines to be used
 * in BOUT++.
 *
 * Should allow to switch between neutral models without changing the
 * code significantly.
 */

#pragma once

#include <string>
#include <bout_types.hxx>
#include <field3d.hxx>

class Solver;
class Mesh;
class Datafile;
class Unit;
class CrossSection;
class Options;
/// The Neutrals class
///
/// Allow to enable or switch between neutral models without changing
/// the rest of physics module code significantly.
///
/// Typical usage requrires to include the neutrals header file:
///
///    #include "neutrals/neutrals.hxx"
///    #include "neutrals/unit.hxx"
///
/// NB: The second include is currently required, as there is not yet
/// a parser to construct the unit section from BOUT.inp. This will
/// hopefully soon be fixed. See issue #2
///
/// Further the neutrals object needs to be stored - the easiest way
/// is to have the pointer be part of the physics class:
///
///    std::unique_ptr<Neutrals> neutrals;
///
/// During initialisation of the Neutrals need to be created. The
/// constructor expects a pointer to the solver, the mesh, and ether
/// the absolute name of the section from the input file where the
/// options for the neutrals are, or a pointer to that section:
///
///     neutrals = NeutralsFactory::create(solver,mesh,"neutrals");
///
/// After that the unit still needs to be set - see issue #2
///
///     BohmUnit * unit = new BohmUnit(B_0,T_e0/e,n_0,m_i/u);
///
/// Further the neutrals need to have certain information about the
/// plasma. What exactly is needed depends on the neutrals
/// model. Further, the staggered quantites are only needed if the
/// velocities are staggered.
///
///     neutrals->setPlasmaDensity(n);
///     neutrals->setIonTemperature(T);
///     neutrals->setElectronTemperature(T);
///     neutrals->setPotential(phi);
///     neutrals->setIonVelocity(U);
///     neutrals->setElectronVelocity(V);
///     neutrals->setPlasmaDensityStag(n_stag);
///
/// After setting the needed a call to init is required, to finalize
/// the initalisation.
///
///     neutrals->init()
///
/// If the neutrals are used during initialisation, e.g. to calculate
/// the potential, the neutrals need to be updated before any
/// quantities of the neutrals can be used. Otherwise it is sufficient
/// to call update in the rhs, after all the fields are valid, and
/// before any quantites from the neutrals are used.
/// After that, the neutrals can be used in the equations of the
/// system, e.g.:
///
///     ddt(vort) += neutrals->getVorticitySource();
///     ddt(n) = - bracket(phi, n, bm, CELL_CENTRE)
///              + neutrals->getDensitySource()
///              .... ;
///
/// Then the neutrals can be configured, using the BOUT++ input file.
/// The section for the neutrals must have the name passed to the
/// factory during the creation of the neutrals object. Thus it is
/// possible to have different section for e.g. different neutrals
/// models, and switch between them.
/// The section can look like this:
///
///     // BOUT.inp section
///     [neutrals]
///     type=parallel
///     use_log_n = true
///     momentum_name = m_n
///
/// For a list of the possible options see the specific models.
/// The fields evolved by the neutrals are reading the boundary
/// conditions etc from the input file. They can be specified as for
/// any other field:
///
///     [neutral_density]
///     function = -10
///     bndry_all = neumann_o2(0)
///
///     [m_n]
///     bndry_all             = neumann_o2(0)
///     bndry_yup             = dirichlet_o4(0)
///
class Neutrals {
public:
  /// Provide the plasma density
  virtual void setPlasmaDensity(const Field3D &n);
  /// Provide the staggered plasma density
  virtual void setPlasmaDensityStag(const Field3D &n_stag);
  /// Provide the electron temperature
  virtual void setElectronTemperature(const Field3D &Te);
  /// Provide the electron temperature. If this is not set, the
  /// electron temperature is assumed.
  virtual void setIonTemperature(const Field3D &Ti);
  /// Provide the parallel electron velocity
  virtual void setElectronVelocity(const Field3D &U);
  /// Provide the parallel ion velocity
  virtual void setIonVelocity(const Field3D &U);
  /// Provide the potential
  virtual void setPotential(const Field3D &phi);
  /// Set the unit system
  virtual void setUnit(const Unit &unit);
  /// scale the extra neutral sources - useful for PID controller
  virtual void scaleSource(BoutReal fac);
  /// Init everything
  virtual void init() {};
  /// dump infos
  virtual void dumpRates(Datafile &dump);
  /// dump more infos
  virtual void dumpMore(Datafile & dump);
  /// Update rates. If needed, this also evolves the neutrals by
  /// e.g. setting the time derivative of the needed neutral
  /// variables.
  virtual void update() = 0;
  /// Source terms (if these terms are actually sinks, the will have a
  /// negative sign)
  /// 
  virtual Field3D getIonVelocitySource() const;
  virtual Field3D getElectronVelocitySource() const;
  virtual Field3D getDensitySource() const;
  // virtual Field3D getElectronEnergySource() const;
  virtual Field3D getVorticitySource() const;
  virtual Field3D getElectronTemperatureSource() const;

  virtual const CrossSection * getCrossSection() const;

  std::string type;

  virtual ~Neutrals();
protected:
  Neutrals(Solver *solver, Mesh *mesh, CrossSection * cs, Options * options);
  void updateMore();
  bool dump_more;
  const Field3D *n;   ///< density
  const Field3D *n_stag; ///< density staggered
  const Field3D *Te;  ///< Electron temperature
  const Field3D *Ti;  ///< Ion Temperature
  const Field3D *Ui;  ///< Ion velocity
  const Field3D *Ue;  ///< Electron velocity
  const Field3D *phi; ///< Electrostatic potential
  Field3D gamma_CX;   ///< charge exchange rate
  Field3D gamma_ion;  ///< ionisation rate
  Field3D gamma_rec;  ///< recombination rate
  Field3D * ionVelocitySource, * densitySource, * electronTemperatureSource;
  const Unit *unit;
  BoutReal mu;
  Mesh *mesh;
  Solver *solver;
  CrossSection * hydrogen;
  Options * options;
};

class NeutralsFactory {
public:
  static std::unique_ptr<Neutrals> create(Solver *solver, Mesh *mesh, Options *option);
  static std::unique_ptr<Neutrals> create(Solver *solver, Mesh *mesh,
                                          std::string option_section);

private:
  NeutralsFactory();
};
