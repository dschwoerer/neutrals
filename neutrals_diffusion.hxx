/* License: GPLv3+
 *
 * Implementation of the neutrals interface for simple diffusion like neutrals
 */

#pragma once
#include "cross_section.hxx"
#include "datafile.hxx"
#include "neutrals.hxx"

class DiffusionNeutrals : public Neutrals {
public:
  DiffusionNeutrals(Solver *solver, Mesh *mesh, CrossSection * cs, Options *options);
  /// update the rates, and also set the time derivative, assuming it
  /// is using the solver to evolve the neutrals. The EIRENE coupling
  /// might rather use monitors.
  virtual void update() override;
  virtual void dumpRates(Datafile &dump);
  virtual void setNeutrals(const Field3D &n_n);
  virtual void scaleSource(const BoutReal) override;
  virtual void init() override;
protected:
  /// neutral density (that we are evolving)
  Field3D n_n;
  /// log of neutral density - only set if evolved
  Field3D l_n_n;
  ///
  bool use_log_n;
  /// name of the density
  std::string density_name;
  /// diffusion `constant`
  Field3D D_neutrals;
  /// recycling flux
  Field3D S_recyc;
  bool doEvolve;   ///< are we evolving the neutrals?
  bool equi_rates; ///< are we using the steady state rates?
  bool onlyion;    ///< Do we exclude CX and recombination?
  bool is_static;  ///< are the neutrals static?
  /// fraction of neutrals that are lost per time unit
  BoutReal loss_fraction;
  /// faction of the target ion flux that is recycled
  BoutReal recycling_fraction;
  /// length over which the neutrals are recycled
  BoutReal recycling_falloff;
  /// to minimum neutral density enforced
  BoutReal lower_density_limit;
  /// to maximum neutral density enforced
  BoutReal higher_density_limit;
  /// routine to update the rates - only to be called if needed
  void calcRates();
  /// set the boundaries
  virtual void setBC();
  /// routine to set the time derivative
  virtual void evolve();
  /// Calculate the diffusion constant
  void calcDiffusion();
  // /// field containing the CX rates
  // Field3D gamma_CX;
  // /// field containing the recombination rates
  // Field3D gamma_rec;
  // /// field containing the ionisation rates
  // Field3D gamma_ion;
  /// fields before multipling with density (if applicable to the model)
  Field3D gamma_CX_over_n;
  Field3D gamma_ion_over_n;
  Field3D gamma_rec_over_n;
  Field3D S0_extra, S_extra;
  /// function returning the recycling profile, for a given target
  /// flux
  virtual Field3D recycle(const FieldPerp &flux);
  void nnsheath_yup();
  /// neutral Temperature;
  BoutReal T_n;
  /// thermal velocity;
  BoutReal v_thermal;
};
