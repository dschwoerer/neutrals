/* License: GPLv3+
 *
 * Implementation of the neutrals interface for simple diffusion like neutrals
 */

#pragma once
#include "neutrals.hxx"
#include "cross_section.hxx"

class DiffusionNeutrals: public Neutrals{
public:
  DiffusionNeutrals(Solver * solver, Options * options);
  /// update the rates, and also set the time derivative, assuming it
  /// is using the solver to evolve the neutrals. The EIRENE coupling
  /// might rather use monitors.
  virtual void update() override;
  virtual void setDensityStag(const Field3D & n_stag);
protected:
  /// neutral density (that we are evolving)
  Field3D n_n;
  /// diffusion `constant`
  Field3D D_neutrals;
  /// pointer to staggered density
  const Field3D * n_stag;
  /// pointer to the solver of the calling physics model
  Solver * solver;
  bool doEvolve; ///< are we evolving the neutrals?
  bool equi_rates; ///< are we using the steady state rates?
  bool onlyion; ///< Do we exclude CX and recombination?
  bool is_static; ///< are the neutrals static?
  /// fraction of neutrals that are lost per time unit
  BoutReal n_n_sink;
  /// faction of the target ion flux that is recycled
  BoutReal recycling_fraction;
  /// length over which the neutrals are recycled
  BoutReal recycling_falloff;
  /// to minimum neutral density enforced
  BoutReal lower_density_limit;
  /// routine to update the rates - only to be called if needed
  void calcRates();
  /// routine to set the time derivative
  void evolve();
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
  /// function returning the recycling profile, for a given target
  /// flux
  virtual Field3D recycle(const FieldPerp & flux);
  CrossSection hydrogen;
};
