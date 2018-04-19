#pragma once

#include "neutrals_diffusion.hxx"

class ParallelNeutrals : public DiffusionNeutrals {
public:
  ParallelNeutrals(Solver *solver, Mesh *mesh, Options *options);
  void setNeutrals(const Field3D &n_n, const Field3D &v_n);
  void setNeutrals(const Field3D &n_n) override {
    throw BoutException("Please also set velocity!");
  }
  void evolve() override;

protected:
  /// parallel neutral momentum (that we are evolving)
  Field3D m_n;
  /// name of the momentum
  std::string momentum_name;
  /// neutral Temperature;
  BoutReal T_n;
  void mnsheath_yup();
};
