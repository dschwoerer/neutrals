#pragma once

#include "neutrals_diffusion.hxx"

class ParallelNeutrals : public DiffusionNeutrals {
public:
  ParallelNeutrals(Solver *solver, Mesh *mesh, CrossSection * cs, Options *options);
  void setNeutrals(const Field3D &n_n, const Field3D &v_n);
  void setNeutrals(const Field3D &n_n) override {
    throw BoutException("Please also set velocity!");
  }
  void evolve() override;
  virtual Field3D getIonVelocitySource() const override;
  virtual Field3D getElectronVelocitySource() const override;

protected:
  /// parallel neutral momentum (that we are evolving)
  Field3D m_n;
  Field3D v_n;
  /// name of the momentum
  std::string momentum_name;
  /// neutral Temperature;
  BoutReal T_n;
  void mnsheath_yup();
  void nnsheath_ydown();
  virtual void setBC() override;
};
