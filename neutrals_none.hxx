#include "neutrals.hxx"

class NoNeutrals: public Neutrals{
public:
  NoNeutrals(Solver *solver, Mesh *mesh, CrossSection * cs):
    Neutrals(solver, mesh,cs){
    gamma_CX=0;
    gamma_rec=0;
    gamma_ion=0;
  };
  void update() override{};
  void Neutrals::dumpRates(Datafile &dump) override {};
  void Neutrals::dumpMore(Datafile &dump) override {};

};
