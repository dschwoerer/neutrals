#include "neutrals.hxx"

class NoNeutrals: public Neutrals{
public:
  NoNeutrals(Solver *solver, Mesh *mesh, CrossSection * cs, Options * options):
    Neutrals(solver, mesh,cs, options){
    gamma_CX=0;
    gamma_rec=0;
    gamma_ion=0;
  };
  void update() override{};
  void dumpRates(Datafile &dump) override {};
  void dumpMore(Datafile &dump) override {};

};
