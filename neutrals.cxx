#include "neutrals.hxx"
#include "neutrals_diffusion.hxx"
#include "interpolation.hxx"
#include "git_version.hxx"

Neutrals::Neutrals(Solver * solver_, Mesh * mesh_):
  n(nullptr), Te(nullptr), Ti(nullptr), Ui(nullptr), Ue(nullptr),
  phi(nullptr),
  unit(nullptr), mu(1850),
  solver(solver_),mesh(mesh_){
  output.write("************************************"
               "**********************************\n") ;
  output.write("\tNeutrals-API Version %s\n",NEUTRALS_GIT_SHA1);
  output.write("************************************"
               "**********************************\n") ;
}

void Neutrals::setPlasmaDensity(const Field3D &n_){
  n=&n_;
}
void Neutrals::setElectronTemperature(const Field3D &Te_){
  Te=&Te_;
}

void Neutrals::setIonTemperature(const Field3D &Ti_){
  Ti=&Ti_;
}

void Neutrals::setIonVelocity(const Field3D &U_){
  Ui=&U_;
}
void Neutrals::setElectronVelocity(const Field3D &U_){
  Ue=&U_;
}

void Neutrals::dumpRates(Datafile & dump){
  SAVE_REPEAT(gamma_CX);
  SAVE_REPEAT(gamma_rec);
  SAVE_REPEAT(gamma_ion);
}

void Neutrals::setUnit(const Unit &unit_){
  unit=&unit_;
}

void Neutrals::setPotential(const Field3D &phi_){
  phi=&phi_;
}

const Field3D & Neutrals::getCXRate() const{
  return gamma_CX;
}

const Field3D & Neutrals::getRecombinationRate() const{
  return gamma_rec;
}

const Field3D & Neutrals::getIonisationRate() const{
  return gamma_ion;
}

/// this code is not staggering aware
// Field3D Neutrals::getIonMomentumSource() const{
//   ASSERT2(Ui        != nullptr);
//   ASSERT2(gamma_CX  != nullptr);
//   ASSERT2(gamma_rec != nullptr);
//   return -(*Ui)*((*gamma_CX)+(*gamma_rec));
// }
  
// Field3D Neutrals::getElectronMomentumSource() const{
//   ASSERT2(Ue        != nullptr);
//   ASSERT2(gamma_rec != nullptr);
//   return -(*Ue)*((*gamma_rec));
// }

Field3D Neutrals::getIonVelocitySource() const{
  ASSERT2(Ui        != nullptr);
  Field3D tmp=gamma_CX+gamma_ion;
  tmp/=*n;//getCXOverN()+getIonOverN();
  return -(*Ui)*(interp_to(tmp,Ui->getLocation()));
}

  
Field3D Neutrals::getElectronVelocitySource() const{
  ASSERT2(Ue        != nullptr);
  Field3D tmp=gamma_ion/(*n);
  return -(*Ue)*(interp_to(tmp,Ue->getLocation()));
}

Field3D Neutrals::getDensitySource() const{
  return gamma_ion - gamma_rec;
}

Field3D Neutrals::getElectronTemperatureSource() const{
  ASSERT2(Ue != nullptr);
  ASSERT2(Te != nullptr);
  ASSERT2(mu > 0);
  Field3D rec_over_n=gamma_rec/(*n);
  Field3D ion_over_n=gamma_ion/(*n);
  return rec_over_n*interp_to(SQ(*Ue),Te->getLocation())/(3.*mu)
    + (*Te)*(rec_over_n-ion_over_n);
}

Field3D Neutrals::getVorticitySource() const{
  ASSERT2(phi != nullptr);
  return -Delp2(*phi)*(gamma_CX+gamma_ion)
    -Grad_perp(*phi)*Grad_perp(gamma_CX+gamma_ion);
}
// Field3D Neutrals::


//   virtual Field3D getIonMomentumSource() const;
//   virtual Field3D getElectronMomentumSource() const;
//   virtual Field3D getIonVelocitySource() const;
//   virtual Field3D getElectronVelocitySource() const;
//   virtual Field3D getDensitySource() const;
//   virtual Field3D getElectronEnergySource() const;
//   virtual Field3D getElectronTemperatureSource() const;

std::unique_ptr<Neutrals>
NeutralsFactory::create(Solver * solver, Mesh * mesh, Options * options){
  std::string type;
  OPTION(options, type, "NotSet");
  if (type == "diffusion"){
    return std::unique_ptr<Neutrals>(new DiffusionNeutrals(solver, mesh, options));
  } else {
    throw BoutException("unknow neutrals model '%s'",type.c_str());
  }
}

std::unique_ptr<Neutrals>
NeutralsFactory::create(Solver * solver, Mesh * mesh, std::string options){
  return create(solver, mesh, Options::getRoot()->getSection(options));
}
