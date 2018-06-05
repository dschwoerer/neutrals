#include "neutrals.hxx"
#include "git_version.hxx"
#include "interpolation.hxx"
#include "cross_section.hxx"
#include "neutrals_diffusion.hxx"
#include "neutrals_parallel.hxx"
#include "neutrals_none.hxx"

#include <difops.hxx>

Neutrals::Neutrals(Solver *solver, Mesh *mesh, CrossSection *cs)
  : dump_more(false), n(nullptr), n_stag(nullptr), Te(nullptr), Ti(nullptr), Ui(nullptr), Ue(nullptr), phi(nullptr),
    unit(nullptr), mu(1850), solver(solver), mesh(mesh), hydrogen(cs) {
  output.write("************************************"
               "**********************************\n");
  output.write("\tNeutrals-API Version %s\n", NEUTRALS_GIT_SHA1);
  output.write("************************************"
               "**********************************\n");
}

Neutrals::~Neutrals() {
  delete hydrogen;
};

void Neutrals::setPlasmaDensity(const Field3D &n_) { n = &n_; }
void Neutrals::setElectronTemperature(const Field3D &Te_) { Te = &Te_; }

void Neutrals::setIonTemperature(const Field3D &Ti_) { Ti = &Ti_; }

void Neutrals::setIonVelocity(const Field3D &U_) { Ui = &U_; }
void Neutrals::setElectronVelocity(const Field3D &U_) { Ue = &U_; }

void Neutrals::setPlasmaDensityStag(const Field3D &n_stag_) {
  n_stag = &n_stag_;
}

const CrossSection * Neutrals::getCrossSection() const {
  return hydrogen;
}

void Neutrals::dumpRates(Datafile &dump) {
  SAVE_REPEAT(gamma_CX);
  SAVE_REPEAT(gamma_rec);
  SAVE_REPEAT(gamma_ion);
}

void Neutrals::dumpMore(Datafile &dump) {
  dumpRates(dump);
  dump_more = true;
  ionVelocitySource = new Field3D;
  densitySource = new Field3D();;
  electronTemperatureSource = new Field3D;
  dump.add(*ionVelocitySource,"ionVelocitySource",true);
  dump.add(*densitySource,"densitySource",true);
  dump.add(*electronTemperatureSource,"electronTemperatureSource",true);
}

void Neutrals::updateMore() {
  if (dump_more) {
    *ionVelocitySource = getIonVelocitySource();
    *densitySource = getDensitySource();
    *electronTemperatureSource = getElectronTemperatureSource();
  }
}


void Neutrals::setUnit(const Unit &unit_) { unit = &unit_; }

void Neutrals::setPotential(const Field3D &phi_) { phi = &phi_; }

const Field3D &Neutrals::getCXRate() const { return gamma_CX; }

const Field3D &Neutrals::getRecombinationRate() const { return gamma_rec; }

const Field3D &Neutrals::getIonisationRate() const { return gamma_ion; }

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

Field3D Neutrals::getIonVelocitySource() const {
  ASSERT2(Ui != nullptr);
  Field3D tmp = gamma_CX + gamma_ion;
  tmp /= *n; // getCXOverN()+getIonOverN();
  return -(*Ui) * (interp_to(tmp, Ui->getLocation()));
}

Field3D Neutrals::getElectronVelocitySource() const {
  ASSERT2(Ue != nullptr);
  Field3D tmp = gamma_ion / (*n);
  return -(*Ue) * (interp_to(tmp, Ue->getLocation()));
}

Field3D Neutrals::getDensitySource() const { return gamma_ion - gamma_rec; }

Field3D Neutrals::getElectronTemperatureSource() const {
  ASSERT2(Ue != nullptr);
  ASSERT2(Te != nullptr);
  ASSERT2(mu > 0);
  Field3D rec_over_n = gamma_rec / (*n);
  Field3D ion_over_n = gamma_ion / (*n);
  return rec_over_n * interp_to(SQ(*Ue), Te->getLocation()) / (3. * mu) +
         (*Te) * (rec_over_n - ion_over_n);
}

Field3D Neutrals::getVorticitySource() const {
  ASSERT2(phi != nullptr);
  return -Delp2(*phi) * (gamma_CX + gamma_ion) -
         Grad_perp(*phi) * Grad_perp(gamma_CX + gamma_ion);
}
// Field3D Neutrals::

//   virtual Field3D getIonMomentumSource() const;
//   virtual Field3D getElectronMomentumSource() const;
//   virtual Field3D getIonVelocitySource() const;
//   virtual Field3D getElectronVelocitySource() const;
//   virtual Field3D getDensitySource() const;
//   virtual Field3D getElectronEnergySource() const;
//   virtual Field3D getElectronTemperatureSource() const;

std::unique_ptr<Neutrals> NeutralsFactory::create(Solver *solver, Mesh *mesh,
                                                  Options *options) {
  std::string type;
  std::unique_ptr<Neutrals> ret;
  OPTION(options, type, "NotSet");
  CrossSection * cs = CrossSectionFactory::create(options);
  if (type == "diffusion") {
    ret = std::unique_ptr<Neutrals>(new DiffusionNeutrals(solver, mesh, cs, options));
  } else if (type == "parallel") {
    ret = std::unique_ptr<Neutrals>(new ParallelNeutrals(solver, mesh, cs, options));
  } else if (type == "none") {
    ret = std::unique_ptr<Neutrals>(new NoNeutrals(solver, mesh,cs));
  } else {
    delete cs;
    throw BoutException("unknow neutrals model '%s'", type.c_str());
  }
  ret->type = type;
  output.write("************************************"
               "**********************************\n");
  return ret;
}

std::unique_ptr<Neutrals> NeutralsFactory::create(Solver *solver, Mesh *mesh,
                                                  std::string options) {
  return create(solver, mesh, Options::getRoot()->getSection(options));
}
