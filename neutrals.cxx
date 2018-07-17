#include "neutrals.hxx"
#include "bout/constants.hxx"
#include "cross_section.hxx"
#include "git_version.hxx"
#include "interpolation.hxx"
#include "neutrals_diffusion.hxx"
#include "neutrals_none.hxx"
#include "neutrals_parallel.hxx"
#include "unit.hxx"

#include <difops.hxx>

Neutrals::Neutrals(Solver *solver, Mesh *mesh, CrossSection *cs, Options *options)
    : dump_more(false), n(nullptr), n_stag(nullptr), Te(nullptr), Ti(nullptr),
      Ui(nullptr), Ue(nullptr), phi(nullptr), unit(nullptr), mu(1850), solver(solver),
      mesh(mesh), hydrogen(cs), options(options) {
  output.write("************************************"
               "**********************************\n");
  output.write("\tNeutrals-API Version %s\n", NEUTRALS_GIT_SHA1);
  output.write("************************************"
               "**********************************\n");
}

Neutrals::~Neutrals() { delete hydrogen; };

void Neutrals::setPlasmaDensity(const Field3D &n_) { n = &n_; }
void Neutrals::setElectronTemperature(const Field3D &Te_) { Te = &Te_; }

void Neutrals::setIonTemperature(const Field3D &Ti_) { Ti = &Ti_; }

void Neutrals::setIonVelocity(const Field3D &U_) { Ui = &U_; }
void Neutrals::setElectronVelocity(const Field3D &U_) { Ue = &U_; }

void Neutrals::setPlasmaDensityStag(const Field3D &n_stag_) { n_stag = &n_stag_; }

void Neutrals::scaleSource(BoutReal fac) { ; }
const CrossSection *Neutrals::getCrossSection() const { return hydrogen; }

void Neutrals::dumpRates(Datafile &dump) {
  SAVE_REPEAT(gamma_CX);
  SAVE_REPEAT(gamma_rec);
  SAVE_REPEAT(gamma_ion);
}

void Neutrals::dumpMore(Datafile &dump) {
  dumpRates(dump);
  dump_more = true;
  ionVelocitySource = new Field3D;
  densitySource = new Field3D();
  ;
  electronTemperatureSource = new Field3D;
  dump.add(*ionVelocitySource, "ionVelocitySource", true);
  dump.add(*densitySource, "densitySource", true);
  dump.add(*electronTemperatureSource, "electronTemperatureSource", true);
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

/*const Field3D &Neutrals::getCXRate() const { return gamma_CX; }

const Field3D &Neutrals::getRecombinationRate() const { return gamma_rec; }

const Field3D &Neutrals::getIonisationRate() const { return gamma_ion; }
*/

// Friction term between neutrals and ions
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
  Field3D result = rec_over_n * interp_to(SQ(*Ue), Te->getLocation()) / (3. * mu);
  // We do not need this term, as we are asking for a change in temperature.
  // The particles that go, change the pressure, but also the density
  // - while keeping the temperature constant.
  // result += (*Te) * (rec_over_n );
  // Ionized particles have the temperature of the background gas.
  // Thus they contribute the temperature difference times the rates
  // over n (i.e. relative importance)
  result += (-(*Te)) * ion_over_n;
  // Charge exchange cooles the plasma as well
  result += (-(*Te)) * gamma_CX / (*n);

  result += -(1.09 * (*Te) - 13.6 * SI::qe / unit->getTemperature()) * rec_over_n -
            30 * SI::qe / unit->getTemperature() * ion_over_n;

  return result;
}

Field3D Neutrals::getVorticitySource() const {
  ASSERT2(phi != nullptr);
  return -Delp2(*phi) * (gamma_CX + gamma_ion) -
         Grad_perp(*phi) * Grad_perp(gamma_CX + gamma_ion);
}

std::unique_ptr<Neutrals> NeutralsFactory::create(Solver *solver, Mesh *mesh,
                                                  Options *options) {
  std::string type;
  std::unique_ptr<Neutrals> ret;
  OPTION(options, type, "NotSet");
  CrossSection *cs = CrossSectionFactory::create(options);
  if (type == "diffusion") {
    ret = std::unique_ptr<Neutrals>(new DiffusionNeutrals(solver, mesh, cs, options));
  } else if (type == "parallel") {
    ret = std::unique_ptr<Neutrals>(new ParallelNeutrals(solver, mesh, cs, options));
  } else if (type == "none") {
    ret = std::unique_ptr<Neutrals>(new NoNeutrals(solver, mesh, cs, options));
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
