#include "neutrals_diffusion.hxx"

#include "bout/constants.hxx"
#include "bout/surfaceiter.hxx"
#include "field_factory.hxx"
#include "helper.hxx"
#include "interpolation.hxx"
#include "radiation.hxx"
#include "radiation_factory.hxx"
#include "unit.hxx"
#include <bout/solver.hxx>

// n_n sheath boundary condition
void DiffusionNeutrals::nnsheath_yup() {
  for (int x = 0; x < mesh->LocalNx; ++x)
    for (int z = 0; z < mesh->LocalNz; ++z) {
      int y = mesh->yend;
      for (int dy = mesh->ystart - 1; dy >= 0; --dy)
        n_n(x, y + 1 + dy, z) = n_n(x, y - dy, z);
    }
  int myg = mesh->ystart;
  for (int x = 0; x < mesh->LocalNx; ++x)
    for (int y = mesh->ystart - 1; y >= 0; --y)
      for (int z = 0; z < mesh->LocalNz; ++z) {
        n_n(x, y, z) = n_n(x, 2 * myg - 1 - y, z);
        ;
      }
}

DiffusionNeutrals::DiffusionNeutrals(Solver *solver, Mesh *mesh, CrossSection *cs,
                                     Options *options)
    : Neutrals(solver, mesh, cs, options), T_n(nan("")), v_thermal(nan("")) {
  OPTION(options, equi_rates, false);
  OPTION(options, recycling_falloff, 4.0);     // m
  OPTION(options, lower_density_limit, 8e10);  // in m^-3
  OPTION(options, higher_density_limit, 4e19); // in m^-3
  OPTION(options, recycling_fraction, 0.9);    // fraction of target flux that is recycled
  OPTION(options, loss_fraction, 1e-5); // fraction of neutrals that is lost per time
  OPTION(options, diffusion_factor, 10.);
  options->get("evolve", doEvolve, true);
  OPTION(options, onlyion, false);
  OPTION(options, use_log_n, false);
  auto extra = options->getSection("density_source");

  if (!extra->isSet("function")) {
    S_extra = 0;
  } else {
    std::string function;
    OPTION(extra, function, "0");
    output.write("\tFound extra functions: %s\n", function.c_str());
    S_extra = FieldFactory::get()->create3D(function, extra, mesh, CELL_CENTRE, 0);
    output.write("\tFunction is in the range [%e,%e]\n", min(S_extra, true),
                 max(S_extra, true));
  }
  S0_extra = S_extra;

  OPTION(options, impurity_fraction, 0);
  if (impurity_fraction > 0) {
    std::string impurity_model_name;
    options->get("impurity_model", impurity_model_name, "hutchinsoncarbonradiation");
    impurity_model = new CrossSection(RadiatedPowerFactory::create(impurity_model_name));
  }

  if (equi_rates && doEvolve) {
    throw BoutException(
        "DiffusionNeutrals:: cannot have equilibrium rates with evolving neutrals!");
  }
  if (doEvolve) {
    n_n = 1e-4;
    OPTION(options, density_name, "neutral_density");
    if (use_log_n) {
      l_n_n = log(n_n);
      solver->add(l_n_n, density_name.c_str());
    } else {
      solver->add(n_n, density_name.c_str());
    }
  }
}

void DiffusionNeutrals::init() {
  // 300 K in units of 40eV
  std::string temperature_unit;
  OPTION(options, temperature_unit, "kb");
  BoutReal t_unit;
  if (temperature_unit == "kb" || temperature_unit == "k") {
    t_unit = SI::kb;
  } else if (temperature_unit == "ev") {
    t_unit = SI::qe;
  } else if (temperature_unit == "default") {
    t_unit = 1;
  } else {
    throw BoutException("Unknown unit - %s", temperature_unit.c_str());
  }
  // 300 K in default unit
  BoutReal val = 300 * SI::kb / t_unit; // /( SI::kb / 40 / SI::qe);
  // printf("%g\n", val);
  OPTION(options, T_n, val);
  T_n /= unit->getTemperature() / t_unit;
  // T_n /= t_unit;
  v_thermal =
      sqrt(2 * T_n * unit->getTemperature() / (SI::M_Deuterium)) / unit->getSpeed();
  output_info.write("\tNeutrals temperature: %.10e eV\n",
                    T_n * unit->getTemperature() / SI::qe);
  output_info.write("\tNeutrals temperature: %.10e kB\n",
                    T_n * unit->getTemperature() / SI::kb);
  output_info.write("\tNeutrals thermal velocity: %.10e\n", v_thermal);
  output_info.write("\tNeutrals thermal velocity: %.10e m/s\n",
                    v_thermal * unit->getSpeed());

  recycling_dist = Field2D(0., n->getMesh());
  if (recycling_falloff) {
    BoutReal falloff = recycling_falloff / unit->getLength();
    Field2D intDy =
        FieldFactory::get()->create2D("realy", nullptr, n->getMesh(), CELL_CENTRE);
    SurfaceIter s(mesh);
    Indices i;
    for (s.first(); !s.isDone(); s.next()) {
      if (s.closed()) {
        continue; // Skip if closed flux surface
      }
      i.x = s.xpos;
      int nproc;
      MPI_Comm comm = s.communicator();
      MPI_Comm_size(comm, &nproc);

      // Get the position of the target
      BoutReal L = intDy(i.x, mesh->yend); // Cell centre
      L += mesh->coordinates()->dy(i.x, mesh->yend) / 2;
      MPI_Bcast(&L, 1, MPI_DOUBLE, nproc - 1, comm);

      // Get the non-normalised values
      BoutReal sum = 0;
      for (i.y = 0; i.y < mesh->LocalNy; ++i.y) {
        recycling_dist[i] = exp(-SQ(L - intDy[i]) / (2 * SQ(falloff)));
        if (i.y >= mesh->ystart && i.y <= mesh->yend) {
          sum += recycling_dist[i];
        }
      }
      // normalise
      MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, comm);
      for (i.y = 0; i.y < mesh->LocalNy; ++i.y) {
        recycling_dist[i] /= sum * mesh->coordinates()->dy[i];
      }
    } // end iteration
  }   // end if recycling fallof
}

void DiffusionNeutrals::scaleSource(BoutReal fac) { S_extra = S0_extra * (fac + 1); }

void DiffusionNeutrals::setNeutrals(const Field3D &n_n_) {
  if (use_log_n) {
    l_n_n = n_n_;
    n_n=exp(l_n_n);
  } else {
    n_n = n_n_;
  }
}

void DiffusionNeutrals::dumpRates(Datafile &dump) {
  if (doEvolve) {
    SAVE_REPEAT(D_neutrals);
    SAVE_REPEAT(S_recyc);
  }
  if (!equi_rates) {
    SAVE_REPEAT(gamma_CX);
    SAVE_REPEAT(gamma_ion);
    SAVE_REPEAT(gamma_rec);
  }
}

void DiffusionNeutrals::update() {
  if (!equi_rates) {
    this->setBC();
    this->calcRates();
    if (doEvolve) {
      this->evolve();
    }
  }
  updateMore();
}

void DiffusionNeutrals::setBC() {
  nnsheath_yup();
  if (doEvolve) {
    mesh->communicate(n_n);
  }
}

void DiffusionNeutrals::evolve() {
  if (n_stag == nullptr) {
    if (n->getLocation() != Ui->getLocation()) {
      throw BoutException("DiffusionNeutrals:: density and velocity at different "
                          "location, but staggered density not given!");
    } else {
      n_stag = n;
    }
  }
  if (use_log_n) {
    n_n = exp(l_n_n);
  }
  FieldPerp flux = sliceXZ(*Ui, mesh->yend + 1) * sliceXZ(*n_stag, mesh->yend + 1);
  S_recyc = recycle(flux);
  calcDiffusion();
  limit_at_least_smooth(D_neutrals, 1e2);
  ddt(n_n) = (+gamma_rec - gamma_ion + S_recyc - n_n * loss_fraction +
              D_neutrals * Laplace(n_n) + S_extra);
  if (use_log_n) {
    ddt(l_n_n) = 1 / n_n * ddt(n_n);
  }

  // ddt(n_n)  += Grad(D_neutrals) * Grad(n_n);
  for (int x = 0; x < mesh->LocalNx; ++x) {
    for (int y = 0; y < mesh->LocalNy; ++y) {
      if (y == mesh->ystart)
        y = mesh->yend + 1;
      for (int z = 0; z < mesh->LocalNz; ++z) {
        ddt(n_n)(x, y, z) = 0;
      }
    }
  }
}

void DiffusionNeutrals::calcDiffusion() {
  // compute D - taken from Bens sim-cat model
  // thermal velocity:
  // http://www.wolframalpha.com/input/?i=sqrt%282*%20300+K+*k_B++%2F%282+u%29%29&a=*MC.K+!*k!_B-_*Unit-
  // 1579 thermal speed of deuterium in m/s @ 300 K
  // sqrt(2 * 300 * SI::kb / 2 / SI::u)
  const BoutReal a0 = PI * SQ(5.29e-11 / unit->getLength()); // normalised
  const BoutReal fac =
      (v_thermal * a0 * (unit->getDensity() * pow(unit->getLength(), 3)));
  Field3D sigma_nn = fac * n_n;
  D_neutrals = SQ(v_thermal) / (sigma_nn + gamma_CX + gamma_ion);
  D_neutrals *= diffusion_factor;
}

void DiffusionNeutrals::calcRates() {
  ASSERT2(unit != nullptr);
  ASSERT2(Ti != nullptr);
  ASSERT2(Te != nullptr);
  ASSERT2(Ui != nullptr);
  if (lower_density_limit > 0) {
    limit_at_least(n_n, lower_density_limit / unit->getDensity());
    limit_at_most(n_n, higher_density_limit / unit->getDensity());
  }
  BoutReal m_i = 2 * SI::Mp;
  Field3D Ti_in_eV = (*Ti) * (unit->getTemperature() / SI::qe);
  Field3D Te_in_eV;
  if (Ti == Te) {
    Te_in_eV = Ti_in_eV;
  } else {
    Te_in_eV = (*Te) * (unit->getTemperature() / SI::qe);
  }
  Field3D nnn0oOmega = n_n * (unit->getDensity() * unit->getTime());
  gamma_ion_over_n = nnn0oOmega * hydrogen->ionisation_rate(Te_in_eV);
  gamma_ion = gamma_ion_over_n * (*n);
  if (!onlyion) {
    Field3D energy = SQ(interp_to(*Ui, CELL_CENTER)) *
                     (m_i / 2 * unit->getSpeed() * unit->getSpeed() / SI::qe); // in eV
    energy += Ti_in_eV;
    gamma_CX_over_n = hydrogen->cx_rate(energy);
    gamma_CX_over_n *= nnn0oOmega;
    gamma_CX = gamma_CX_over_n * (*n);

    gamma_rec_over_n = (*n) *
                       hydrogen->recombination_rate((*n) * unit->getDensity(), Te_in_eV) *
                       (unit->getDensity() * unit->getTime()); // guard cells not set
    gamma_rec = (*n) * gamma_rec_over_n;
  } else {
    gamma_CX = 0;
    gamma_CX_over_n = 0;
    gamma_rec = 0;
    gamma_rec_over_n = 0;
  }
  radiation_loss = 0;
  if (impurity_fraction > 0) {
    radiation_loss +=
        impurity_model->power((*Ti) * (unit->getTemperature() / SI::qe),
                              (*n) * unit->getDensity(),
                              (*n) * (unit->getDensity() * impurity_fraction)) /
        (unit->getTemperature() * unit->getDensity() / unit->getTime());
  }
}

Field3D DiffusionNeutrals::recycle(const FieldPerp &flux) {
  Field3D result;
  result.allocate();
#if CHECK > 1
  for (auto i : result) {
    result[i] = nan("");
  }
#endif
  // Set a recycling source in first and last processors
  // BoutReal L;	//Parallel length coordinate
  // BoutReal TargetFluxInner,TargetFluxOuter;
  SurfaceIter s(mesh);
  MPI_Comm comm;
  int nproc;
  for (s.first(); !s.isDone(); s.next()) {
    int i = s.xpos;
    if (s.closed()) {
      // todo: change loops
      for (int j = 0; j < mesh->LocalNy; j++) {
        for (int k = 0; k < mesh->LocalNz; k++) {
          result(i, j, k) = 0;
        }
      }
      continue; // Skip if closed flux surface
    }
    Coordinates *coord = mesh->coordinates();
    if (recycling_falloff) {
      comm = s.communicator();
      MPI_Comm_size(comm, &nproc);
      BoutReal TargetFluxOuter[mesh->LocalNz];
      for (int k = 0; k < mesh->LocalNz; k++) {
        // TargetFluxOuter = Flux(i,mesh->yend+1,k)/mesh->Bxy(i,mesh->yend+1);
        // todo: replace abs with min(,0)?
        TargetFluxOuter[k] = abs(flux(i, k) * recycling_fraction);
        ASSERT2(TargetFluxOuter[k] != 0)
      }
      MPI_Bcast(TargetFluxOuter, mesh->LocalNz, MPI_DOUBLE, nproc - 1, comm);
      // printf("%d %g\n",mesh->GlobalNy,tmpsum);
      for (int j = 0; j < mesh->LocalNy; j++) {
        for (int k = 0; k < mesh->LocalNz; k++) {
          // if(TargetFluxOuter[k] == 0){throw BoutException("maeh");}
          result(i, j, k) = recycling_dist(i, j) * TargetFluxOuter[k];
        }
      }
    } else {
      for (int j = 0; j < mesh->LocalNy; j++) {
        for (int k = 0; k < mesh->LocalNz; k++) {
          result(i, j, k) = 0;
        }
      }
      if (s.lastY()) {
        for (int k = 0; k < mesh->LocalNz; k++) {
          result(i, mesh->yend, k) = flux(i, k) / coord->dy(i, mesh->yend);
        }
      }
    }
  }
  // limit_at_least(result,1e-30,true);
  return result;
}

Field3D DiffusionNeutrals::getElectronTemperatureSource() const {
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
  result += (T_n - (*Te)) * ion_over_n;
  // Charge exchange cooles the plasma as well
  result += (T_n - (*Te)) * gamma_CX / (*n);
  // Radiation terms
  // Based on Ben's SD1D model

  result += -(1.09 * (*Te) - 13.6 * (SI::qe / unit->getTemperature())) * rec_over_n -
            30 * (SI::qe / unit->getTemperature()) * ion_over_n;

  result -= 2. / 3. * radiation_loss / (*n);

  return result;
}
