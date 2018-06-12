#include "neutrals_diffusion.hxx"

#include "bout/constants.hxx"
#include "bout/surfaceiter.hxx"
#include "helper.hxx"
#include "interpolation.hxx"
#include "radiation.hxx"
#include "radiation_factory.hxx"
#include "unit.hxx"
#include <bout/solver.hxx>
#include "field_factory.hxx"

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

DiffusionNeutrals::DiffusionNeutrals(Solver *solver, Mesh *mesh, CrossSection * cs, Options *options)
  : Neutrals(solver, mesh, cs) {
  OPTION(options, equi_rates, false);
  OPTION(options, recycling_falloff, 4.0);    // m
  OPTION(options, lower_density_limit, 8e10); // in m^-3
  OPTION(options, higher_density_limit, 4e19); // in m^-3
  OPTION(options, recycling_fraction, 0.9);   // fraction of target flux that is recycled
  OPTION(options, loss_fraction, 1e-5); // fraction of neutrals that is lost per time
  options->get("evolve", doEvolve, true);
  OPTION(options, onlyion, false);
  OPTION(options, use_log_n, false);
  auto extra = options->getSection("density_source");
  
  if (!extra->isSet("function")){
    S_extra=0;
  } else {
    std::string function;
    OPTION(extra,function,"0");
    output.write("\tFound extra functions: %s\n",function.c_str());
    S_extra = FieldFactory::get()->create3D("function",extra,mesh, CELL_CENTRE,0);
    output.write("\tFunction is in the range [%e,%e]\n",min(S_extra,true),max(S_extra,true));
  }
  S0_extra=S_extra;

  if (equi_rates && doEvolve) {
    throw BoutException(
        "DiffusionNeutrals:: cannot have equilibrium rates with evolving neutrals!");
  }
  if (doEvolve) {
    n_n = 1e-4;
    OPTION(options, density_name, "neutral_density");
    if (use_log_n) {
      l_n_n=log(n_n);
      solver->add(l_n_n, density_name.c_str());
    } else {
      solver->add(n_n, density_name.c_str());
    }
  } else {
    throw BoutException("we really should get the density!");
  }
}

void DiffusionNeutrals::scaleSource(BoutReal fac) {
  S_extra = S0_extra*fac;
}

void DiffusionNeutrals::setNeutrals(const Field3D &n_n_) {
  n_n = n_n_;
  if (use_log_n) {
    l_n_n = log(n_n);
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
    n_n =  exp(l_n_n);
  }
  FieldPerp flux = sliceXZ(*Ui, mesh->yend + 1) * sliceXZ(*n_stag, mesh->yend + 1);
  S_recyc = recycle(flux);
  calcDiffusion();
  limit_at_least_smooth(D_neutrals, 1e2);
  ddt(n_n) = (
              + gamma_rec - gamma_ion + S_recyc
              - n_n * loss_fraction
              + D_neutrals * Laplace(n_n)
              + S_extra
              );
  if (use_log_n) {
    ddt(l_n_n) = 1/n_n * ddt(n_n);
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
  const BoutReal thermal_speed_neut = 1579 / unit->getSpeed();
  const BoutReal a0 = PI * SQ(5.29e-11 / unit->getLength()); // normalised
  const BoutReal fac =
      (thermal_speed_neut * a0 * (unit->getDensity() * pow(unit->getLength(), 3)));
  Field3D sigma_nn = fac * n_n;
  D_neutrals = SQ(thermal_speed_neut) / (sigma_nn + gamma_CX + gamma_ion);
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
      BoutReal falloff = recycling_falloff / unit->getLength();
      comm = s.communicator();
      MPI_Comm_size(comm, &nproc);
      BoutReal TargetFluxOuter[mesh->LocalNz];
      for (int k = 0; k < mesh->LocalNz; k++) {
        // if(s.lastY())
        // TargetFluxOuter = Flux(i,mesh->yend+1,k)/mesh->Bxy(i,mesh->yend+1);
        // todo: replace abs with min(,0)?
        TargetFluxOuter[k] = abs(flux(i, k) * recycling_fraction);
        // if (TargetFluxOuter[k] == 0){throw BoutException("fuuu");}
      }
      MPI_Bcast(TargetFluxOuter, mesh->LocalNz, MPI_DOUBLE, nproc - 1, comm);
      // if(s.firstY())		//On inner target boundary, broadcast target flux
      // TargetFluxInner = Flux(i,2,k)/mesh->Bxy(i,2);
      // MPI_Bcast(&TargetFluxInner,1,MPI_DOUBLE,0,comm);
      // BoutReal sum=0;
      // BoutReal den;
      BoutReal tmpsum = 0;
      for (int j = 2; j < mesh->GlobalNy - 2; j++) {
        BoutReal L =
            (mesh->GlobalNy - 4.5 - mesh->YGLOBAL(j)) * coord->dy(i, j) *
            sqrt(coord->g_22(i, j)); // Make zero point between guard cell and first cell
        // printf("%d %g  ",mesh->YGLOBAL(j),L);
        tmpsum += (2 / (falloff * sqrt(TWOPI))) * exp(-SQ(L) / (2 * SQ(falloff))) *
                  coord->dy(i, j);
      }
      // printf("%d %g\n",mesh->GlobalNy,tmpsum);
      for (int j = 0; j < mesh->LocalNy; j++) {
        // TODO: BUG: L_para should not depend on j?
        // BoutReal L_para=(mesh->GlobalNy-4)*coord->dy(i,j)*sqrt(coord->g_22(i,j));
        BoutReal L =
            (mesh->GlobalNy - 4.5 - mesh->YGLOBAL(j)) * coord->dy(i, j) *
            sqrt(coord->g_22(i, j)); // Make zero point between guard cell and first cell
        // TODO: replace abs(?) with min(?,0)?
        // BoutReal tmp=1/falloff*exp(-L/falloff)/(1-exp(-L_para/falloff));
        BoutReal tmp = (2 / (falloff * sqrt(TWOPI))) * exp(-SQ(L) / (2 * SQ(falloff)));
        tmp /= tmpsum;
        // if(tmp == 0){throw BoutException("maeh");}
        for (int k = 0; k < mesh->LocalNz; k++) {
          // if(TargetFluxOuter[k] == 0){throw BoutException("maeh");}
          result(i, j, k) = tmp * TargetFluxOuter[k];
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
