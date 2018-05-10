#include "interpolation.hxx"

#include "helper.hxx"
#include "neutrals_parallel.hxx"
#include "derivs.hxx"
#include "unit.hxx"
#include "bout/solver.hxx"

ParallelNeutrals::ParallelNeutrals(Solver *solver, Mesh *mesh, CrossSection * cs, Options *options)
  : DiffusionNeutrals(solver, mesh, cs, options) {
  // 300 K in units of 40eV
  BoutReal val = 300 * SI::kb / 40 / SI::qe;
  //printf("%g\n", val);
  OPTION(options, T_n, val);
  OPTION(options, momentum_name, "m_n");
  if (doEvolve) {
    solver->add(m_n, momentum_name.c_str());
  }
}

void ParallelNeutrals::setNeutrals(const Field3D &n_n_, const Field3D &m_n_) {
  n_n = n_n_;
  m_n = m_n_;
}

void ParallelNeutrals::evolve() {
  if (n_stag == nullptr) {
    if (n->getLocation() != Ui->getLocation()) {
      throw BoutException("DiffusionNeutrals:: density and velocity at different "
                          "location, but staggered density not given!");
    } else {
      n_stag = n;
    }
  }
  CELL_LOC mylow = Ui->getLocation();
  //ASSERT1(mylow==CELL_YLOW);
  m_n.setLocation(mylow);
  n_n.applyBoundary();
  m_n.applyBoundary();
  //nnsheath_yup();
  //mnsheath_yup();
  //nnsheath_ydown();
  FieldPerp flux = sliceXZ(*Ui, mesh->yend + 1) * sliceXZ(*n_stag, mesh->yend + 1);
  S_recyc = recycle(flux);
  ddt(n_n) = - Div_par(m_n, CELL_CENTRE);
  //ddt(n_n) += gamma_rec - gamma_ion + S_recyc;
  //ddt(n_n) += -n_n * loss_fraction;
  
  //ddt(n_n) -= 10*D4DY4(n_n);

  auto v_n = m_n / interp_to(n_n, mylow);
  auto tmp = gamma_CX + gamma_ion;
  //ddt(m_n) = 0;
  ddt(m_n) = -Vpar_Grad_par(v_n, m_n);
  //ddt(m_n) += -Vpar_Grad_par(m_n, v_n);
  ddt(m_n) += - Grad_par(n_n * T_n, mylow);
  //ddt(m_n) += - 10*D4DY4(m_n);
    //+((*Ui) - v_n) * interp_to(tmp, mylow);

  // compute D - taken from Bens sim-cat model
  // thermal velocity:
  if (false) {
    const BoutReal thermal_speed_neut =
      sqrt(2 * T_n * unit->getTemperature() / (2 * SI::Mp));
    const BoutReal a0 = PI * SQ(5.29e-11 / unit->getLength()); // normalised
    const BoutReal fac =
      (thermal_speed_neut * a0 * (unit->getDensity() * pow(unit->getLength(), 3)));
    Field3D sigma_nn = fac * n_n;
    D_neutrals = SQ(thermal_speed_neut) / (sigma_nn + gamma_CX + gamma_ion);
    limit_at_least_smooth(D_neutrals, 1e2);
    ddt(n_n) += D_neutrals * Laplace(n_n);
  } else {
    D_neutrals = 0;
  }
  
  // // ddt(n_n)  += Grad(D_neutrals) * Grad(n_n);
  // for (int x = 0; x < mesh->LocalNx; ++x) {
  //   for (int y = 0; y < mesh->LocalNy; ++y) {
  //     if (y == mesh->ystart)
  //       y = mesh->yend + 1;
  //     for (int z = 0; z < mesh->LocalNz; ++z) {
  //       ddt(n_n)(x, y, z) = 0;
  //     }
  //   }
  // }
}

void ParallelNeutrals::mnsheath_yup() {
  for (int x = 0; x < mesh->LocalNx; ++x) {
    // for (int y=mesh->yend+1;y<mesh->LocalNy;++y){
    int y = mesh->yend + 1;
    for (int z = 0; z < mesh->LocalNz; ++z) {
      m_n(x, y, z) = 0;
    }
    y += 1;
    if (y < mesh->LocalNy) {
      for (int z = 0; z < mesh->LocalNz; ++z) {
        m_n(x, y, z) = -m_n(x, y - 2, z);
      }
    }
  }
}

void ParallelNeutrals::nnsheath_ydown() {
  for (int x = 0; x < mesh->LocalNx; ++x) {
    int y = 1;
    for (int z = 0; z < mesh->LocalNz; ++z) {
      n_n(x, 1, z) = 1;
    }
    if (y < mesh->LocalNy) {
      for (int z = 0; z < mesh->LocalNz; ++z) {
        n_n(x, 0, z) = 1-n_n(x,  2, z);
      }
    }
  }
}
