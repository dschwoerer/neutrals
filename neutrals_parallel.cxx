#include "interpolation.hxx"

#include "bout/solver.hxx"
#include "derivs.hxx"
#include "helper.hxx"
#include "neutrals_parallel.hxx"
#include "unit.hxx"

ParallelNeutrals::ParallelNeutrals(Solver *solver, Mesh *mesh, CrossSection *cs,
                                   Options *options)
    : DiffusionNeutrals(solver, mesh, cs, options) {

  OPTION(options, momentum_name, "m_n");
  if (doEvolve) {
    // m_n.setLocation(CELL_YLOW);
    solver->add(m_n, momentum_name.c_str());
  }
}

void ParallelNeutrals::setNeutrals(const Field3D &n_n_, const Field3D &m_n_) {
  if (use_log_n) {
    l_n_n = n_n_;
  } else {
    n_n = n_n_;
  }
  m_n = m_n_;
}

void ParallelNeutrals::setBC() {
  CELL_LOC mylow = Ui->getLocation();
  m_n.setLocation(mylow);
  m_n.applyBoundary();
  if (use_log_n) {
    l_n_n.applyBoundary();
  } else {
    n_n.applyBoundary();
  }
  if (doEvolve) {
    if (use_log_n) {
      mesh->communicate(l_n_n, m_n);
    } else {
      mesh->communicate(n_n, m_n);
    }
  }
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
  // nnsheath_yup();
  // mnsheath_yup();
  // nnsheath_ydown();
  // Works only for aiolos mesh, but this does the correct interpolation
  FieldPerp n_sheath, u_sheath;
  if (n_stag) {
    n_sheath = sliceXZ(*n_stag, mesh->yend + 1);
    u_sheath = sliceXZ(*Ui, mesh->yend + 1);
  } else {
    // fallback for bout mesh - only works for uniform meshes
    n_sheath = (-sliceXZ(*n, mesh->yend - 1) + 9 * sliceXZ(*n, mesh->yend) +
                9 * sliceXZ(*n, mesh->yend + 1) - sliceXZ(*n, mesh->yend + 2)) /
               16.;
    throw BoutException("Not implemented");
    u_sheath = nan("");
  }
  FieldPerp flux = u_sheath * n_sheath;
  S_recyc = recycle(flux);

  if (use_log_n) {
    n_n = exp(l_n_n);
  }
  calcDiffusion();
  // checkData(n_n,RGN_ALL);
  // Field3D D_stag = interp_to(D_neutrals,CELL_YLOW);
  if (use_log_n) {
    auto tmp = D_neutrals * DDY(n_n);
    mesh->communicate(tmp);
    tmp.applyBoundary("dirichlet_o4(0)");
    ddt(l_n_n) = (-Div_par(m_n, CELL_CENTRE) + gamma_rec - gamma_ion + S_recyc +
                  S_extra
                  //+ D_neutrals * Laplace(n_n)
                  //+ DDY(D_neutrals) * DDY(n_n)
                  //+ DDY(D_stag * DDY(n_n,CELL_YLOW),CELL_CENTRE)
                  + DDY(tmp)) /
                     n_n
                 //+ 100*D2DY2(n_n)
                 - loss_fraction;
  } else {
    ddt(n_n) = (
        //- Div_par(m_n, CELL_CENTRE)
        +gamma_rec - gamma_ion + S_recyc
        //+ S_extra
        //+ 100*D2DY2(n_n)
        //+ D_neutrals * Laplace(n_n)
        //- n_n * loss_fraction
    );
  }

  v_n = m_n / interp_to(n_n, mylow);
  // ddt(m_n) = 0;
  ddt(m_n) =
      (-Vpar_Grad_par(v_n, m_n, DIFF_U3) - Vpar_Grad_par(m_n, v_n, DIFF_U3) -
       Grad_par(n_n * T_n, mylow) + interp_to(gamma_CX, mylow) * (*Ui - v_n) +
       interp_to(gamma_rec, mylow) * (*Ui) + interp_to(gamma_ion, mylow) * (-v_n) +
       100 * D2DY2(m_n, CELL_DEFAULT, DIFF_C2) - interp_to(S_recyc, mylow) * v_thermal);
  // ddt(m_n) += - 10*D4DY4(m_n);
  //+((*Ui) - v_n) * interp_to(tmp, mylow);

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

Field3D ParallelNeutrals::getIonVelocitySource() const {
  ASSERT2(Ui != nullptr);
  Field3D tmp = gamma_CX + gamma_ion;
  tmp /= *n; // getCXOverN()+getIonOverN();
  return (v_n - *Ui) * (interp_to(tmp, Ui->getLocation()));
}

Field3D ParallelNeutrals::getElectronVelocitySource() const {
  ASSERT2(Ue != nullptr);
  Field3D tmp = gamma_ion / (*n);
  return (v_n - *Ue) * (interp_to(tmp, Ue->getLocation()));
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
        n_n(x, 0, z) = 1 - n_n(x, 2, z);
      }
    }
  }
}
