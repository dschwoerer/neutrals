#include "interpolation.hxx"

#include "helper.hxx"
#include "neutrals_parallel.hxx"
#include "derivs.hxx"
#include "unit.hxx"
#include "bout/solver.hxx"

ParallelNeutrals::ParallelNeutrals(Solver *solver, Mesh *mesh, CrossSection * cs, Options *options)
  : DiffusionNeutrals(solver, mesh, cs, options) {

  OPTION(options, momentum_name, "m_n");
  if (doEvolve) {
    //m_n.setLocation(CELL_YLOW);
    solver->add(m_n, momentum_name.c_str());
  }
}

void ParallelNeutrals::setNeutrals(const Field3D &n_n_, const Field3D &m_n_) {
  n_n = n_n_;
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
      mesh->communicate(l_n_n,m_n);
    } else {
      mesh->communicate(n_n,m_n);
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
  //ASSERT1(mylow==CELL_YLOW);
  //nnsheath_yup();
  //mnsheath_yup();
  //nnsheath_ydown();
  //[-0.0625  0.5625  0.5625 -0.0625]
  // [-1.,  9.,  9., -1.] / 16
  //[ 0.3125  0.9375 -0.3125  0.0625]
  //n_sheath=sliceXZ(*n_stag, mesh->yend + 1);
  auto n_sheath=(-sliceXZ(*n, mesh->yend -1)+9*sliceXZ(*n, mesh->yend )+9*sliceXZ(*n, mesh->yend + 1)-sliceXZ(*n, mesh->yend + 2))/16.;
  FieldPerp flux = sliceXZ(*Ui, mesh->yend + 1) * n_sheath;
  //printf("\t%g\n",n_sheath(mesh->xstart,0));//(*n_stag)(mesh->xstart,mesh->yend+1,0));
  S_recyc = recycle(flux);
  calcDiffusion();
  if (use_log_n) {
    n_n=exp(l_n_n);
  }
  // checkData(n_n,RGN_ALL);
  ddt(n_n) = (
              - Div_par(m_n, CELL_CENTRE)
              + gamma_rec - gamma_ion + S_recyc
              //+ 100*D2DY2(n_n)
              + D_neutrals * Laplace(n_n)
              -n_n * loss_fraction
              + S_extra
              );
  if (use_log_n) {
    ddt(l_n_n) = ddt(n_n)/n_n;
  }

  v_n = m_n / interp_to(n_n, mylow);
  auto tmp = gamma_CX + gamma_ion;
  //ddt(m_n) = 0;
  ddt(m_n) = (
              - Vpar_Grad_par(v_n, m_n)
           // - Vpar_Grad_par(m_n, v_n);
              - Grad_par(n_n * T_n, mylow)
              + interp_to(gamma_CX,  mylow) * (*Ui - v_n)
              + interp_to(gamma_rec, mylow) * (*Ui)
              + interp_to(gamma_ion, mylow) * ( - v_n)
              + 1000 * D2DY2(m_n,CELL_DEFAULT, DIFF_C2)
              - interp_to(S_recyc, mylow) * v_thermal
             );
  //ddt(m_n) += - 10*D4DY4(m_n);
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
        n_n(x, 0, z) = 1-n_n(x,  2, z);
      }
    }
  }
}
