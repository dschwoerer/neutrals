#include "neutrals_diffusion.hxx"

#include "bout/surfaceiter.hxx"
#include "bout/constants.hxx"
/*

  OPTION(options, neutrals,              false) ;
  OPTION(options, neutrals_nabla_gamma,  false) ;
  OPTION(options, neutrals_diffusion  ,  false) ;
  OPTION(options, neutrals_diffusion_spacial, false) ;
  OPTION(options, mu_neutrals         ,  2e-5 ) ;   // m^2/s                                                                                                            
  OPTION(options, neutrals_falloff    ,  4.0  ) ;   // m                                                                                                                
  OPTION(options, neutrals_limit      ,  1e-8 ) ;   // relative                                                                                                         
  OPTION(options, recycling_coeff     ,  0.1  ) ;
  OPTION(options, neutrals_equi       ,  false) ;
  if (!neutrals_equi){
    OPTION(options, neutrals_static     ,  false) ;
  } else {
    output.write("\tneutrals are static, as equilibrium rates are used\n");
    neutrals_static=true;
  }
  OPTION(options, neutrals_onlyion    ,  false) ;
  OPTION(options, vort_neut_correct   ,  true ) ;


 */

DiffusionNeutrals::DiffusionNeutrals(Solver * solver_, Options * options):
  //Neutrals(solver_,options),
  solver(solver_){
  options->get("evolve",doEvolve,true);
  OPTION(options, equi_rates       ,  false) ;
  if (equi_rates && doEvolve){
    throw BoutException("DiffusionNeutrals:: cannot have equilibrium rates with evolving neutrals!");
  }
  if (doEvolve){
    solver->add(n_n,"neutrals");
  }
  n_n=1.;
}

void
DiffusionNeutrals::update(){
  if (!equi_rates){
    this->calcRates();
  }
  if (doEvolve){
    this->evolve();
  }
}


void DiffusionNeutrals::evolve(){
  if (n_stag == nullptr){
    if (n->getLocation() != Ui->getLocation()){
      throw BoutException("DiffusionNeutrals:: density and velocity at different location, but staggered density not given!");
    } else {
      n_stag=n;
    }
  }
  FieldPerp flux=sliceXZ(*Ui,mesh->yend+1)*sliceXZ(*n_stag,mesh->yend+1);
  Field3D S_recyc = recycle(flux);
  //debug2=S_recyc;
  ddt(n_n) = gamma_rec - gamma_ion
    + S_recyc;
  //debug3=gamma_rec - gamma_ion;
  ddt(n_n) += - n_n * n_n_sink;
}

void DiffusionNeutrals::setDensityStag(const Field3D & n_stag_){
  n_stag=&n_stag_;
}

void DiffusionNeutrals::calcRates(){
  #warning not added
#if 0
  
  Field3D T_in_eV=T*T_0;
  Field3D nnn0oOmega = n_n*(n_0/Omega_i);
  gamma_ion_over_n   = nnn0oOmega*hydrogen.ionisation_rate(T_in_eV);
  gamma_ion          = gamma_ion_over_n*n;
  if (!neutrals_onlyion){
    Field3D energy=SQ(interp_to(U,CELL_CENTER))*(m_i/2*c_s*c_s/e); // in eV
    energy+=T_in_eV;
    //limit_at_most(energy,1e4); //10 keV
    gamma_CX_over_n    = hydrogen.cx_rate(energy);
    gamma_CX_over_n   *= nnn0oOmega;
    gamma_CX           = gamma_CX_over_n*n;
    gamma_rec_over_n   = n * hydrogen.recombination_rate(n*n_0,T_in_eV)
      *(n_0/Omega_i); // guad cells not set
    gamma_rec          = n * gamma_rec_over_n;
  } else {
    gamma_CX         = 0;
    gamma_CX_over_n  = 0;
    gamma_rec        = 0;
    gamma_rec_over_n = 0;
  }
#endif
}


Field3D
DiffusionNeutrals::recycle(const FieldPerp &flux){
  Field3D result;
  result.allocate();
  //for (int i=0;i<mesh->LocalNx*mesh->LocalNy*mesh->LocalNz;++i){
  for (auto i:result){
    result[i]=nan("");
  }
  //Set a recycling source in first and last processors
  //BoutReal L;	//Parallel length coordinate
  //BoutReal TargetFluxInner,TargetFluxOuter;
  SurfaceIter s(mesh);
  MPI_Comm comm;
  int nproc;
  for(s.first(); !s.isDone(); s.next()){
    int i = s.xpos;
    if(s.closed()){
      // todo: change loops
      for(int j=0; j<mesh->LocalNy; j++){
	for(int k=0; k <mesh->LocalNz; k++) {
	  result(i,j,k)=0;
	}
      }
      throw BoutException("What???\n");
      continue;	//Skip if closed flux surface
    }
    Coordinates *coord = mesh->coordinates();
    //COORD(coord);
    if (recycling_falloff){
      comm = s.communicator();
      MPI_Comm_size(comm,&nproc);
      BoutReal TargetFluxOuter[mesh->LocalNz];
      for(int k=0; k <mesh->LocalNz; k++) {
	//if(s.lastY())
	//TargetFluxOuter = Flux(i,mesh->yend+1,k)/mesh->Bxy(i,mesh->yend+1);
	// todo: replace abs with min(,0)?
	TargetFluxOuter[k] = abs(flux(i,k)*recycling_fraction);
        //if (TargetFluxOuter[k] == 0){throw BoutException("fuuu");}
      }
      MPI_Bcast(TargetFluxOuter,mesh->LocalNz,MPI_DOUBLE,nproc-1,comm);
      //if(s.firstY())		//On inner target boundary, broadcast target flux
      //TargetFluxInner = Flux(i,2,k)/mesh->Bxy(i,2);
      //MPI_Bcast(&TargetFluxInner,1,MPI_DOUBLE,0,comm);
      //BoutReal sum=0;
      //BoutReal den;
      BoutReal tmpsum=0;
      for(int j=2; j<mesh->GlobalNy-2; j++){
        BoutReal L = (mesh->GlobalNy-4.5-mesh->YGLOBAL(j))*coord->dy(i,j)*sqrt(coord->g_22(i,j)); //Make zero point between guard cell and first cell
        //printf("%d %g  ",mesh->YGLOBAL(j),L);
        tmpsum+=(2/(recycling_falloff*sqrt(TWOPI)))*exp(-SQ(L)/(2*SQ(recycling_falloff)))*coord->dy(i,j);
      }
      //printf("%d %g\n",mesh->GlobalNy,tmpsum);
      for(int j=0; j<mesh->LocalNy; j++){
	// TODO: BUG: L_para should not depend on j?
	//BoutReal L_para=(mesh->GlobalNy-4)*coord->dy(i,j)*sqrt(coord->g_22(i,j));
	BoutReal L = (mesh->GlobalNy-4.5-mesh->YGLOBAL(j))*coord->dy(i,j)*sqrt(coord->g_22(i,j)); //Make zero point between guard cell and first cell
	// TODO: replace abs(?) with min(?,0)?
	//BoutReal tmp=1/falloff*exp(-L/falloff)/(1-exp(-L_para/falloff));
        BoutReal tmp=(2/(recycling_falloff*sqrt(TWOPI)))*exp(-SQ(L)/(2*SQ(recycling_falloff)));
        tmp/=tmpsum;
        //if(tmp == 0){throw BoutException("maeh");}
	for (int k=0;k<mesh->LocalNz;k++){
          //if(TargetFluxOuter[k] == 0){throw BoutException("maeh");}
	  result(i,j,k) = tmp*TargetFluxOuter[k];
	}
      }
    } else {
      for(int j=0; j<mesh->LocalNy; j++){
	for(int k=0; k <mesh->LocalNz; k++) {
	  result(i,j,k)=0;
	}
      }
      if(s.lastY()){
        for(int k=0; k <mesh->LocalNz; k++) {
          result(i,mesh->yend,k)=flux(i,k)/coord->dy(i,mesh->yend);
        }
      }
    }
  }
  //limit_at_least(result,1e-30,true);
  return result;
}
