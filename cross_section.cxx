#include "cross_section.hxx"
#include <boutexception.hxx>
#include <globals.hxx>
// Here should be the crosssections();

CrossSection::CrossSection(RadiatedPower *atom) : atom(atom) { ; }

Field3D CrossSection::ionisation_rate(const Field3D &T_e) {
  Field3D result;
  result.allocate();
#pragma omp parallel
  for (auto i = T_e.begin(); !i.done(); ++i) {
    const BoutReal Te = T_e[i];
    result[i] = atom->ionisation(Te);
  }
  return result;
}

// BoutReal CrossSection::ionisation_rate(const BoutReal T_e) {
//   return atom.ionisation(T_e);
// }

Field3D CrossSection::recombination_rate(const Field3D &n, const Field3D &T_e) {
  Field3D result;
  result.allocate();
  for (auto i : result) {
    result[i] = atom->recombination(n[i], T_e[i]);
  }
  return result;
}

Field3D CrossSection::power(const Field3D &T, const Field3D &n, const Field3D &imp) {
  Field3D result;
  result.allocate();
  for (auto i : result) {
    result[i] = atom->power(T[i], n[i], imp[i]);
  }
  return result;
}

// // TODO: only compute for actual domain
// Field3D CrossSection::recombination_rate(const Field3D &n, const BoutReal T_e) {
//   Field3D result;
//   // result.allocate();
//   for (auto i: result){
//     result[i]= atom.recombination(n[i],T_e);
//   }
//   //#endif*/
//   return result;
// }
/*
Field3D CrossSection::recombination_rate(const BoutReal n, const Field3D & T_e){
  Field3D result;
  result.allocate();

  for(int i=0;i<mesh->LocalNx;i++)
    for(int j=0;j<mesh->LocalNy;j++)
      for(int k=0;k<mesh->LocalNz;k++){
        BoutReal Te=T_e(i,j,k);
        result(i,j,k) = atom.recombination(n,Te);
      }

  return result;
}

BoutReal CrossSection::recombination_rate(const BoutReal n, const BoutReal T_e){
  return  atom.recombination(n,T_e);
}
*/
Field3D CrossSection::cx_rate(const Field3D &T_e) {
  Field3D result;
  result.allocate();
  for (auto i : T_e) {
    result[i] = atom->chargeExchange(T_e[i]);
  }
  //#endif
  return result;
}

// BoutReal CrossSection::cx_rate(const BoutReal T_e) {
//   // std::cerr << atom.chargeExchange(T_e)<<"\t";
//   return atom.chargeExchange(T_e);
// }
