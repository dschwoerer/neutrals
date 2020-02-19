#include "radiation.hxx"

int main(){
  UpdatedRadiatedPower h {};
  for (double n=1e18; n < 5e24;n*=1e2){
    for (double e=1e-0;e<1e3;e*=1.02){
      printf("%.3e\t%.8e\t%.8e\t%.8e\t%.8e\n",n,e,h.ionisation(e),h.recombination(n,e),h.chargeExchange(e));
    }
    printf("\n");
    printf("\n");
  }
  return 0;
}
    
