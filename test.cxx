#include <bout/physicsmodel.hxx>

#include "neutrals.hxx"
#include "neutrals_diffusion.hxx"

class MonitorExample : public PhysicsModel {
protected:
  // Initialisation
  int init(bool restarting) {
    solver->add(f, "f");
    neut = NeutralsFactory::create(solver,Options::getRoot()->getSection("neutrals") );
    neut->setPlasmaDensity(f);
    neut->setElectronTemperature(f);
    neut->setIonTemperature(f);
    neut->setElectronVelocity(f);
    neut->setIonVelocity(f);
    return 0;
  }
  
  // Calculate time-derivatives
  int rhs(BoutReal t) {
    ddt(f) = -f;
    neut->update();
    return 0;
  } 
  
  // This called every output timestep
  int outputMonitor(BoutReal simtime, int iter, int NOUT) {
    output.write("\nOutput monitor, time = %e, step %d of %d\n", 
                 simtime, iter, NOUT);
    return 0;
  }
  
  // This called every timestep
  int timestepMonitor(BoutReal simtime, BoutReal dt) {
    output.write("\nTimestep monitor, time = %e, dt = %e\n", simtime, dt);
    return 0;
  }
  
private:
  Field2D f;
  std::unique_ptr<Neutrals> neut;
};

// Create a default main()
BOUTMAIN(MonitorExample);
