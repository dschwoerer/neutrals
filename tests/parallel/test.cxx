#include <bout/physicsmodel.hxx>
#include <field_factory.hxx>

#include "unit.hxx"
#include "neutrals.hxx"
#include "neutrals_diffusion.hxx"

class ParallelTest : public PhysicsModel {
protected:
  // Initialisation
  int init(bool restarting) {
    unit = new SIUnit();//BohmUnit(.5, 40, 8e18, 2);
    //g.setLocation(CELL_YLOW);
    solver->add(g, "g");
    solver->add(f, "f");
    //g.setLocation(CELL_YLOW);
    neut = NeutralsFactory::create(solver, mesh, "neutrals");
    neut->setUnit(*unit);
    neut->setPlasmaDensity(f);
    neut->setPlasmaDensityStag(g);
    neut->setElectronTemperature(f);
    neut->setIonTemperature(f);
    neut->setElectronVelocity(g);
    neut->setIonVelocity(g);
    null=0;
    //src=FieldFactory::get()->create3D(".1*(y-pi)", nullptr, mesh);
    TestCrossSection * tcs = (TestCrossSection*) dynamic_cast<const TestCrossSection*>(neut->getCrossSection());
    tcs->myion=null;
    tcs->myrec=null;
    tcs->mycx=null;
    //g.setLocation(CELL_YLOW);
    
    return 0;
  }

  // Calculate time-derivatives
  int rhs(BoutReal t) {
    ddt(f) = 0;
    ddt(g) = 0;
    neut->update();
    return 0;
  }


private:
  Field3D f,g;
  Field3D src, null;
  std::unique_ptr<Neutrals> neut;
  Unit *unit;
};

// Create a default main()
BOUTMAIN(ParallelTest);
