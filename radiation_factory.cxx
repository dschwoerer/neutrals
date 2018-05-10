#include "radiation_factory.hxx"
#include "radiation.hxx"


RadiatedPower * RadiatedPowerFactory::create(Options * options){
  std::string name;
  options->get("RadiationType", name , "Not Set");

  if (name == "hydrogenradiatedpower") {
    return new HydrogenRadiatedPower();
  } else if (name == "updatedradiatedpower") {
    return new UpdatedRadiatedPower();
  } else if (name == "testingpower") {
    return new TestingPower(options->getSection(name));
  } else {
    throw BoutException("Cannot handle requested CrossSectionType.\n"
                        "Requested %s - but we only support:\n"
                        "	* HydrogenRadiatedPower\n"
                        "	* UpdatedRadiatedPower\n"
                        "	* TestingPower\n",name.c_str());
  }
}
