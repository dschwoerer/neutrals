#include "radiation_factory.hxx"
#include "radiation.hxx"


RadiatedPower * RadiatedPowerFactory::create(Options * options){
  std::string name;
  options->get("RadiationType", name , "updatedradiatedpower");
  return create(name, options);
}
RadiatedPower * RadiatedPowerFactory::create(const std::string & name, Options * options){

  if (name == "radiatedpower") {
    return new RadiatedPower();
  } else if (name == "hutchinsoncarbonradiation") {
    return new HutchinsonCarbonRadiation();
  } else if (name == "hydrogenradiatedpower") {
    return new HydrogenRadiatedPower();
  } else if (name == "updatedradiatedpower") {
    return new UpdatedRadiatedPower();
  } else if (name == "testingpower") {
    return new TestingPower(options->getSection(name));
  } else {
    throw BoutException("Cannot handle requested CrossSectionType.\n"
                        "Requested %s - but we only support:\n"
                        "	* RadiatedPower\n"
                        "	* HutchinsonCarbonRadiation\n"
                        "	* HydrogenRadiatedPower\n"
                        "	* UpdatedRadiatedPower\n"
                        "	* TestingPower\n",name.c_str());
  }
}
