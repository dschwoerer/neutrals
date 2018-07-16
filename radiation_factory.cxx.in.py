#!/usr/bin/env python3


models=[
    "RadiatedPower",
    # Needs an argument - file to read
    #"InterpRadiatedPower",
    "HutchinsonCarbonRadiation",
    "HydrogenRadiatedPower",
    "UpdatedRadiatedPower",
    ## commented out:
    # "PostJensen"
    # Testing only:
    "TestingPower"
]
takesOptions={}
for model in models:
    takesOptions[model]=False
takesOptions["TestingPower"]=True
              

def printl(*args):
    print(*args,end='')

printl("""#include "radiation_factory.hxx"
#include "radiation.hxx"


RadiatedPower * RadiatedPowerFactory::create(Options * options){
  std::string name;
  options->get("RadiationType", name , "updatedradiatedpower");
  return create(name, options);
}
RadiatedPower * RadiatedPowerFactory::create(const std::string & name, Options * options){

  """)
for model in models:
    printl("""if (name == "%s") {
    return new %s(%s);
  } else """%(model.lower(), model,"options->getSection(name)" if takesOptions[model] else ""))
printl("""{
    throw BoutException("Cannot handle requested CrossSectionType.\\n"
                        "Requested %s - but we only support:\\n"\
""")
for model in models:
    printl('\n                        "\t* '+model+'\\n"')
print(""",name.c_str());
  }
}""")
