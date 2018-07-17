#pragma once

#include "radiation.hxx"
#include <options.hxx>

class RadiatedPowerFactory {
public:
  static RadiatedPower *create(Options *options);
  static RadiatedPower *create(const std::string &name, Options *options = nullptr);
};
