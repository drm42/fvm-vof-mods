// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _THERMALFIELDS_H_
#define _THERMALFIELDS_H_


#include "Field.h"
#include "FieldLabel.h"

struct ThermalFields
{
  ThermalFields(const string baseName);

  Field temperature;
  Field temperatureN1;
  Field temperatureN2;
  Field specificHeat;
  Field density;
  Field enthalpy;
  Field enthalpyN1;
  Field enthalpyInverse;
  Field dHdT;
  Field meltFrac;
  Field continuityResidual;
  Field heatFlux;
  Field temperatureGradient;
  Field conductivity;
  Field source;
  Field convectionFlux;

  Field zero;                     //used to fill in continuityResidual
  Field one;                      //used to fill in density
};

#endif
