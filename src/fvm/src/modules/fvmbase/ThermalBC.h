// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "FloatVarDict.h"
#include "AMG.h"

template<class T>
struct ThermalBC : public FloatVarDict<T>
{
  ThermalBC()
  {
      this->defineVar("specifiedTemperature",T(300.0));
      this->defineVar("specifiedHeatFlux",T(0.0));
      this->defineVar("convectiveCoefficient", T(0.0));
      this->defineVar("farFieldTemperature", T(300.0));
      this->defineVar("surfaceEmissivity",T(1.0));
  }
  string bcType;
};

template<class T>
struct ThermalVC : public FloatVarDict<T>
{
  ThermalVC()
  {
      this->defineVar("thermalConductivity",T(1.0));
      this->defineVar("density", T(1.0));
      this->defineVar("specificHeat", T(1.0));
      	
  }
  string vcType;
};

template<class T>
struct ThermalModelOptions : public FloatVarDict<T>
{
  ThermalModelOptions()
  {
    this->defineVar("initialTemperature",T(300.0));
    this->defineVar("timeStep", T(1e-7));
    this->defineVar("latentHeat", T(1e5));
    this->defineVar("solidTemp", T(400));
    this->defineVar("liquidTemp", T(410));
    this->relativeTolerance=1e-8;
    this->absoluteTolerance=1e-16;
    this->linearSolver = 0;
    this->useCentralDifference=false;
    this->transient=false;
    this->enthalpyModel=false;
    this->timeDiscretizationOrder = 1;
    this->polynomialCp=false;
  }
  double relativeTolerance;
  double absoluteTolerance;
  bool useCentralDifference;
  LinearSolver *linearSolver;
  bool transient;
  bool enthalpyModel;
  int timeDiscretizationOrder;
  bool polynomialCp;
  string cpFile;
#ifndef SWIG
  LinearSolver& getLinearSolver()
  {
    if (this->linearSolver == 0)
    {
        LinearSolver* ls(new  AMG());
        ls->relativeTolerance = 1e-1;
        ls->nMaxIterations = 20;
        ls->verbosity=0;
        this->linearSolver = ls;
    }
    return *this->linearSolver ;
  }
#endif
};

