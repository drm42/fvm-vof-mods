#ifndef _ENTHALPYFUNCTION_H_
#define _ENTHALPYFUNCTION_H_

#include <vector>
#include <iostream>
#include "Polynomial.h"

template <class T>
class EnthalpyFunction {
  typedef std::vector<T> TVector;
  typedef std::vector<TVector> VecVector;
  typedef std::vector<int> intVector;
  typedef Polynomial<T> TPoly;
  typedef std::vector<TPoly> PolyVector;

  int _numRanges; 
  PolyVector _polynomials;
  TVector _enthalpyRanges;
  TVector _tempRanges;
  T _solidTemp;
  T _liquidTemp;
  T _latentHeat;

public:
  EnthalpyFunction(): 
    _numRanges(0),
    _polynomials(),
    _enthalpyRanges(),
    _tempRanges(),
    _solidTemp(0),
    _liquidTemp(1),
    _latentHeat(0)
  {
    
  }

  void setFunctionValues(int numRanges, TVector& ranges, intVector& orders, 
		   VecVector& polys, T solidTemp, T liquidTemp, T latentHeat, 
		   T Tref)
  {
    _solidTemp = solidTemp;
    _liquidTemp = liquidTemp;
    _latentHeat = latentHeat;
    T max(0);
    T min(0);
    for(int i=0; i<numRanges; i++){
      max = ranges.at(i);
      TVector tmpcoeff(orders.at(i)+2);
      T counter(0.);
      for(int k=0; k<=orders.at(i); k++){
	tmpcoeff.at(k+1) = 1/(counter+1)*polys.at(i).at(k);
	counter += 1;
      }

      if ((solidTemp <= min) && (liquidTemp > min)){
	tmpcoeff.at(1) += latentHeat/(liquidTemp-solidTemp);
      }
      
      if ((solidTemp > min) && (solidTemp < max)){
	TPoly tmpPoly(orders.at(i)+1,tmpcoeff);
	if (_numRanges != 0){
	  tmpPoly.setCoefficient(0,_enthalpyRanges.at(_numRanges-1)-
				 tmpPoly.eval(min));
	}
	_enthalpyRanges.push_back(tmpPoly.eval(solidTemp));
	_tempRanges.push_back(solidTemp);
	_polynomials.push_back(tmpPoly);
	_numRanges++;
	min = solidTemp;
	tmpcoeff.at(1) += latentHeat/(liquidTemp-solidTemp);
      }

      if ((liquidTemp > min) && (liquidTemp < max)){
	TPoly tmpPoly(orders.at(i)+1,tmpcoeff);
	tmpPoly.setCoefficient(0,_enthalpyRanges.at(_numRanges-1)-
			       tmpPoly.eval(min));
	_enthalpyRanges.push_back(tmpPoly.eval(liquidTemp));
	_tempRanges.push_back(liquidTemp);
	_polynomials.push_back(tmpPoly);
	_numRanges++;
	min = liquidTemp;
	tmpcoeff.at(1) -= latentHeat/(liquidTemp-solidTemp);
      }
      
      TPoly tmpPoly(orders.at(i)+1,tmpcoeff);
      if (_numRanges != 0){
	tmpPoly.setCoefficient(0,_enthalpyRanges.at(_numRanges-1)-
			       tmpPoly.eval(min));
      }
      _enthalpyRanges.push_back(tmpPoly.eval(max));
      _tempRanges.push_back(max);
      _polynomials.push_back(tmpPoly);
      _numRanges++;
      min = max;
    }

    T refEnthalpy(findEnthalpy(Tref));
    for (int i=0; i<_numRanges; i++){
      _polynomials.at(i).incCoefficient(0,-refEnthalpy);
      _enthalpyRanges.at(i) -= refEnthalpy;
    }
  }

  T findEnthalpy(T temp){
    for (int i=0; i<_numRanges; i++){
      if(temp <= _tempRanges.at(i)){
	return _polynomials.at(i).eval(temp);
      }
    }
    return -1;
  }

  T findEnthalpyInverse(T enthalpy){
    for (int i=0; i<_numRanges; i++){
      if(enthalpy <= _enthalpyRanges.at(i)){
	T min(0.);
	if(i!=0) {
	  min = _tempRanges.at(i-1);
	}
	T max(_tempRanges.at(i));
	//return _polynomials.at(i).solve(min,max,enthalpy,
	//				(_liquidTemp-_solidTemp)/100,100);
	return _polynomials.at(i).solve(min,max,enthalpy,1.0e-10,100);
      }
    }
    return -1;
  }

  T findDHdT(T temp){
    for (int i=0; i<_numRanges; i++){
      if(temp <= _tempRanges.at(i)){
	return _polynomials.at(i).derivativeEval(temp);
      }
    }
    return -1;
  }

  void printEnthalpyFunction(){
    std::cout << "Number of Ranges :" << _numRanges << std::endl;
    std::cout << "Solid Temp: " << _solidTemp << " Liquid Temp: " 
	      << _liquidTemp << " Latent Heat: " << _latentHeat << std::endl;
    for(int i=0; i<_numRanges; i++){
      std::cout << "Range: " << i << " Max Temp/Enthalpy ";
      std::cout << _tempRanges.at(i) << "/" << _enthalpyRanges.at(i); 
      std::cout << std::endl;
      _polynomials.at(i).print();
    }
  }

};

#endif
