#ifndef _POLYNOMIAL_H_
#define _POLYNOMIAL_H_

#include <vector>
#include <iostream>
#include <cmath>

template <class T>
class Polynomial {

  typedef std::vector<T> TVector;

  TVector _coefficients;
  T _order;

 public:

 Polynomial(int poly_order, TVector coeff) :
  _order(poly_order),
  _coefficients(poly_order+1)
  {
    _coefficients = coeff;
  }

  int getOrder(){
    return _order;
  }

  T eval(T x){
    T value(0);
    for(int i=0; i<=_order; i++){
      value = value+std::pow(x,i)*_coefficients.at(i);
    }
    return value;
  }

  T derivativeEval(T x){
    T value(0);
    for(int i=1; i<=_order; i++){
      value = value+std::pow(x,i-1)*i*_coefficients.at(i);
    }
    return value;
  }

  T findZero(T lowerBound, T uppderBound, T tol, int maxIter){
    T a(lowerBound);
    T b(uppderBound);
    T fa(eval(a));
    T fb(eval(b));
    T swap(0);
    T s(0);
    T fs(0);
    T d(0);
    
    if (fa*fb >= 0){
      if (fa < fb){
	return a;
      }
      else {
	return b;
      }
    } 

    if (std::abs(fa) < std::abs(fb)){
      swap = a;
      a = b;
      b = swap;
      swap = fa;
      fa = fb;
      fb = swap;
    }
    T c(a);
    T fc(fa);
    bool mflag(true);
    int i(0);

    while ( (fb != 0) && (std::abs(b-a) > tol)){
      if ((fa != fc) && (fb != fc)) {
	s = a*fb*fc/(fa-fb)/(fa-fc)+b*fa*fc/(fb-fa)/(fb-fc)+
	  c*fa*fb/(fc-fa)/(fc-fb);
      }
	
      else {
	s = b-fb*(b-a)/(fb-fa);
      }

      if ((!(((s>(3*a+b)/4) && (s<b))||((s<(3*a+b)/4)&&(s>b)))) ||
	  ((mflag == true) && (std::abs(s-b) >= std::abs(b-c)/2)) ||
	  ((mflag == false) && (std::abs(s-b) >= std::abs(c-d)/2)) ||
	  ((mflag == true) && (std::abs(b-c) < tol)) ||
	  ((mflag == false) && (std::abs(c-d) < tol))){
	s = (a+b)/2;
	mflag = true;
      }
      else {
	mflag = false;
      }

      fs = eval(s);
      d = c;
      c = b;
      fc = fb;
      
      if (fa*fs < 0) {
	b = s;
	fb = fs;
      }
      else {
	a = s;
	fa = fs;
      }

      if(std::abs(fa) < std::abs(fb)){
	swap = a;
	a = b;
	b = swap;
	swap = fa;
	fa = fb;
	fb = swap;
      }
      i++;
      if(i > maxIter){
	return b;
      }
    } 

    return b;

  }

  T solve(T lowerBound,T upperBound, T rightHandSide, T tol,
	  int maxIters){
    TVector new_coeff(_order);
    new_coeff = _coefficients;
    new_coeff.at(0) -= rightHandSide;
    Polynomial<T> eqn(_order, new_coeff);
    return eqn.findZero(lowerBound, upperBound, tol, maxIters);
  }

  int setCoefficient(int i, T value){
    if ((i > _order) || (i < 0)){
      return 1;
    }
    else{
      _coefficients.at(i) = value;
      return 0;
    }
  }

  int incCoefficient(int i, T value){
    if ((i > _order) || (i < 0)){
      return 1;
    }
    else{
      _coefficients.at(i) += value;
      return 0;
    }
  }

  void print(){
    std::cout << "Order: " << _order << " Coefficients: ";
    for (int i=0; i<=_order; i++){
      std::cout << _coefficients.at(i) << " ";
    }
    std::cout << std::endl;
  }

};

#endif
