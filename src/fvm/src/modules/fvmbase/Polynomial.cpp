#include "Array.h"

template <class T>
class Polynomial {

  typedef Array<T> TArray;

  shared_ptr<TArray> _coefficients;
  T _order;

 public:

 Polynomial(int poly_order, TArray coeff) :
  _order(poly_order),
  _coefficients(new TArray(poly_order + 1))
    {
      for (int i=0; i<=poly_order; i++){
	_coefficients[i] = coeff[i];
      }
    }

  int getOrder(){
    return _order;
  }

  T eval(T x){
    T value(0.);
    for(int i=0; i<=_order; i++){
      T power(1.);
      for(int k=0; k<i; k++){
	power = power*x;
      } 
      value = value + power*_coefficients[i];
    }
    return value;
  }
};

int main() {
  typedef Array<double> dArray;
  typedef Polynomial<double> dPolynomial;

  double ans;
  double order;

  shared_ptr<dArray> coef(new dArray(2));
  coef[0] = 1;
  coef[1] = 1;
  dPolynomial poly(1,coef);
  coef[1] = 2;

  ans = poly.eval(1);
  order = poly.getOrder();

  return 0;
  
}
