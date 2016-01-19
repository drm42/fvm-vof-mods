// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.


#ifndef _TIMEDERIVATIVEENTHALPYDISCRETIZATION_H_
#define _TIMEDERIVATIVEENTHALPYDISCRETIZATION_H_

#include "CRMatrix.h"
#include "Field.h"
#include "MultiField.h"
#include "MultiFieldMatrix.h"
#include "Mesh.h"
#include "Discretization.h"
#include "StorageSite.h"
#include <iostream>

template<class X, class Diag, class OffDiag>
class TimeDerivativeEnthalpyDiscretization : public Discretization
{
public:
  typedef typename NumTypeTraits<X>::T_Scalar T_Scalar;
  
  typedef CRMatrix<Diag,OffDiag,X> CCMatrix;
  typedef typename CCMatrix::DiagArray DiagArray;
  typedef typename CCMatrix::PairWiseAssembler CCAssembler;

  typedef Array<X> XArray;
  typedef Array<T_Scalar> TArray;
  typedef Array<int> IntArray;

  TimeDerivativeEnthalpyDiscretization(const MeshList& meshes,
                               const GeomFields& geomFields,
                               Field& varField,
                               Field& varN1Field,
                               Field& varN2Field,
                               const Field& rhoCpField,
			       const Field& densityField,
			       const Field& enthalpyField,
			       const Field& enthalpyN1Field,
			       const Field& enthalpyInverseField,
			       const Field& dHdTField,
                               const T_Scalar dT,
			       const T_Scalar latentHeat,
			       const T_Scalar solidTemp,
			       const T_Scalar liquidTemp,
			       const T_Scalar Tref) :
      Discretization(meshes),
      _geomFields(geomFields),
      _varField(varField),
      _varN1Field(varN1Field),
      _varN2Field(varN2Field),
      _rhoCpField(rhoCpField),
      _densityField(densityField),
      _enthalpyField(enthalpyField),
      _enthalpyN1Field(enthalpyN1Field),
      _enthalpyInverseField(enthalpyInverseField),
      _dHdTField(dHdTField),
      _dT(dT),
      _latentHeat(latentHeat),
      _solidTemp(solidTemp),
      _liquidTemp(liquidTemp),
      _Tref(Tref)
      {}
  
  void discretize(const Mesh& mesh, MultiFieldMatrix& mfmatrix,
                  MultiField& xField, MultiField& rField)
  {

    const StorageSite& cells = mesh.getCells();

    const TArray& rhoCp =
      dynamic_cast<const TArray&>(_rhoCpField[cells]);

    const TArray& density =
      dynamic_cast<const TArray&>(_densityField[cells]);

    const TArray& enthalpy =
      dynamic_cast<const TArray&>(_enthalpyField[cells]);

    const TArray& h_old =
      dynamic_cast<const TArray&>(_enthalpyN1Field[cells]);

    const TArray& enthalpyInverse =
      dynamic_cast<const TArray&>(_enthalpyInverseField[cells]);

    const TArray& dHdT =
      dynamic_cast<const TArray&>(_dHdTField[cells]);

    const TArray& cellVolume =
      dynamic_cast<const TArray&>(_geomFields.volume[cells]);
    
    const MultiField::ArrayIndex cVarIndex(&_varField,&cells);
    CCMatrix& matrix =
      dynamic_cast<CCMatrix&>(mfmatrix.getMatrix(cVarIndex,cVarIndex));

    DiagArray& diag = matrix.getDiag();

    const XArray& x = dynamic_cast<const XArray&>(_varField[cells]);
    const XArray& xN1 = dynamic_cast<const XArray&>(_varN1Field[cells]);
    
    XArray& rCell = dynamic_cast<XArray&>(rField[cVarIndex]);

    
    const int nCells = cells.getSelfCount();

   
    if (_varN2Field.hasArray(cells))
    {
        // second order
      throw CException("Enthalpy model not implemented for 2nd order or higher time discretization");
        /*const XArray& xN2 = dynamic_cast<const XArray&>(_varN2Field[cells]);

        T_Scalar onePointFive(1.5);
        T_Scalar two(2.0);
        T_Scalar pointFive(0.5);
        if (_geomFields.volumeN1.hasArray(cells))
	{
	    const TArray& cellVolumeN1 = 
	      dynamic_cast<const TArray&>(_geomFields.volumeN1[cells]);
	    const TArray& cellVolumeN2 = 
              dynamic_cast<const TArray&>(_geomFields.volumeN2[cells]);
            for(int c=0; c<nCells; c++)
	    {
                const T_Scalar rhoVbydT = rhoCp[c]*cellVolume[c]/_dT;
                const T_Scalar rhobydT = rhoCp[c]/_dT;
		const T_Scalar term1 = onePointFive*cellVolume[c];
                const T_Scalar term2 = two*cellVolumeN1[c];
                const T_Scalar term3 = pointFive*cellVolumeN2[c];
                rCell[c] -= rhobydT*(term1*x[c]- term2*xN1[c]
                                      + term3*xN2[c]);
                diag[c] -= rhoVbydT*onePointFive;
	    }
	}
	else
	{
            for(int c=0; c<nCells; c++)
            {
                const T_Scalar rhoVbydT = rhoCp[c]*cellVolume[c]/_dT;
                rCell[c] -= rhoVbydT*(onePointFive*x[c]- two*xN1[c]
                                      + pointFive*xN2[c]);
                diag[c] -= rhoVbydT*onePointFive;
	    }
	    }*/
    }
    else
    {
        if (_geomFields.volumeN1.hasArray(cells))
	{
	    throw CException("Enthalpy model not implemented for dynamic meshes");
	    /*const TArray& cellVolumeN1 =
	      dynamic_cast<const TArray&>(_geomFields.volumeN1[cells]);
	    for(int c=0; c<nCells; c++)
            {	    
	        const T_Scalar rhoVbydT = rhoCp[c]*cellVolume[c]/_dT;
		const T_Scalar rhobydT = rhoCp[c]/_dT;
	        rCell[c] -= rhobydT*(cellVolume[c]*x[c] 
                                     - cellVolumeN1[c]*xN1[c]);
		
	        diag[c] -= rhoVbydT;
		
		}*/
	}
        else if (_geomFields.ibTypeN1.hasArray(cells))
	{
	  throw CException("Enthalpy model not implemented for interface boundaries");
	  /*const IntArray& ibTypeN1 =
	      dynamic_cast<const IntArray&>(_geomFields.ibTypeN1[cells]);
	    const IntArray& ibType =
	      dynamic_cast<const IntArray&>(_geomFields.ibType[cells]);
	    for(int c=0; c<nCells; c++)
            {	    
	        const T_Scalar rhoVbydT = rhoCp[c]*cellVolume[c]/_dT;

                if (ibTypeN1[c] == Mesh::IBTYPE_FLUID &&
                    (ibType[c] == Mesh::IBTYPE_FLUID))
                  rCell[c] -= rhoVbydT*(x[c]  - xN1[c]);
                else if (ibType[c] == Mesh::IBTYPE_FLUID)
                  rCell[c] -= rhoVbydT*(x[c] );
	        diag[c] -= rhoVbydT;
		}*/
	}
        else
	{
	  for(int c=0; c<nCells; c++)
            {
	      /*if (c == 0) {
		std::cout << h_old[c] << ' ' << enthalpy[c] << std::endl;
		std::cout << dHdT[c] << ' ' << enthalpyInverse[c] << std::endl;
		std::cout << x[c] << std::endl;
		}*/
	      const T_Scalar VbydT = cellVolume[c]/_dT;
	      rCell[c] += VbydT*(h_old[c]-enthalpy[c])+VbydT*dHdT[c]*
		enthalpyInverse[c]-VbydT*dHdT[c]*x[c];
	      diag[c] -= VbydT*dHdT[c];
	
	    }
	}
    }
  
  }

  

private:
  const GeomFields& _geomFields;
  const Field& _varField;
  const Field& _varN1Field;
  const Field& _varN2Field;
  const Field& _rhoCpField;
  const Field& _densityField; 
  const Field& _enthalpyField;
  const Field& _enthalpyN1Field;
  const Field& _enthalpyInverseField;
  const Field& _dHdTField;
  const T_Scalar _dT;
  const T_Scalar _latentHeat;
  const T_Scalar _solidTemp;
  const T_Scalar _liquidTemp;
  const T_Scalar _Tref;
};

#endif
