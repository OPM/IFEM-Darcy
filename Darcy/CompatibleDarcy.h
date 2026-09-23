// $Id$
//==============================================================================
//!
//! \file CompatibleDarcy.h
//!
//! \date Sep 21 2026
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Integrand for Darcy flow with a div-compatible basis.
//!
//==============================================================================

#ifndef _COMPATIBLE_DARCY_H_
#define _COMPATIBLE_DARCY_H_

#include "Darcy.h"


/*!
  \brief Darcy flow integrand using a divergence-compatible flux basis.
  \details The Darcy flux has one scalar basis for each spatial component,
  while pressure uses the final scalar basis.
*/

class CompatibleDarcy : public Darcy
{
  //! \brief Element right-hand-side vector indices.
  enum ResidualVectors { Fqx = 1, Fqy = 2, Fqz = 3, Fp = 4, NVEC = 5 };

  //! \brief Element tangent matrix indices.
  enum TangentMatrices {
    qxqx =  1, qyqy =  2, qzqz =  3, pp =  4,
    qxqy =  5, qxqz =  6, qxp  =  7,
    qyqx =  8, qyqz =  9, qyp  = 10,
    qzqx = 11, qzqy = 12, qzp  = 13,
    pqx  = 14, pqy  = 15, pqz  = 16,
    NMAT = 17
  };

  using IntVec  = const std::vector<int>;    //!< Convenience type alias
  using UintVec = const std::vector<size_t>; //!< Convenience type alias

public:
  //! \brief The constructor initializes the integrand.
  explicit CompatibleDarcy(unsigned short int n, int torder = 0);

  using Darcy::getLocalIntegral;
  //! \brief Returns a mixed local integral container.
  LocalIntegral* getLocalIntegral(const UintVec& nen,
                                  size_t, bool neumann) const override;

  using Darcy::initElement;
  //! \brief Extracts element vectors from a mixed-basis solution.
  bool initElement(const IntVec& MNPC,
                   const UintVec& elem_sizes, const UintVec& basis_sizes,
                   LocalIntegral& elmInt) override;
  //! \brief Extracts element vectors from a mixed-basis solution.
  bool initElementBou(const IntVec& MNPC,
                      const UintVec& elem_sizes, const UintVec& basis_sizes,
                      LocalIntegral& elmInt) override;

  using Darcy::evalIntMx;
  //! \brief Evaluates the mixed interior integrand.
  bool evalIntMx(LocalIntegral& elmInt, const MxFiniteElement& fe,
                 const Vec3& X) const override;

  using Darcy::evalBouMx;
  //! \brief Evaluates a prescribed normal Darcy flux.
  bool evalBouMx(LocalIntegral& elmInt, const MxFiniteElement& fe,
                 const Vec3& X, const Vec3& normal) const override;

  //! \brief Evaluates secondary fields at a result point.
  bool evalSol(Vector& s, const MxFiniteElement& fe,
               const Vec3& X, const IntVec& MNPC,
               const UintVec& elem_sizes,
               const UintVec& basis_sizes) const override;

  //! \brief Returns the number of primary/secondary fields.
  size_t getNoFields(int fld) const override;

  //! \brief Returns the name of a primary field.
  std::string getField1Name(size_t i, const char* prefix) const override;
  //! \brief Returns the name of a secondary field.
  std::string getField2Name(size_t i, const char* prefix) const override;

  //! \brief No compatible flux-pressure norm integrand is available yet.
  NormBase* getNormIntegrand(AnaSol*) const override { return nullptr; }

protected:
  //! \brief Extracts element solution vectors for current element.
  bool getElementSol(const IntVec& MNPC,
                     const UintVec& elem_sizes, const UintVec& basis_sizes,
                     Vectors& eV, size_t nLevels = 1) const;

private:
  size_t scalarBasis; //!< Index of basis used for pressure
};

#endif
