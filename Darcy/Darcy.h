// $Id$
//==============================================================================
//!
//! \file Darcy.h
//!
//! \date Mar 27 2015
//!
//! \author Yared Bekele
//!
//! \brief Integrand implementations for Darcy flow problems.
//!
//==============================================================================

#ifndef _DARCY_H_
#define _DARCY_H_

#include "DarcyBase.h"

#include <memory>

class Field;
class RealFunc;
class TractionFunc;
class VecFunc;


/*!
  \brief Class representing the integrand of the Darcy problem.
*/

class Darcy : public DarcyBase
{
protected:
  int pp = 0; //!< Block for pressure
  int cc = 0; //!< Block for concentration
  int cp = 0; //!< Block coupling concentration to pressure

public:
  //! \brief The constructor initializes all pointers to zero.
  explicit Darcy(unsigned short int n, int torder = 0);
  //! \brief The destructor deletes \ref reacInt and \ref bodyforce.
  virtual ~Darcy();

  //! \brief Parses a data section from an XML element.
  bool parse(const tinyxml2::XMLElement* elem) override;

  //! \brief Returns the body force vector at a given point.
  Vec3 getBodyForce(const Vec3& X) const;

  //! \brief Defines the concentration source function.
  virtual void setCSource(RealFunc*) {}

  //! \brief Sets the concentration field from dependent simulation.
  void setNamedField(const std::string& name, Field* field) override;

  //! \brief Defines a scalar flux function.
  void setFlux(RealFunc* f) { flux = f; }
  //! \brief Defines a vectorial flux function.
  void setFlux(VecFunc* f) { vflux = f; }
  //! \brief Defines a vectorial flux function.
  void setFlux(TractionFunc* f) { tflux = f; }

  //! \brief Evaluates the boundary fluid flux (if any) at specified point.
  double getFlux(const Vec3& X, const Vec3& normal) const;
  //! \brief Evaluates the potential source (if any) at specified point.
  double getPotential(const Vec3& X) const;

  //! \brief Defines the global integral for calculating reaction forces only.
  void setSecondaryInt(GlobalIntegral* gq) override;
  //! \brief Returns the system quantity to be integrated by \a *this.
  GlobalIntegral& getGlobalInt(GlobalIntegral* gq) const override;

  using DarcyBase::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] time Time stepping parameters
  //! \param[in] X Cartesian coordinates of current integration point
  bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
               const TimeDomain& time, const Vec3& X) const override;

  using DarcyBase::evalBou;
  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  bool evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
               const Vec3& X, const Vec3& normal) const override;

  //! \brief Evaluates the secondary solution at a result point.
  //! \param[out] s Array of solution field values at current point
  //! \param[in] eV Element solution vectors
  //! \param[in] fe Finite element data at current point
  //! \param[in] X Cartesian coordinates of current point
  bool evalSol2(Vector& s, const Vectors& eV,
                const FiniteElement& fe, const Vec3& X) const override;

  //! \brief Returns the number of primary/secondary solution field components.
  //! \param[in] fld which field set to consider (1=primary, 2=secondary)
  size_t getNoFields(int fld) const override { return fld > 1 ? 2*nsd+2 : 1; }

  //! \brief Returns the name of the primary solution field.
  //! \param[in] prefix Name prefix
  std::string getField1Name(size_t, const char* prefix) const override;

  //! \brief Returns the name of a secondary solution field component.
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  std::string getField2Name(size_t i, const char* prefix) const override;

  //! \brief Returns a pointer to an Integrand for solution norm evaluation.
  //! \param[in] asol Pointer to analytical solution (optional)
  //!
  //! \note The Integrand object is allocated dynamically and has to be deleted
  //! manually when leaving the scope of the pointer variable receiving the
  //! returned pointer value.
  NormBase* getNormIntegrand(AnaSol* asol) const override;

  //! \brief Returns pressure in a point.
  //! \param[out] eV Element solution vectors
  //! \param[in] fe Finite element data at current point
  //! \param[in] level Time level to evaluate at
  virtual double pressure(const Vectors& eV, const FiniteElement& fe,
                          size_t level) const;

  //! \brief Returns pressure gradient in a point.
  //! \param[in] eV Element solution vectors
  //! \param[in] fe Finite element data at current point
  Vec3 pressureGradient(const Vectors& eV, const FiniteElement& fe) const;

  //! \brief Returns the darcy velocity in a point.
  //! \param[out] q Darcy velocity vector
  //! \param[in] eV Element solution vectors
  //! \param[in] fe Finite element data at current point
  //! \param[in] X Cartesian coordinates of current point
  bool evalDarcyVel(Vector& q, const Vectors& eV,
                    const FiniteElement& fe, const Vec3& X) const;

  //! \brief Returns the fluid mass density in a point.
  //! \param[in] fe Finite element data at current point
  double getDensity(const FiniteElement& fe) const;

protected:
  VecFunc*      bodyforce; //!< Body force function
  RealFunc*     flux;      //!< Flux function
  VecFunc*      vflux;     //!< Flux vector function
  TractionFunc* tflux;     //!< Flux traction function

  std::unique_ptr<RealFunc> source; //!< Source function
  std::unique_ptr<Field>    cField; //!< Tracer concentration field
  Vector                    cVec;   //!< Tracer concentration values

  GlobalIntegral* reacInt; //!< Reaction-forces-only integral

public:
  char extEner; //!< If \e true, external energy is to be computed
};


/*!
  \brief Class representing the integrand of Darcy energy norms.
*/

class DarcyNorm : public NormBase
{
public:
  //! \brief Enumeration of regular norm entries
  enum NormEntries {
    H1_Ph = 0,
    EXT_ENERGY,
    H1_Ch,
    H1_P,
    H1_E_Ph,
    H1_C,
    H1_E_Ch,
    TOTAL_NORM_E
  };

  //! \brief Enumeration of recovery norm entries
  enum RecoveryEntries {
    H1_Pr = 0,
    H1_Pr_Ph,
    H1_Cr,
    H1_Cr_Ch,
    TOTAL_NORM_REC,
    H1_E_Pr,
    H1_E_Cr,
    TOTAL_E_REC,
    EFF_REC_Ph,
    EFF_REC_Ch,
    EFF_REC_TOTAL
  };

  //! \brief The only constructor initializes its data members.
  //! \param[in] p The Darcy problem to evaluate norms for
  //! \param[in] a The analytical darcy flux (optional)
  explicit DarcyNorm(Darcy& p, VecFunc* a = nullptr);
  //! \brief Empty destructor.
  virtual ~DarcyNorm();

  using NormBase::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
               const Vec3& X) const override;

  using NormBase::evalBou;
  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite Element quantities
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  bool evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
               const Vec3& X, const Vec3& normal) const override;

  using NormBase::finalizeElement;
  //! \brief Finalizes the element norms after the numerical integration.
  //! \param elmInt The local integral object to receive the contributions
  bool finalizeElement(LocalIntegral& elmInt) override;

  //! \brief Returns whether this norm has explicit boundary contributions.
  bool hasBoundaryTerms() const override { return true; }

  //! \brief Returns the number of norm groups or size of a specified group.
  //! \param[in] group The norm group to return the size of
  //! (if zero, return the number of groups)
  size_t getNoFields(int group) const override;

  //! \brief Returns the name of a norm quantity.
  //! \param[in] i The norm group (one-based index)
  //! \param[in] j The norm number (one-based index)
  //! \param[in] prefix Common prefix for all norm names
  std::string getName(size_t i, size_t j, const char* prefix) const override;

  //! \brief Returns whether a norm quantity stores element contributions.
  //! \param[in] i The norm group (one-based index)
  //! \param[in] j The norm number (one-based index)
  bool hasElementContributions(size_t i, size_t j) const override;

protected:
  VecFunc* anasol; //!< Analytical darcy flux field
};

#endif
