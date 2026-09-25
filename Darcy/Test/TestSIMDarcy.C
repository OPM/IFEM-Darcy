//==============================================================================
//!
//! \file TestSIMDarcy.C
//!
//! \date April 28 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for driver for NURBS-based FEM analysis of of Darcy Flow.
//!
//==============================================================================

#include "Darcy.h"
#include "SIMDarcy.h"
#include "SIM2D.h"
#include "Vec3.h"

#include "Catch2Support.h"


TEST_CASE("TestSIMDarcy.Parse")
{
  Darcy itg(2);
  SIMDarcy<SIM2D> sim(itg);
  REQUIRE(sim.read("Wavefront_k10_p2_b20.xinp"));

  const Darcy& darcy = static_cast<const Darcy&>(*sim.getProblem());

  sim.init();

  Matrix Kinv;
  REQUIRE(darcy.getInvPermeability(Vec3(),Kinv));
  REQUIRE(Kinv.cols() == 2);
  REQUIRE(Kinv.rows() == 2);
  CHECK_THAT(Kinv(1,1), WithinRel(0.1));
  CHECK_THAT(Kinv(2,2), WithinRel(1.0));
  Vec3 body = darcy.getBodyForce(Vec3());
  CHECK_THAT(body[0], WithinAbs(0.0, 1e-14));
  CHECK_THAT(body[1], WithinAbs(1.0, 1e-14));
  CHECK_THAT(body[2], WithinAbs(0.0, 1e-14));
  double flux = darcy.getFlux(Vec3(),Vec3());
  CHECK_THAT(flux, WithinAbs(0.0, 1e-14));
  double src = darcy.getPotential(Vec3());
  CHECK_THAT(src, WithinRel(1.5515174, 1e-7));
  src = darcy.getPotential(Vec3(0.25, 0.25, 0.0));
  CHECK_THAT(src, WithinRel(156.157, 1e-6));
}
