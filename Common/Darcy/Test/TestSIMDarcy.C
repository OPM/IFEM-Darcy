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

#include "Function.h"
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

  Vec3 perm = darcy.getMaterial().getPermeability(Vec3());
  REQUIRE_THAT(perm[0], WithinRel(98.1));
  REQUIRE_THAT(perm[1], WithinRel(9.81));
  REQUIRE_THAT(perm[2], WithinAbs(0.0, 1e-14));
  Vec3 body = darcy.getBodyForce(Vec3());
  REQUIRE_THAT(body[0], WithinAbs(0.0, 1e-14));
  REQUIRE_THAT(body[1], WithinAbs(0.0, 1e-14));
  REQUIRE_THAT(body[2], WithinAbs(0.0, 1e-14));
  double flux = darcy.getFlux(Vec3(),Vec3());
  REQUIRE_THAT(flux, WithinAbs(0.0, 1e-14));
  double src = darcy.getPotential(Vec3());
  REQUIRE_THAT(src, WithinRel(1.5515174, 1e-7));
  src = darcy.getPotential(Vec3(0.25, 0.25, 0.0));
  REQUIRE_THAT(src, WithinRel(156.157, 1e-6));
}
