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

#include <string>

#include "gtest/gtest.h"

TEST(TestSIMDarcy, Parse)
{
  Darcy itg(2);
  SIMDarcy<SIM2D> sim(itg);
  EXPECT_TRUE(sim.read("Wavefront_k10_p2_b20.xinp"));

  const Darcy& darcy = static_cast<const Darcy&>(*sim.getProblem());

  Vec3 perm = darcy.getPermeability(Vec3());
  ASSERT_FLOAT_EQ(perm[0], 98.1);
  ASSERT_FLOAT_EQ(perm[1], 9.81);
  ASSERT_FLOAT_EQ(perm[2], 0.0);
  Vec3 body = darcy.getBodyForce(Vec3());
  ASSERT_FLOAT_EQ(body[0], 0.0);
  ASSERT_FLOAT_EQ(body[1], 0.0);
  ASSERT_FLOAT_EQ(body[2], 0.0);
  double flux = darcy.getFlux(Vec3(),Vec3());
  ASSERT_FLOAT_EQ(flux, 0.0);
  double src = darcy.getPotential(Vec3());
  ASSERT_FLOAT_EQ(src, 1.5515174);
  src = darcy.getPotential(Vec3(0.25, 0.25, 0.0));
  ASSERT_FLOAT_EQ(src, 156.157);
}
