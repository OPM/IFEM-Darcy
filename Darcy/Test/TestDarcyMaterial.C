//==============================================================================
//!
//! \file TestDarcyMaterial.C
//!
//! \date May 13 2015
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Tests for darcy material models.
//!
//==============================================================================

#include "DarcyMaterial.h"
#include "Vec3.h"
#include "Vec3Oper.h"
#include "MatVec.h"

#include "Catch2Support.h"
#include <tinyxml2.h>


TEST_CASE("TestDarcyMaterial.Parse")
{
  Vec3 X(1.0,2.0,3.0), K;
  Matrix Kmat;
  tinyxml2::XMLDocument doc;

  doc.Parse("<materialdata>"
            "  <permeability>1.23</permeability>"
	    "</materialdata>");
  DarcyMaterial mat1(doc.RootElement());
  K = mat1.getPermeability(X);
  std::cout <<"K1 = "<< K << std::endl;
  CHECK(K.equal(Vec3(1.23,1.23,1.23)));

  doc.Parse("<materialdata>"
            "  <permeability type=\"diag\">1.23|4.56</permeability>"
	    "</materialdata>");
  DarcyMaterial mat2(doc.RootElement());
  K = mat2.getPermeability(X);
  std::cout <<"K2 = "<< K << std::endl;
  CHECK(K.equal(Vec3(1.23,4.56)));

  doc.Parse("<materialdata>"
            "  <permeability type=\"diag\">1.23*x|4.56*y|7.89*z</permeability>"
	    "</materialdata>");
  DarcyMaterial mat3(doc.RootElement());
  K = mat3.getPermeability(X);
  std::cout <<"K3 = "<< K << std::endl;
  CHECK(K.equal(Vec3(1.23,9.12,23.67)));

  doc.Parse("<materialdata>"
            "  <permeability type=\"matrix\">"
            "     10.0  2.1"
            "      2.3 30.0"
            "  </permeability>"
	    "</materialdata>");
  DarcyMaterial mat4(doc.RootElement());
  REQUIRE(mat4.getPermeability(&Kmat));
  CHECK_THAT(Kmat(1,1), WithinRel(10.0));
  CHECK_THAT(Kmat(1,2), WithinRel(2.1));
  CHECK_THAT(Kmat(2,1), WithinRel(2.3));
  CHECK_THAT(Kmat(2,2), WithinRel(30.0));

  doc.Parse("<materialdata>"
            "  <permeability type=\"matrix\">"
            "     10.0  2.1  0.0"
            "      2.3 30.0  0.2"
            "      0.3  4.5 50.0 99.99"
            "  </permeability>"
	    "</materialdata>");
  DarcyMaterial mat5(doc.RootElement());
  REQUIRE(mat5.getPermeability(&Kmat));
  CHECK_THAT(Kmat(1,1), WithinRel(10.0));
  CHECK_THAT(Kmat(1,2), WithinRel(2.1));
  CHECK_THAT(Kmat(1,3), WithinRel(0.0));
  CHECK_THAT(Kmat(2,1), WithinRel(2.3));
  CHECK_THAT(Kmat(2,2), WithinRel(30.0));
  CHECK_THAT(Kmat(2,3), WithinRel(0.2));
  CHECK_THAT(Kmat(3,1), WithinRel(0.3));
  CHECK_THAT(Kmat(3,2), WithinRel(4.5));
  CHECK_THAT(Kmat(3,3), WithinRel(50.0));
}
