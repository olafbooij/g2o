// g2o - General Graph Optimization
// Copyright (C) 2011 R. Kuemmerle, G. Grisetti, W. Burgard
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright
//   notice, this list of conditions and the following disclaimer in the
//   documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
// TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "gtest/gtest.h"

#include "unit_test/test_helper/evaluate_jacobian.h"

#include "g2o/types/sclam2d/edge_se2_sensor_calib.h"
#include "g2o/stuff/os_specific.h"

using namespace std;
using namespace g2o;
using namespace Eigen;

static SE2 randomSE2()
{
  return SE2(Vector3d::Random());
}

TEST(Sclam2D, EdgeSE2SensorCalibJacobian)
{
  VertexSE2 v1;
  v1.setId(0); 

  VertexSE2 v2;
  v2.setId(1); 

  VertexSE2 v3;
  v3.setId(2); 

  EdgeSE2SensorCalib e;
  e.setVertex(0, &v1);
  e.setVertex(1, &v2);
  e.setVertex(2, &v3);
  e.setInformation(EdgeSE2::InformationType::Identity());

  JacobianWorkspace jacobianWorkspace;
  JacobianWorkspace numericJacobianWorkspace;
  numericJacobianWorkspace.updateSize(&e);
  numericJacobianWorkspace.allocate();

  //for (int k = 0; k < 10000; ++k) {
  for (int k = 0; k < 1; ++k) {
    v1.setEstimate(randomSE2());
    v2.setEstimate(randomSE2());
    v3.setEstimate(randomSE2());
    e.setMeasurement(randomSE2());

    evaluateJacobian(e, jacobianWorkspace, numericJacobianWorkspace);
  }
}
