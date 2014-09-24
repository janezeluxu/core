/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef AWRRHS_H
#define AWRRHS_H

#include "apf.h"
#include "apfMesh.h"
#include "apfDynamicMatrix.h"
#include "Teuchos_ParameterList.hpp"

namespace awr {

class RHS
{
  public:
    RHS(apf::Mesh* m, const Teuchos::ParameterList& p);
    virtual ~RHS() {};
    void assemble();
    virtual void 
    evaluateElementRHS(apf::MeshEntity* e,
                       int integration_order,
                       apf::DynamicMatrix& k) = 0;
  protected:
    apf::Mesh* mesh_;
    Teuchos::ParameterList params_;
  private:
    RHS(const RHS&);
    RHS& operator=(const RHS&);
};

}

#endif