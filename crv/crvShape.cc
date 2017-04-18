/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include <PCU.h>
#include "crv.h"
#include "crvAdapt.h"
#include "crvShape.h"
#include "crvTables.h"
#include <maCoarsen.h>
#include <maEdgeSwap.h>
#include <maOperator.h>
#include <pcu_util.h>

/* This is similar to maShape.cc, conceptually, but different enough
 * that some duplicate code makes sense */

namespace crv {

bool isBoundaryEntity(apf::Mesh* m, apf::MeshEntity* e)
{
  return m->getModelType(m->toModel(e)) < m->getDimension();
}
/** \brief checks if any entity has two entities of
 * dimension on the boundary
 * \details this is useful for some shape correction assessments,
 * and in general, curved elements with multiple entities on the boundary
 * are at risk for poor quality, since this strongly constrains
 * their shape */
static bool hasTwoEntitiesOnBoundary(apf::Mesh* m, apf::MeshEntity* e, int dimension)
{
  apf::Downward down;
  int count = 0;
  int nd = m->getDownward(e,dimension,down);
  for (int i = 0; i < nd; ++i){
    if(isBoundaryEntity(m,down[i]))
      ++count;
    if(count == 2)
      return true;
  }
  return false;
}

/* Mark Edges based on the invalidity code the element has been
 * tagged with.
 */
static int markEdges(ma::Mesh* m, ma::Entity* e, int tag,
    ma::Entity* edges[6])
{
  if ( tag <= 1 ) // if its valid, or not checked, don't worry about it
    return 0;
  int dim = (tag-2)/6;
  int index = (tag-2) % 6;
  int n = 0;
  int md = m->getDimension();

  switch (dim) {
    case 0:
    {
      // if we have an invalid vertex, operate on its edges
      ma::Downward ed;
      m->getDownward(e,1,ed);
      n = md;
      if(md == 2){
        edges[0] = ed[index];
        edges[1] = ed[(index+2) % 3];
      } else {
        PCU_ALWAYS_ASSERT(index < 4);
        edges[0] = ed[vertEdges[index][0]];
        edges[1] = ed[vertEdges[index][1]];
        edges[2] = ed[vertEdges[index][2]];
      }
    }
    break;
    case 1:
    {
      // if we have a single invalid edge, operate on it
      ma::Downward ed;
      m->getDownward(e,1,ed);
      edges[0] = ed[index];
      n = 1;
    }
    break;
    case 2:
    {
      // if we have an invalid face, operate on its edges
      ma::Downward ed, faces;
      m->getDownward(e,2,faces);
      m->getDownward(faces[index],1,ed);
      n = 3;
      edges[0] = ed[0];
      edges[1] = ed[1];
      edges[2] = ed[2];
    }
    break;
    case 3:
      m->getDownward(e,1,edges);
      n = 6;
      break;
    default:
      fail("invalid quality tag in markEdges\n");
      break;
  }

  return n;
}

class EdgeSwapper : public ma::Operator
{
public:
  EdgeSwapper(Adapt* a)
  {
    adapter = a;
    mesh = a->mesh;
    edges[0] = edges[1] = edges[2] = edges[3] = edges[4] = edges[5] = 0;
    simplex = 0;
    edgeSwap = ma::makeEdgeSwap(a);
    md = mesh->getDimension();
    ne = ns = 0;
  }
  virtual ~EdgeSwapper()
  {
    delete edgeSwap;
  }
  virtual int getTargetDimension() {return md;}
  virtual bool shouldApply(ma::Entity* e)
  {
    int tag = crv::getTag(adapter,e);
    ne = markEdges(mesh,e,tag,edges);
    simplex = e;
    return (ne > 0);
  }
  virtual bool requestLocality(apf::CavityOp* o)
  {
    return o->requestLocality(edges,ne);
  }
  virtual void apply()
  {
    for (int i = 0; i < ne; ++i){
      if (edgeSwap->run(edges[i])){
        ns++;
        crv::clearTag(adapter,simplex);
        ma::clearFlag(adapter,edges[i],ma::COLLAPSE | ma::BAD_QUALITY);
        break;
      }
    }
  }
private:
  Adapt* adapter;
  ma::Mesh* mesh;
  ma::Entity* simplex;
  ma::Entity* edges[6];
  ma::EdgeSwap* edgeSwap;
  int md;
  int ne;
public:
  int ns;
};

class EdgeReshaper : public ma::Operator
{
public:
  EdgeReshaper(Adapt* a)
  {
    adapter = a;
    mesh = a->mesh;
    edges[0] = edges[1] = edges[2] = edges[3] = edges[4] = edges[5] = 0;
    simplex = 0;
    md = mesh->getDimension();
    ne = nr = 0;
    qual = makeQuality(mesh,2);
  }
  virtual ~EdgeReshaper()
  {
    delete qual;
  }
  virtual int getTargetDimension() {return md;}
  virtual bool shouldApply(ma::Entity* e)
  {
    int tag = crv::getTag(adapter,e);
    ne = markEdges(mesh,e,tag,edges);
    simplex = e;
    return (ne > 0);
  }
  virtual bool requestLocality(apf::CavityOp* o)
  {
    return o->requestLocality(edges,ne);
  }
  virtual void apply()
  {
    for (int i = 0; i < ne; ++i){
      if (!isBoundaryEntity(mesh,edges[i]) &&
          repositionEdge(edges[i])){
        nr++;
        crv::clearTag(adapter,simplex);
        ma::clearFlag(adapter,edges[i],ma::COLLAPSE | ma::BAD_QUALITY);
        break;
      }
    }
  }
private:
  /** \brief reposition second order edge control point based on XJ Luo's
      thesis and bezier.tex in SCOREC/docs repo, only works for second order */
  bool repositionEdge(ma::Entity* edge)
  {
    // lets assume we have an edge we want to fix
    // only support second order for now
    int P = mesh->getShape()->getOrder();
    if (P != 2) return false;

    ma::Entity* verts[4];
    ma::Entity* edges[6];

    ma::Vector pivotPoint;
    ma::Vector edgeVectors[3];
    mesh->getDownward(simplex,0,verts);
    mesh->getDownward(simplex,1,edges);

    // pick a pivotVert, the vertex with the worse jacobian determinant
    ma::Entity* pivotVert;
    int pivotIndex;
    {
      apf::MeshElement* me = apf::createMeshElement(mesh,simplex);

      ma::Entity* edgeVerts[2];
      mesh->getDownward(edge,0,edgeVerts);
      apf::Matrix3x3 J;
      pivotIndex = apf::findIn(verts,4,edgeVerts[0]);
      PCU_ALWAYS_ASSERT(pivotIndex >= 0);

      ma::Vector xi = crv::elem_vert_xi[apf::Mesh::TET][pivotIndex];
      apf::getJacobian(me,xi,J);

      double j = apf::getJacobianDeterminant(J,3);
      pivotVert = edgeVerts[0];

      int index = apf::findIn(verts,4,edgeVerts[1]);
      PCU_ALWAYS_ASSERT(index >= 0);
      xi = crv::elem_vert_xi[apf::Mesh::TET][index];
      apf::getJacobian(me,xi,J);
      if (apf::getJacobianDeterminant(J,3) < j){
        pivotVert = edgeVerts[1];
        pivotIndex = index;
      }
      apf::destroyMeshElement(me);
    }

    mesh->getPoint(pivotVert,0,pivotPoint);

    // local, of edges around vert, [0,2]
    int edgeIndex = 0;

    for (int i = 0; i < 3; ++i){
      // theres only one point, so reuse this...
      edgeVectors[i] = ma::getPosition(mesh,edges[vertEdges[pivotIndex][i]])
                     - pivotPoint;
      if (edges[vertEdges[pivotIndex][i]] == edge)
        edgeIndex = i;
    }
    PCU_ALWAYS_ASSERT(edgeIndex >= 0);
    ma::Vector normal = apf::cross(edgeVectors[(1+edgeIndex) % 3],
        edgeVectors[(2+edgeIndex) % 3]);
    double length = normal.getLength();
    double validity = edgeVectors[edgeIndex]*normal;

    if(validity > 1e-10)
      return false;
    ma::Vector oldPoint = ma::getPosition(mesh,edge);
    apf::Adjacent adjacent;
    mesh->getAdjacent(edge,3,adjacent);

    // places the new point at a 20 degree angle with the plane
    double angle = apf::pi/9.;
    ma::Vector newPoint = edgeVectors[edgeIndex] + pivotPoint
        + normal/length*(-validity/length +
            edgeVectors[edgeIndex].getLength()*sin(angle));

    mesh->setPoint(edge,0,newPoint);

    for (std::size_t i = 0; i < adjacent.getSize(); ++i){
      if (qual->checkValidity(adjacent[i]) > 0){
        mesh->setPoint(edge,0,oldPoint);
        return false;
      }
    }

    return true;
  }
  Adapt* adapter;
  ma::Mesh* mesh;
  Quality* qual;
  ma::Entity* simplex;
  ma::Entity* edges[6];
  int md;
  int ne;
public:
  int nr;
};

static bool isCornerTriAngleLarge(crv::Adapt *a,
    ma::Entity* tri, int index)
{
  ma::Mesh* m = a->mesh;
  apf::Element* elem = apf::createElement(m->getCoordinateField(),tri);
  apf::NewArray<apf::Vector3> nodes;
  apf::getVectorNodes(elem,nodes);
  apf::destroyElement(elem);

  ma::Vector normal = ma::getTriNormal(m,tri);

  int P = m->getShape()->getOrder();
  int r = index*(P-1)+3; // index to the right
  int l = ((index+2) % 3)*(P-1)+3+P-2; // index to the left

  ma::Vector cornerNormal = apf::cross((nodes[r]-nodes[index]),
      (nodes[l]-nodes[index]));

  // this statement is not exactly a fair comparison, but at least gives
  // some level of control over what is considered an invalid angle
  // an angle is "too large" if the dot product between the corner triangle
  // and the triangle formed by its vertices is negative
  if (cornerNormal*normal < a->input->validQuality)
    return true;

  // one final check to handle the odd case where one of the two tets shared
  // by the angle is invalid at the vertex between the two edges,
  // but all of the faces do not have large angles
  // check both its tets are okay
  if (m->getDimension() == 3 && !isBoundaryEntity(m,tri)){
    ma::Entity* triVerts[3];
    m->getDownward(tri,0,triVerts);
    apf::Up up;
    m->getUp(tri,up);

    apf::Matrix3x3 J;
    for (int i = 0; i < up.n; ++i){
      apf::MeshElement* me = apf::createMeshElement(m,up.e[i]);
      ma::Entity* verts[4];
      m->getDownward(up.e[i],0,verts);

      int tetIndex = apf::findIn(verts,4,triVerts[index]);
      ma::Vector xi = crv::elem_vert_xi[apf::Mesh::TET][tetIndex];
      apf::getJacobian(me,xi,J);
      apf::destroyMeshElement(me);

      if(apf::getJacobianDeterminant(J,3) < a->input->validQuality){
        return true;
      }
    }
  }
  return false;
}

/* Checks if an angle of a triangle is large (>= 180 degrees)
 * which can be caused by two edges on the boundary curved to it
 *
 */

static ma::Entity* isLargeAngleTri(crv::Adapt* a, ma::Entity* e)
{
  ma::Mesh* m = a->mesh;
  ma::Entity* edges[3];
  m->getDownward(e,1,edges);
  for (int i = 0; i < 3; ++i)
  {
    ma::Entity* e0 = edges[i];
    ma::Entity* e1 = edges[(i+1) % 3];
    if(isBoundaryEntity(m,e0) && isBoundaryEntity(m,e1))
    {
      if(isCornerTriAngleLarge(a,e,(i+1) % 3)){
        ma::Entity* edge = edges[(i+2) % 3];
        if(!ma::getFlag(a,edge,ma::SPLIT) && !isBoundaryEntity(m,edge)){
          return edge;
        }
      }
    }
  }
  return 0;
}

/* Checks if an angle of a tet is large (>= 180 degrees)
 * which can be caused by two edges on the boundary curved to it
 *
 * An analytic approach, looking at the control point net of points
 * by comparing surface normals of each pair of control net points
 * adjacent to the triangle is an incredibly complex ordering exercise
 * Rather than attempt to do the ordering, sampling the Jacobian at
 * P+1 points is used. A validity check on this edge could also be used
 */

static ma::Entity* isLargeAngleTet(crv::Adapt* a, ma::Entity* e)
{
  ma::Mesh* m = a->mesh;

  ma::Entity* faces[4];
  m->getDownward(e,2,faces);

  int index = -1;
  int P = m->getShape()->getOrder();
  // find edge that matters
  for (int i = 0; i < 6; ++i){
    if(isBoundaryEntity(m,faces[edgeFaces[i][0]]) &&
        isBoundaryEntity(m,faces[edgeFaces[i][1]])){
      index = i;
      break;
    }
  }
  if(index < 0) return 0;

  apf::FieldShape* fs = m->getShape();
  ma::Entity* edges[6];
  m->getDownward(e,1,edges);
  if(!isBoundaryEntity(m,edges[index])) return 0;
  ma::Entity* edge = 0;

  // lets do a sampling approach. At each point on the edge
  int bt = apf::Mesh::EDGE;
  ma::Entity* leftFace = faces[edgeFaces[index][0]];
  ma::Entity* rightFace = faces[edgeFaces[index][1]];
  apf::MeshElement* leftMe = apf::createMeshElement(m,leftFace);
  apf::MeshElement* rightMe = apf::createMeshElement(m,rightFace);
  apf::Vector3 leftXi,rightXi;

  apf::NewArray<apf::Vector3> nodeXi(P+1);
  nodeXi[0] = apf::Vector3(-1,0,0);
  for (int i = 0; i < P-1; ++i)
    fs->getNodeXi(bt,i,nodeXi[i+1]);
  nodeXi[P] = apf::Vector3(1,0,0);

  for (int i = 0; i < P+1; ++i){
    leftXi = apf::boundaryToElementXi(m,edges[index],leftFace,nodeXi[i]);
    rightXi = apf::boundaryToElementXi(m,edges[index],rightFace,nodeXi[i]);
    apf::Matrix3x3 leftJ,rightJ;

    apf::getJacobian(leftMe,leftXi,leftJ);
    apf::getJacobian(rightMe,rightXi,rightJ);
    apf::Vector3 leftN = (apf::cross(leftJ[0],leftJ[1])).normalize();
    apf::Vector3 rightN = (apf::cross(rightJ[0],rightJ[1])).normalize();
    // jacobian has two rows, each with a vector
    // these two vectors form the plane
    // compare their directions to see if they are close to coplanar)
    // 10 degree tolerance
    if(std::fabs(leftN*rightN) > 0.9){
      edge = edges[oppEdges[index]];
      break;
    }
  }
  apf::destroyMeshElement(leftMe);
  apf::destroyMeshElement(rightMe);

  return edge;
}

// BAD_QUALITY flag is used on edges to identify them
// as splits for quality, rather than for size refinement
// These two functions handle two seperate situations

// First, triangles are looked at to see if they have an angle > 180 degrees
// There is also a check for the weird situation described above
// Second, tets are looked at to see if they have two faces on the boundary
// where the jacobian determinant is negative along the shared edge,
// indicative of a large angle (curving around a cylinder or sphere).

static int markEdgesOppLargeAnglesTri(Adapt* a)
{
  int count = 0;
  int prev_count;

  ma::Mesh* m = a->mesh;
  ma::Entity* e;

  do {
    ma::Iterator* it = m->begin(2);
    prev_count = count;
    while ((e = m->iterate(it)))
    {
      if(!hasTwoEntitiesOnBoundary(m,e,1)) continue;
      ma::Entity* edge = isLargeAngleTri(a,e);
      if (edge)
      {
        PCU_ALWAYS_ASSERT(m->getType(edge) == 1);
        ma::setFlag(a,edge,ma::SPLIT);
        ma::setFlag(a,edge,ma::BAD_QUALITY);
        if (a->mesh->isOwned(edge))
          ++count;
      }
    }
    m->end(it);
  } while(count > prev_count);
  return PCU_Add_Long(count);
}

static int markEdgesOppLargeAnglesTet(Adapt* a)
{
  int count = 0;
  int prev_count;

  ma::Mesh* m = a->mesh;
  ma::Entity* e;

  do {
    ma::Iterator* it = m->begin(3);
    prev_count = count;
    while ((e = m->iterate(it)))
    {
      ma::Entity* edge = isLargeAngleTet(a,e);
      if (edge && !ma::getFlag(a,edge,ma::SPLIT))
      {
        PCU_ALWAYS_ASSERT(m->getType(edge) == 1);
        ma::setFlag(a,edge,ma::SPLIT);
        ma::setFlag(a,edge,ma::BAD_QUALITY);
        if (a->mesh->isOwned(edge))
          ++count;
      }
    }
    m->end(it);
  } while(count > prev_count);
  return PCU_Add_Long(count);
}

/* The whole idea is to do the quality check once,
 * and then use the results to mark edges, etc for
 * fixing
 */
static int markEdgesToFix(Adapt* a, int flag)
{
  // do a invalidity check first
  int invalid = markInvalidEntities(a);
  if ( !invalid )
    return 0;
  int count = 0;

  ma::Mesh* m = a->mesh;
  ma::Entity* e;
  ma::Entity* edges[3];
  ma::Iterator* it = m->begin(m->getDimension());
  while ((e = m->iterate(it)))
  {
    int tag = crv::getTag(a,e);
    int n = markEdges(m,e,tag,edges);
    for (int i = 0; i < n; ++i){
      ma::Entity* edge = edges[i];
      PCU_ALWAYS_ASSERT(edge);
      if (edge && !ma::getFlag(a,edge,flag))
      {
        ma::setFlag(a,edge,flag);
        if (a->mesh->isOwned(edge))
          ++count;
      }
    }
  }
  m->end(it);

  return PCU_Add_Long(count);
}

int fixLargeBoundaryAngles(Adapt* a)
{
  double t0 = PCU_Time();
  int count = markEdgesOppLargeAnglesTet(a);
  count += markEdgesOppLargeAnglesTri(a);

  if ( ! count){
    return 0;
  }
  splitEdges(a);
  double t1 = PCU_Time();
  ma::print("split %d boundary edges with "
      "large angles in %f seconds",count,t1-t0);
  return 0;
}

static void collapseInvalidEdges(Adapt* a)
{
  double t0 = PCU_Time();
  ma::Mesh* m = a->mesh;
  int maxDimension = m->getDimension();
  PCU_ALWAYS_ASSERT(checkFlagConsistency(a,1,ma::COLLAPSE));
  int successCount = 0;
  for (int modelDimension=1; modelDimension <= maxDimension; ++modelDimension)
  {
    checkAllEdgeCollapses(a,modelDimension);
    findIndependentSet(a);
    successCount += ma::collapseAllEdges(a, modelDimension);
  }
  successCount = PCU_Add_Long(successCount);
  double t1 = PCU_Time();
  ma::print("Collapsed %d bad edges "
      "in %f seconds",successCount, t1-t0);
}

static void swapInvalidEdges(Adapt* a)
{
  double t0 = PCU_Time();
  EdgeSwapper es(a);
  ma::applyOperator(a,&es);
  double t1 = PCU_Time();
  ma::print("Swapped %d bad edges "
      "in %f seconds",es.ns, t1-t0);
}

static void repositionInvalidEdges(Adapt* a)
{
  double t0 = PCU_Time();
  EdgeReshaper es(a);
  ma::applyOperator(a,&es);
  double t1 = PCU_Time();
  ma::print("Repositioned %d bad edges "
      "in %f seconds",es.nr, t1-t0);
}

int fixInvalidEdges(Adapt* a)
{
  int count = markEdgesToFix(a,ma::BAD_QUALITY | ma::COLLAPSE );
  if (! count){
    return 0;
  }

  if(a->mesh->getShape()->getOrder() == 2)
    repositionInvalidEdges(a);
  collapseInvalidEdges(a);
  swapInvalidEdges(a);
  return count;
}

}
