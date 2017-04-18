#include <apf.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfZoltan.h>
#include <gmi_mesh.h>
#include <parma.h>
#include <PCU.h>
#include <pcu_util.h>

int main(int argc, char** argv)
{
  PCU_ALWAYS_ASSERT(argc == 4);
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_register_mesh();
  //load model and mesh
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);
  apf::MeshTag* weights = Parma_WeighByMemory(m);
  apf::Balancer* balancer = makeZoltanBalancer(m, apf::GRAPH, apf::REPARTITION);
  balancer->balance(weights, 1.10);
  delete balancer;
  apf::removeTagFromDimension(m, weights, m->getDimension());
  Parma_PrintPtnStats(m, "");
  m->destroyTag(weights);
  m->writeNative(argv[3]);
  // destroy mds
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}
