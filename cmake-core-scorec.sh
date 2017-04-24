
#flags=" -g -O2"
flags=" -g -O0 -Wall -Wextra -Werror"

cmake .. \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_C_COMPILER="mpicc" \
  -DCMAKE_CXX_COMPILER="mpicxx" \
  -DCMAKE_C_FLAGS="${flags}" \
  -DCMAKE_CXX_FLAGS="${flags}" \
  -DSIM_MPI="mpich3.1.2" \
  -DENABLE_ZOLTAN=ON \
  -DENABLE_SIMMETRIX=True \
  -DSIM_PARASOLID=True \
  -DMDS_ID_TYPE=long \
  -DMDS_SET_MAX:STRING="1024" \
  -DIS_TESTING=True \
  -DPCU_COMPRESS=True \
  -DMESHES="/fasttmp/liuc11/pumi-meshes" \
  ..

