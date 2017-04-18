ZOLTAN_DIR=/lore/seol/mpich3-gcc4.9.2-install
PARMETIS_DIR=/lore/esyoon/openlib/petsc-3.7.5-gcc-4.9.2
PREFIX=/lore/seol/mpich3-gcc4.9.2-install

cmake .. \
  -DCMAKE_C_COMPILER="/usr/local/mpich3/latest/bin/mpicc" \
  -DCMAKE_CXX_COMPILER="/usr/local/mpich3/latest/bin/mpicxx" \
  -DCMAKE_C_FLAGS=" -g -O2 -Wall -Wextra -Werror" \
  -DCMAKE_CXX_FLAGS=" -g -O2 -Wall -Wextra -Werror" \
  -DENABLE_ZOLTAN=ON \
  -DZOLTAN_INCLUDE_DIR="$ZOLTAN_DIR/include" \
  -DZOLTAN_LIBRARY="$ZOLTAN_DIR/lib/libzoltan.a" \
  -DPARMETIS_INCLUDE_DIR="$PARMETIS_DIR/include" \
  -DPARMETIS_LIBRARY="$PARMETIS_DIR/lib/libparmetis.so" \
  -DMETIS_LIBRARY="$PARMETIS_DIR/lib/libmetis.so" \
  -DCMAKE_INSTALL_PREFIX=$PREFIX \
  -DCMAKE_BUILD_TYPE=Debug \
  -DBUILD_EXES=ON \
  -DIS_TESTING=ON \
  -DMESHES=/users/seol/develop/meshes \
  -DMPIRUN=/usr/local/mpich3/latest/bin/mpirun

