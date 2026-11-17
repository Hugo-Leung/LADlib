#!/bin/bash
set -ex 
cd /opt
export PATH=/opt/CLHEP/bin:$PATH

version="v11.4.0"
fname="geant4-${version}"

url="https://gitlab.cern.ch/geant4/geant4/-/archive/"${version}"/"${fname}".tar.gz"

wget ${url}

tar zxfp "${fname}.tar.gz"

cd ${fname}
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=/opt/geant4 \
 -DGEANT4_INSTALL_DATA=ON \
 -DGEANT4_USE_QT=ON \
 -DGEANT4_USE_OPENGL_X11=ON \
 -DGEANT4_INSTALL_EXAMPLES=OFF \
 -DGEANT4_BUILD_MULTITHREADED=ON \
 -DGEANT4_USE_GDML=ON \
 -DGEANT4_USE_SYSTEM_CLHEP=ON \
 /opt/${fname}


njobs=5
ncpu=$(nproc)
if [[ ${ncpu} -lt ${njobs} ]]; then
	njobs=${cpu}
fi
make -j${njobs} install

cd /opt
rm -rf ${fname}
rm "${fname}.tar.gz"
echo /opt/geant4 >> /etc/ld.so.conf
ldconfig
