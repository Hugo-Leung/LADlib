#!/bin/bash
set -e 

cd /opt

version="2.4.7.1"
fname="clhep-"${version}
url="https://proj-clhep.web.cern.ch/proj-clhep/dist1/"${fname}".tgz"

wget ${url}
tar zxf "${fname}.tgz"

cd ${version}

mkdir build 
cd build 

cmake -DCMAKE_INSTALL_PREFIX=/opt/CLHEP ../CLHEP

njobs=5
ncpu=$(nproc)
if [[ ${ncpu} -lt ${njobs} ]]; then
	njobs=${cpu}
fi

make -j${njobs} install

cd /opt	
rm -rf ${version}
rm "${fname}.tgz"

echo /opt/CLHEP >> /etc/ld.so.conf
ldconfig
