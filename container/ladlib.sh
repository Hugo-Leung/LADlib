#!/bin/bash
set -ex 
cd /opt

export PATH=/opt/hcana/bin:$PATH
export LD_LIBRARY_PATH=/opt/hcana/lib64:$LD_LIBRARY_PATH

git clone --depth=1 https://github.com/JeffersonLab/LADlib.git LADlib-src
cd LADlib-src
mkdir build 

njobs=5
ncpu=$(nproc)
if [[ ${ncpu} -lt ${njobs} ]]; then
	njobs=${cpu}
fi

cmake -B build -S . -DCMAKE_INSTALL_PREFIX=/opt/LADlib
cmake --build build -j${njobs}
cmake --install build

echo /opt/LADlib/lib64 >> /etc/ld.so.conf
ldconfig

