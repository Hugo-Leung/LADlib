From docker.io/library/almalinux:9-minimal
ARG DEF_HASH=""

COPY  ./dnf.sh /opt/dnf.sh
COPY ./packages /opt/packages
COPY ./hcana.sh /opt/hcana.sh
COPY ./ladlib.sh /opt/ladlib.sh

RUN <<EOF
set -ex
mkdir /w
mkdir /volatile
mkdir /cache
mkdir /scratch
mkdir /cvmfs
ln -s /w/work /work
EOF

ARG CACHEBUST=0
RUN <<EOF
curl -fSL https://pki.jlab.org/JLabCA.crt -o  /etc/pki/ca-trust/source/anchors/JLabCA.crt
update-ca-trust
bash /opt/dnf.sh
bash /opt/hcana.sh
bash /opt/ladlib.sh

rm /opt/dnf.sh  /opt/hcana.sh /opt/ladlib.sh /opt/packages 
EOF

ENV PATH=/opt/hcana/bin:$PATH
ENV LD_LIBRARY_PATH=/opt/LADlib/lib64/:$LD_LIBRARY_PATH
ENV ROOT_INCLUDE_PATH=/opt/LADlib/include:/opt/hcana/include:$ROOT_INCLUDE_PATH


LABEL DEF_HASH="${DEF_HASH}"

RUN set -ex; hcana -v && echo "hcana is running successfully." && \
    hcana -l -b -q -e "gSystem->Load(\"libLAD\")" && \
    echo "LADlib is loaded successfully."

ENTRYPOINT ["hcana"]


