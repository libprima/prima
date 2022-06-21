#!/bin/bash

# do the job in /tmp
cd /tmp || exit 1

# use wget to fetch the Intel repository public key
wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB

# add to your apt sources keyring so that archives signed with this key will be trusted.
apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB

# optionally, remove the public key
rm GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB

# the installation
echo "deb https://apt.repos.intel.com/oneapi all main" | tee /etc/apt/sources.list.d/oneAPI.list
add-apt-repository "deb https://apt.repos.intel.com/oneapi all main"
apt-get update
#apt install intel-basekit intel-hpckit  # Instead of this line, the following line seems to suffice
apt-get install -y intel-oneapi-common-vars intel-oneapi-compiler-fortran

# source `setvars.sh`
if [[ -f /opt/intel/oneapi/setvars.sh ]] ; then
    source /opt/intel/oneapi/setvars.sh > /dev/null ; env | grep intel | sed "s/:\/User.*$//" \
        | sed 's/^/export /' > "/tmp/setenv.sh"
    source "/tmp/setenv.sh"
else
    exit 2
fi

IFORT="$(find /opt/*ntel* -type f -executable -name ifort -print -quit)"
IFX="$(find /opt/*ntel* -type f -executable -name ifx -print -quit)"
ln -s "$IFORT" /usr/bin
ln -s "$IFX" /usr/bin
