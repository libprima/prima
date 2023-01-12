#!/bin/bash
# This script installs the Fortran compilers provided in Intel OneAPI.
# Usage: sudo bash install_oneapi_linux.sh

# do the job in the temporary directory of the system
cd /tmp || exit 42

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
installer_exit_code=$?

# Run the script that sets the environment variables.
source /opt/intel/oneapi/setvars.sh

# Show the result of the installation.
echo "The latest ifort installed is:"
ifort --version
echo "The path to ifort is:"
command -v ifort
echo "The latest ifx installed is:"
ifx --version
echo "The path to ifx is:"
command -v ifx

# Exit with the installer exit code.
exit $installer_exit_code
