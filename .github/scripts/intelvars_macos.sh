source /opt/intel/oneapi/setvars.sh
sudo mkdir -p /usr/local/bin
sudo ln -s $(which ifort) /usr/local/bin
sudo ln -s $(which ifx) /usr/local/bin
echo "ONEAPI_ROOT=/opt/intel/oneapi" >> "$GITHUB_ENV"
echo "IFORT_COMPILER21=$(dirname $(dirname $(dirname $(which ifort))))" >> "$GITHUB_ENV"
echo "PATH=$PATH:$(dirname $(which ifort)):$(dirname $(which ifx))"
