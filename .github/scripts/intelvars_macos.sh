source /opt/intel/oneapi/setvars.sh
sudo mkdir -p /usr/local/bin
sudo ln -s "$(command -v ifort)" /usr/local/bin
{
echo ONEAPI_ROOT="/opt/intel/oneapi"
echo IFORT_COMPILER21="$(dirname "$(dirname "$(dirname "$(command -v ifort)")")")"
echo PATH='$PATH':"$(dirname "$(command -v ifort)")"
} >> "$GITHUB_ENV"
