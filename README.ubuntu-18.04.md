# Build instructions for Ubuntu 18.04
First get gfortran9
```
sudo add-apt-repository -y ppa:ubuntu-toolchain-r/test
sudo apt update
sudo apt -y install gfortran-9
```
Then make gfortran-9 be the default gfortran (so we don't have to change the makefile)
```
sudo update-alternatives --install /usr/bin/gfortran gfortran /usr/bin/gfortran-7 7
sudo update-alternatives --install /usr/bin/gfortran gfortran /usr/bin/gfortran-9 9
```
Install dependencies
```
sudo apt -y install nco
sudo apt -y install netcdf-bin libnetcdff-dev pgplot5
```
Build & test
```
make -f makefile_gfortran64 clean
make -f makefile_gfortran64
ctest
```

