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
sudo apt -y install netcdf-bin
```
Build
```
make -f makefile_gfortran64 clean
make -f makefile_gfortran64
```
Run some tests
```
cd 00_Genray_Regression_Tests
mkdir test-CMod-LH-edge
cd test-CMod-LH-edge
cp ../genray.in_CMod_LH_edge genray.in
cp ../g1060728011.01100 .
../../xgenray
ncdump gold-genray.nc > gold.txt
ncdum genray.nc > run.txt
tail -n +2 gold.txt > gold0.txt
tail -n +2 run.txt > run0.txt
diff gold0.txt run0.txt
```

