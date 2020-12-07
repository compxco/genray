#!/bin/bash
# July 17, 2007
# This script will compare a results of present test runs with
# saved results in this directory.
# You can execute this script, but it will run all the tests and
# give a jumble of overlapping plots.  
# To see individual test comparisons, scrape in to a terminal the defn 
# of XGENRAY, then scrape in commands for each test section.

# Edit following line to point to genray executable to be tested.
#XGENRAY='/home/bobh/cql3d/genray/genray_cvs/code/gr_v6_070720/xgenray'
#XGENRAY='/home/bobh/cql3d/genray/genray_cvs/code/genray_v6-0/xgenray'
#XGENRAY='/home/bobh/cql3d/genray/genray_cvs/code/gr070121/xgenray'
#XGENRAY='/home/bobh/cql3d/genray/genray_cvs/code/gr061106/xgenray'
#XGENRAY='/home/bobh/cql3d/genray/genray_cvs/code/gr060211/xgenray'
#XGENRAY='/home/bobh/cql3d/genray/genray_cvs/code/gr060421/xgenray'
#XGENRAY='/home/bobh/cql3d/genray/genray_cvs/code/gr_v6_070720_25_1/xgenray'
#XGENRAY='/home/bobh/cql3d/genray/genray_cvs/code/genray_v6-1_070726/xgenray'
#XGENRAY='/home/bobh/cql3d/genray/genray_cvs/code/genray_v6-3_070801/xgenray'
#XGENRAY='/home/bobh/cql3d/genray/genray_cvs/code/genray_v7-8_080801/xgenray'
#XGENRAY='/home/bobh/cql3d/genray/genray_cvs/code/genray_v7-11_081007/xgenray'
#XGENRAY='/home/bobh/cql3d/genray/genray_cvs/code/gr081202/xgenray'
#XGENRAY='/home/bobh/cql3d/genray/genray_cvs/code/genray_v7-17_090527/xgenray_nowrite'
#XGENRAY='/home/bobh/cql3d/genray/genray_cvs/code/gr090901/xgenray_nowrite'
#XGENRAY='/home/bobh/cql3d/genray/genray_cvs/code/gr090903_090921/xgenray_debug_BOUND_CHECK'
#XGENRAY='/home/bobh/cql3d/genray/genray_cvs/code/genray_cswim_svn/trunk/xgenray'
#XGENRAY='/home/bobh/cql3d/genray/genray_cvs/code/genray_cswim_svn/trunk/xgenray_nowrite'
#XGENRAY='/home/bobh/cql3d/genray/genray_cvs/code/gr110314/code/xgenray'
#XGENRAY='/home/bobh/cql3d/genray/genray_cvs/code/gr110503/xgenray'
#XGENRAY='/home/bobh/cql3d/genray/genray_cvs/code/genray_v9_0_110801/xgenray'
#XGENRAY='/home/bobh/cql3d/genray/genray_cvs/code/genray_v10.0_111031/xgenray'
#XGENRAY='/home/bobh/cql3d/genray/genray_cvs/code/gr111031_rz_launch_rz_dens_alloc_zeff/xgenray'
#XGENRAY='/home/bobh/cql3d/genray/genray_cvs/code/genray_cswim_svn/trunk/xgenray'
#XGENRAY='/home/bobh/cql3d/genray/genray_cvs/code/genray_cswim_svn/trunk/xgenray_nowrite'
#XGENRAY='/home/bobh/cql3d/genray/genray_cvs/code/genray_cswim_svn/Bertelli_genray_last_version_16July2014_mom_cons/xgenray'
#XGENRAY='/home/bobh/cql3d/genray/genray_cvs/code/genray_cswim_svn/Bertelli_genray_last_version_16July2014_mom_cons/xgenray_nowrite'
#XGENRAY='../../xgenray_intel.edison'   #At NERSC, edison
#XGENRAY='/home/bobh/cql3d/genray/genray_cvs/code/genray_cswim_svn/trunk/xgenray'
#XGENRAY='/home/bobh/cql3d/genray/genray_cvs/code/genray_cswim_svn/genray_170301/xgenray'
#XGENRAY='/home/bobh/cql3d/genray/genray_cvs/code/genray_cswim_svn/genray_170301_updtd170725/xgenray'
#XGENRAY='/home/bobh/cql3d/genray/genray_cvs/code/genray_cswim_svn/genray_v10.11_170725/xgenray'
#XGENRAY='../$PWD/xgenray'
#XGENRAY='/global/homes/u/u650/genray_cswim_svn/genray_v10.11_170725/xgenray_intel.edison'   #At NERSC, edison.

#
#define GENRAY executable location
#
XGENRAY='../../myxgenray'

#
#define JXDRAW="yes" to execute XDRAW else "no"
#define XDRAW executable location 
#
JXDRAW="no"
XDRAW1="../xdraw/xdraw"
XDRAW="../$XDRAW1"


#
#define execution of TEST_XX
#TEST_XX="do" to execute or TEST_XX="skip" to skip it
#
TEST_1="do"
TEST_2="do"
TEST_3="do"
TEST_4="do"
TEST_4_1="do"
TEST_5="do"
TEST_5_1="do"
TEST_6="do"
TEST_7="do"
TEST_7_1="do"
TEST_7_EDGE_DENS_RZ="do"
TEST_7_MPI="skip"
TEST_8="do"
TEST_8_1="do"
TEST_8_2="do"
TEST_8_3="do"
TEST_9="do"
TEST_10="do"
TEST_10_1="skip"

if [ "$TEST_1" = "do" ]; then
mkdir test1
cp genray.dat_SS_Xmode_180_emission test1/genray.dat
cp Scen4_bn2.57_129x129 test1/
cp *.in xdraw.ini idl_emission_spectrum_colorbar.pro test1/
cd test1
time $XGENRAY > log_test1   #1m30s on bob8  (1m48s, compx1 w bounds-check)
                            #23s on CompX2
rm  con1 eps* genray.doc gr3d* *.sap zrn.dat
$XDRAW genr &
cd ..
cp genray.bin_test1 genray.bin
$XDRAW1 genr &
cd test1
$XDRAW em &
cd ..
cp emis.bin_test1 emis.bin
$XDRAW1 em &
fi

#gs xdraw_genr_SS_Xmode_180_emission.ps&
#gs xdraw_em_SS_Xmode_180_emission.ps&
#echo "test1: Compare results with xdraw_genr_SS_Xmode_180_emission.ps and xdraw_em_SS_Xmode_180_emission.ps"
#sleep 10

if [ "$TEST_2" = "do" ]; then
mkdir test2
cp genray.dat_HPRT_test_case_061219 test2/genray.dat
cp equilib.dat_HPRT_test_case_061219 drawgenr.in test2/
cd test2
time $XGENRAY > log_test2    #3.04s on compx1
if [ -f "con1" ]; then
rm  con1 eps* genray.doc gr3d* *.sap zrn.dat
fi
if [ "$JXDRAW" = "yes" ]; then
$XDRAW genr &
fi
cd ..
cp genray.bin_test2 genray.bin
if [ "$JXDRAW" = "yes" ]; then
$XDRAW1 genr &
#echo "test2: Compare results with idl_105830p193r1.bob_GR030417.4_gr061120.ps"
fi
fi

if [ "$TEST_3" = "do" ]; then
mkdir test3
cp genray.dat_113544.00325_mod test3/genray.dat
cp g113544.00325_mod drawgenr.in drawem.in test3/
cd test3/
time $XGENRAY > log_test3    #2m59s on compx1 (3m23s w bounds-check) 
                             #30.2s on CompX2
if [ -f "con1" ]; then
rm  con1 eps* genray.doc gr3d* *.sap zrn.dat
fi
if [ "$JXDRAW" = "yes" ]; then
$XDRAW genr &
fi
cd ..
cp genray.bin_test3 genray.bin
if [ "$JXDRAW" = "yes" ]; then
$XDRAW1 genr &
fi
cd test3
if [ "$JXDRAW" = "yes" ]; then
$XDRAW em &
fi
cd ..
cp emis.bin_test3 emis.bin
if [ "$JXDRAW" = "yes" ]; then
$XDRAW1 em &
fi
#cd ..
#echo "test3: Compare results with emis.bin_113544.00325_mod and genray.bin_113544.00325_mod"
fi

if [ "$TEST_4" = "do" ]; then
mkdir test4
cp genray.dat_GT_PoP2005_i_ox=1.1_updtd test4/genray.dat
cp g113544.00325_mod drawgenr.in test4/
cd test4/
time $XGENRAY > log_test4    #13.2s on compx1, 4.53s on CompX2(2016)
emacs ECcone_optimal.dat &
cd ..
emacs ECcone_optimal.dat &
echo "test4: Compare output text files ECcone_optimal.dat"
fi

if [ "$TEST_4_1" = "do" ]; then
mkdir test4.1
cp genray.dat_GT_PoP2005_id6_i_ox2_090720 test4.1/genray.dat
cp g113544.00325_mod drawgenr.in drawfreqelec.in test4.1/
cd test4.1/
time $XGENRAY > log_test4.1    #27s on compx1
if [ "$JXDRAW" = "yes" ]; then
$XDRAW genr &
fi
cd ..
echo "test4.1: Compare genray.bin_test4.1 with test4/genray.bin using xdraw"
cp genray.bin_test4.1 genray.bin
if [ "$JXDRAW" = "yes" ]; then
$XDRAW1 genr &
fi
fi

if [ "$TEST_5" = "do" ]; then
mkdir test5
cp genray.dat_CANONICAL_2004_ITER_TEST_one_ray  test5/genray.dat
cp  g521022.01000 drawgenr.in test5/
cd test5/
time $XGENRAY > log_test5    #1.6s on compx1
if [ "$JXDRAW" = "yes" ]; then
$XDRAW genr &
fi
cd ..
echo "test5: Compare genray.bin_CANONICAL_2004_ITER_TEST_one_ray with test5/genray.bin using xdraw"
cp genray.bin_test5 genray.bin
if [ "$JXDRAW" = "yes" ]; then
$XDRAW1 genr &
fi
fi

if [ "$TEST_5_1" = "do" ]; then
mkdir test5.1
cp genray.dat_id10iabsorp1  test5.1/genray.dat
cp  g521022.01000 drawgenr.in test5.1/
cd test5.1/
time $XGENRAY > log_test5.1    #1m52s on compx1, 17.8s on CompX2
if [ "$JXDRAW" = "yes" ]; then
$XDRAW genr &
fi
cd ..
echo "test5.1: Compare  genray.bin_id10iabsorp1 with test5.1/genray.bin using xdraw"
cp genray.bin_test5.1 genray.bin
if [ "$JXDRAW" = "yes" ]; then
$XDRAW1 genr &
fi
fi

if [ "$TEST_6" = "do" ]; then
mkdir test6
cp genray.dat_shot106270_multi_ray_multi_cone test6/genray.dat
cp g106270.02500 drawgenr.in test6/
cd test6/
time $XGENRAY > log_test6    #15.7s on compx1, 3.2 s on CompX2
if [ "$JXDRAW" = "yes" ]; then
$XDRAW genr &
fi
cd ..
echo "test6: Compare genray.bin_shot106270 with test6/genray.bin using xdraw"
cp genray.bin_test6 genray.bin
if [ "$JXDRAW" = "yes" ]; then
$XDRAW1 genr &
fi
fi

if [ "$TEST_7" = "do" ]; then
echo "test 7"
mkdir test7
cp genray.in_CMod_LH_edge  test7/genray.in
cp g1060728011.01100 drawgenr.in test7/
cd test7/
time $XGENRAY > log_test7    #1m54s on compx1 w xgenr_nw,  18.6s on CompX2
if [ "$JXDRAW" = "yes" ]; then
$XDRAW genr &
fi
cd ..
echo "test7: Compare genray.bin_CMod_LH_edge with test7/genray.bin using xdraw"
cp genray.bin_test7 genray.bin
if [ "$JXDRAW" = "yes" ]; then
$XDRAW1 genr &
fi
fi

#An MPI test, if at edison.nersc.gov
#if [ "$TEST_7_MPI" = "do" ]; then
#echo "test 7 MPI"
#mkdir test7.mpi
#cp genray.in_CMod_LH_edge  test7.mpi/genray.in
#cp g1060728011.01100 drawgenr.in test7.mpi/
#cp edison_batchscript_mpi.slurm test7.mpi/
#cd test7.mpi/
#sbatch edison_batchscript_mpi.slurm
#Wait for batch job to complete
#if [ "$JXDRAW" = "yes" ]; then
#$XDRAW genr &
#fi
#cd ..
#echo "test7.mpi: Compare genray.bin_CMod_LH_edge with test7/genray.bin using xdraw"
#cp genray.bin_test7 genray.bin
#if [ "$JXDRAW" = "yes" ]; then
#$XDRAW1 genr &
#fi
#fi

if [ "$TEST_7_1" = "do" ]; then
echo "test 7.1"
mkdir test7.1  #Bonoli-Englade LH dispersion with warm plasma term
cp genray.in_CMod_LH_edge_id16  test7.1/genray.in
cp g1060728011.01100 drawgenr.in drawonet1.in test7.1/
cd test7.1/
time $XGENRAY > log_test7.1    #25.8s on CompX2
if [ "$JXDRAW" = "yes" ]; then
$XDRAW genr &
fi
cd ..
echo "test7.1: "
echo "test7.1: Compare genray.bin_CMod_LH_edge_id16 with test7/genray.bin using xdraw"
echo "test7.1: Warm plasma gives substantial broadening of heating and current drive (BH180126)."
echo "test7.1: NEEDS further investigation1"
echo "test7.1: "
cp genray.bin_test7 genray.bin
if [ "$JXDRAW" = "yes" ]; then
$XDRAW1 genr &
fi
cd test7.1/
if [ "$JXDRAW" = "yes" ]; then
$XDRAW onet1 &
fi
cd ..
cp onetwo1.bin_test7.1 onetwo1.bin
if [ "$JXDRAW" = "yes" ]; then
$XDRAW1 onet1 &
fi
cp onetwo1.bin_test7 onetwo1.bin
if [ "$JXDRAW" = "yes" ]; then
$XDRAW1 onet1 &
fi
fi

if [ "$TEST_7_EDGE_DENS_RZ" = "do" ]; then
echo "test7-edge_dens_rz"
mkdir test7_edge_dens_rz
cp genray.in_i_edge_dens_rz_mesh=1 test7_edge_dens_rz
cp genray.in_i_edge_dens_rz_mesh=2 test7_edge_dens_rz
cp drawgenr.in test7_edge_dens_rz
cd test7_edge_dens_rz
ln -s ../g1060728011.01100 g1060728011.01100
cp genray.in_i_edge_dens_rz_mesh=1 genray.in
time $XGENRAY > log_test_i_edge_dens_rz_mesh=1  #210 secs
ln -s ./dens_temp_rz_out.dat dens_temp_rz_in.dat
rm psi_rho_rz_out.dat
gzip log_test_i_edge_dens_rz_mesh=1
cp genray.in_i_edge_dens_rz_mesh=2 genray.in
time $XGENRAY > log_test_i_edge_dens_rz_mesh=1  #140 secs,  12.7s on CompX2
gzip log_test_i_edge_dens_rz_mesh=2
if [ "$JXDRAW" = "yes" ]; then
$XDRAW genr &
fi
cd ..
echo "test7_edge_dens_rz: Comparing with saved result, and with"
echo "test7_edge_dens_rz: result of test7, obtained without rz_grid"
echo "test7_edge_dens_rz: [Should be close to test7]."
cp genray.bin_test7_edge_dens_rz genray.bin
if [ "$JXDRAW" = "yes" ]; then
$XDRAW1 genr &
fi
cp genray.bin_test7 genray.bin
if [ "$JXDRAW" = "yes" ]; then
$XDRAW1 genr &
fi
fi

#test8.n  compares various methods of calculation of ITER central CD,
# for both accuracy and speed.  test8 method (adj) is more general,
# but test8.1 with ieffic=4, mom_cons, is much faster for this first
# harmonic, O-mode EC case, and result is about the same.

if [ "$TEST_8" = "do" ]; then
mkdir test8
cp genray.dat_EC_ITER_Central_CD  test8/genray.dat
cp  equilib.dat_EC_ITER_Central_CD test8/equilib.dat
cp drawonet2.in test8
cp drawgenr.in test8
cd test8/
time $XGENRAY > log_test8    #2 min 53s on CompX2 w xgenray_nowrite,  
                             #144s on CompX2
#                            #4 secs if use already computed adj files.
grep runtime log_test8
grep 'total toroidal current toroidal_cur_total' log_test8  #25.8 kA
if [ "$JXDRAW" = "yes" ]; then
$XDRAW genr &
$XDRAW onet2 &
fi
cd ..
echo "test8: Compare CD with Prater_NF2008_Fig12.pdf and EC_ITER_Central_CD.ps (prior)"
cp onetwo1.bin_EC_ITER_Central_CD onetwo1.bin
if [ "$JXDRAW" = "yes" ]; then
$XDRAW1 onet2 &
fi
acroread Prater_NF2008_Fig12.pdf &
fi

if [ "$TEST_8_1" = "do" ]; then
mkdir test8.1
cp genray.dat_EC_ITER_Central_CD_ieffic4_mom_cons  test8.1/genray.dat
cp  equilib.dat_EC_ITER_Central_CD test8.1/equilib.dat
cp drawonet2.in test8.1
cd test8.1/
time $XGENRAY > log_test8.1    #0.60 secs on CompX2 w xgenray_nowrite
echo 'ieffic=4, ieffic_mom_cons=1'
grep runtime log_test8.1
grep 'total toroidal current toroidal_cur_total' log_test8.1   #25.5 kA
if [ "$JXDRAW" = "yes" ]; then
$XDRAW onet2 &
fi
cd ..
fi

if [ "$TEST_8_2" = "do" ]; then
mkdir test8.2
cp genray.dat_EC_ITER_Central_CD_ieffic4  test8.2/genray.dat
cp  equilib.dat_EC_ITER_Central_CD test8.2/equilib.dat
cp drawonet2.in test8.2
cd test8.2/
time $XGENRAY > log_test8.2    #0.60 secs on CompX2 w xgenray_nowrite
echo 'ieffic=4'
grep runtime log_test8.2
grep 'total toroidal current toroidal_cur_total' log_test8.2   #13.0 kA
if [ "$JXDRAW" = "yes" ]; then
$XDRAW onet2 &
fi
cd ..
fi

if [ "$TEST_8_3" = "do" ]; then
mkdir test8.3
cp genray.dat_EC_ITER_Central_CD_ieffic3  test8.3/genray.dat
cp  equilib.dat_EC_ITER_Central_CD test8.3/equilib.dat
cp drawonet2.in test8.3
cd test8.3/
time $XGENRAY > log_test8.3    #0.59 secs on CompX2 w xgenray_nowrite
echo 'ieffic=3'
grep runtime log_test8.3
grep 'total toroidal current toroidal_cur_total' log_test8.3   #10.1 kA
if [ "$JXDRAW" = "yes" ]; then
$XDRAW onet2 &
fi
cd ..
fi

if [ "$TEST_9" = "do" ]; then
mkdir test9
cp genray.in_CMod_LH_edge_LSC  test9/genray.in
cp g1060728011.01100 drawgenr.in test9/
cd test9/
time $XGENRAY > log_test9    #1m58s on compx1 w xgenr_nw
if [ "$JXDRAW" = "yes" ]; then
$XDRAW genr &   #LSC damping, but on the first interation, with Maxwl distn.
fi
cp ../drawdelpwr_lsc.in .
if [ "$JXDRAW" = "yes" ]; then
$XDRAW delpwr_lsc &  #Gives QL damping with LH plateau
fi
cd ..
echo "test9: Compare genray.bin_CMod_LH_edge_LSC with test7/genray.bin using xdraw"
echo "test9: Here, QL damping on Maxwellian with LSC-like option is very close to linear damping in test7"
echo "test9:  QL damping with LH plateau reduces damping for all rays, in this case"
cp genray.bin_test7 genray.bin    #You can see slight differences from above LSC damping on a Maxwl.
cp delpwr.bin_test9 delpwr.bin
if [ "$JXDRAW" = "yes" ]; then
$XDRAW1 genr &
$XDRAW1 delpwr_lsc &
fi
fi

if [ "$TEST_10_1" = "do" ]; then
mkdir test10.1
cp genray.in_test10.1 test10.1
cp eqdsk_MASTU test10.1/
cd test10.1/
time $XGENRAY > log_test10.1
emacs ECcone_optimal.dat &
cd ..
emacs ECcone_optimal.dat_test10.1 &
echo "test10.1: Compare output text files ECcone_optimal.dat"
fi

if [ "$TEST_10" = "do" ]; then
mkdir test10
cp genray.in_test10 test10/genray.in
cp eqdsk_MASTU drawgenr.in drawfreqelec.in test10/
cd test10/
time $XGENRAY > log_test10    #32s on compx2
if [ "$JXDRAW" = "yes" ]; then
$XDRAW genr &
fi
cd ..
echo "test10: Compare genray.bin_test10 with test4/genray.bin using xdraw"
cp genray.bin_test10 genray.bin
if [ "$JXDRAW" = "yes" ]; then
$XDRAW1 genr &
fi
fi

#rm -rf test1 test2 test3 test4 test4.1 test5 test5.1 test6 test7 test7_edge_dens_r test8 test8.1 test8.2 test8.3 test9 test10 test10.1


echo "After finishing viewing results, can remove all subdirectories with a:"
echo  "rm -rf test1 test2 test3 test4 test4.1 test5 test5.1 test6 test6.1 test7 test7.1 test7_ec test7_edge_dens_rz test7_id16 test8 test8.1 test8.2 test8.3 test9 emis.bin genray.bin"
