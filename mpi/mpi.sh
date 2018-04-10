#!/bin/sh
#cmm comments out all write(*,*) and pause statements, for MPI runs
#auto makes MPI insertions in the *.f files at CMPIINSERTPOSITION points,
#  insertions specified in mpi.ins.

mkdir TMP_GENRAY
cp *.f *.i *.inc mpi/mpi.ins TMP_GENRAY
cd TMP_GENRAY
echo Start preprocessing
../mpi/cmm -a *.f
../mpi/cmm -p WRITE "write(*," *.f
../mpi/cmm -p PAUSE "pause" *.f
echo auto
../mpi/auto genray.f mpi.ins >genray_.f
mv genray_.f genray.f
echo Start building
$1 -c *.f
$1 -o $2 *.o $3 $4 
cd ..
cp TMP_GENRAY/$2 .
rm -rf TMP_GENRAY
