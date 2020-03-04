#!/bin/bash
../../xgenray &> logfile
ncdump gold-genray.nc > gold.txt
ncdump genray.nc > run.txt
tail -n +2 gold.txt > gold0.txt
tail -n +2 run.txt > run0.txt
diff gold0.txt run0.txt
