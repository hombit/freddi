#!/bin/bash
# Run from python/test/data/ directory
# FREDDI=../../../cmake-build-debug/freddi ./generate_test_data.sh

FREDDI="${FREDDI-freddi} --precision=6"

$FREDDI --prefix=default_args

$FREDDI --prefix=sinusF_Thot --Thot=1e4 --initialcond=sinusF

$FREDDI --prefix=gaussF_Thot --Thot=1e4 --initialcond=gaussF

$FREDDI --prefix=quasistat_Mdisk0 --initialcond=quasistat --Mdisk0=1e26

$FREDDI --prefix=quasistat_Thot_Cirr_Tirr --Cirr=2e-4 --angulardistdisk=isotropic --boundcond=Tirr --Thot=1e4 --initialcond=quasistat

$FREDDI --prefix=quasistat_Thot_Cirr_Tirr_Qirr2Qvis --Cirr=2e-4 --angulardistdisk=isotropic --boundcond=Tirr --Thot=1e4 --initialcond=quasistat --Qirr2Qvishot=1.0 --time=100

$FREDDI --prefix=lambdas --lambda=8000 --lambda=5000 --lambda=3000

# Passbands files are from
# http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?mode=browse&gname=Swift&asttype=
$FREDDI --prefix=passbands --time=10 --tau=1 --passband=passbands/Swift_B.dat --passband=passbands/Swift_V.dat
