#!/bin/bash
# Run from test/data/ directory
# FREDDI=../../../cmake-build-debug/freddi ./generate_test_data.sh

FREDDI=${FREDDI-freddi}

$FREDDI --prefix=default_args

$FREDDI --prefix=sinusF_Thot --Thot=1e4 --initialcond=sinusF

$FREDDI --prefix=gaussF_Thot --Thot=1e4 --initialcond=gaussF

$FREDDI --prefix=quasistat_Mdisk0 --initialcond=quasistat --Mdisk0=1e26

$FREDDI --prefix=quasistat_Thot_Cirr_Tirr --Cirr=1e-3 --boundcond=Tirr --initialcond=quasistat

$FREDDI --prefix=lambdas --lambda=8000 --lambda=5000 --lambda=3000
