#!/bin/bash

INTEGRATION_DIR="$1"
SCRIPT_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

if [[ "$INTEGRATION_DIR"x == "x" ]] ; then
    echo "usage:  $0 <Directory to use for checking>"
    exit 1
fi

mkdir -p $INTEGRATION_DIR

if [ "$(ls -A $INTEGRATION_DIR)" ] ; then
    echo "$INTEGRATION_DIR is not empty!"
    exit 1
fi

cd $INTEGRATION_DIR
INTEGRATION_DIR=$(pwd)

cd $SCRIPT_DIR
cd ../..
HASE_ROOT=$(pwd)


## Prepare the folders

mkdir "${INTEGRATION_DIR}/graybat"
mkdir "${INTEGRATION_DIR}/mpi"
mkdir "${INTEGRATION_DIR}/threaded"

for i in {graybat,mpi,threaded}
do
    cp -a example/matlab_example/* "${INTEGRATION_DIR}/$i/"
    cd "${INTEGRATION_DIR}/$i"
    sed  -ir "s/^parallelMode *=.\+/parallelMode = '${i}'/" laserPumpCladdingExample.m
    cd -
done

cp "${SCRIPT_DIR}/integration_testing.m" ${INTEGRATION_DIR}
cp "${SCRIPT_DIR}/gain_analysis.m"       ${INTEGRATION_DIR}
cp "${SCRIPT_DIR}/gain_line_baseline.txt" ${INTEGRATION_DIR}

# create the submit-script

cd $INTEGRATION_DIR
cat "${SCRIPT_DIR}/submit_base.sh"                        > submit.sh
echo "cd ${INTEGRATION_DIR}"                             >> submit.sh
echo "matlab -r integration_testing"                     >> submit.sh
