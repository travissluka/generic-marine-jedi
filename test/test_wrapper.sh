#!/bin/bash

# simple script to combine into a single ctest the 1) running of executable and
# 2)comparing its answers to reference answers.

# run the executable
echo ""
echo "==============================================================================="
echo "Running test executable"
echo "==============================================================================="
${MPI_CMD} $1 $2 testoutput/${COMPARE_TESTNAME}.log
e=$?
if [[ "$e" -gt 0 ]]; then
    echo -e "Failed to run executable. Error code: $e \n"
    exit $e
fi

# run compare, if needed
if [[ "$SKIP_COMPARE" == "FALSE" ]]; then
    echo ""
    echo "==============================================================================="
    echo "Running compare script"
    echo "==============================================================================="

    $COMPARE_SCRIPT testoutput/${COMPARE_TESTNAME}.log \
                    testref/${COMPARE_TESTNAME}.ref \
                    ${COMPARE_TOL_F} ${COMPARE_TOL_I}
    e=$?
    if [[ "$e" -gt 0 ]]; then
        echo -e "Failed in the COMPARE step. Error code: $e \n"
        exit $e
    fi
    echo -e "PASSED \n"
fi
