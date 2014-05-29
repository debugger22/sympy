#! /usr/bin/env bash

# Exit on error
set -e
# Echo each command
set -x

if [[ "${TEST_SPHINX}" == "true" ]]; then
    cd doc
    make html-errors
    make clean
    make latex
    cd _build/latex
    export LATEXOPTIONS="-interaction=nonstopmode"
    make all
elif [[ "${TEST_SAGE}" == "true" ]]; then
    sage -v
    sage -python bin/test sympy/external/tests/test_sage.py
else
    # We change directories to make sure that we test the installed version of
    # sympy.
    mkdir empty
    cd empty

    if [[ "${TEST_DOCTESTS}" == "true" ]]; then
        cat << EOF | python
import sympy
if not sympy.doctest():
    raise Exception('Tests failed')
EOF
        cd ..
        bin/doctest doc/
    elif [[ "${TEST_SLOW}" == "true" ]]; then
        cat << EOF | python
# Pep8 tests
pep8 --ignore=E501,E502,E121,E122,E124,E123,E125,E126,E127,E128,E129,E131,E262,E265,E226 ./sympy/
import sympy
if not sympy.test(split='${SPLIT}', slow=True, timeout=180):
    # Travis times out if no activity is seen for 10 minutes. It also times
    # out if the whole tests run for more than 50 minutes.
    raise Exception('Tests failed')
EOF
    elif [[ "${TEST_THEANO}" == "true" ]]; then
        cat << EOF | python
import sympy
if not sympy.test('*theano*'):
    raise Exception('Tests failed')
EOF
    else
        cat << EOF | python
import sympy
if not sympy.test(split='${SPLIT}'):
    raise Exception('Tests failed')
EOF
        fi
fi
