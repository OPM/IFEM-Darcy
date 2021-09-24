#!/bin/bash

function clone_ifem {
  # Clone IFEM
  if ! test -d ${WORKSPACE}/deps/IFEM
  then
    pushd .
    mkdir -p $WORKSPACE/deps/IFEM
    cd $WORKSPACE/deps/IFEM
    git init .
    git remote add origin https://github.com/OPM/IFEM
    git fetch --depth 1 origin $IFEM_REVISION:branch_to_build
    test $? -eq 0 || exit 1
    git checkout branch_to_build
    popd
  fi
}

declare -a sidestreams
sidestreams=(IFEM-AdvectionDiffusion)

declare -A sidestreamRev
sidestreamRev[IFEM-AdvectionDiffusion]=master

declare -a downstreams
downstreams=(IFEM-CoSTA)

declare -A downstreamRev
downstreamRev[IFEM-CoSTA]=master

IFEM_REVISION=master
if grep -qi "ifem=" <<< $ghprbCommentBody
then
  IFEM_REVISION=pull/`echo $ghprbCommentBody | sed -r 's/.*ifem=([0-9]+).*/\1/g'`/merge
fi

clone_ifem

source $WORKSPACE/deps/IFEM/jenkins/build-ifem-module.sh

parseRevisions
printHeader IFEM-Darcy

build_module_and_upstreams IFEM-Darcy

test $? -eq 0 || exit 1

# If no downstream builds we are done
if ! grep -q "with downstreams" <<< $ghprbCommentBody
then
  exit 0
fi

clone_sidestreams Darcy
build_downstreams IFEM-Darcy

test $? -eq 0 || exit 1
