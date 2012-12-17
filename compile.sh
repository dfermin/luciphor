#!/bin/sh

PWIZ_ROOT=$PWD/pwiz_root
PWIZ_BJAM=$PWIZ_ROOT/libraries/boost-build/engine/bin/bjam

BOOST_BUILD_PATH=$PWIZ_ROOT/libraries/boost-build 
export BOOST_BUILD_PATH

toolset="gcc"
if [ `uname -s` == "Darwin" ]; then
	toolset="darwin"
fi

cmd="$PWIZ_BJAM toolset=$toolset $@"

#touch src/luciphor.cpp 

clear

echo -e "\n$cmd\n";

$cmd 

GCCVERSION=$(gcc -dumpversion);

if [ $toolset == "gcc" ]; then
	mv ./bin/gcc-$GCCVERSION/release/link-static/threading-multi/luciphor .
else
	mv ./bin/darwin-*/release/link-static/threading-multi/luciphor .
fi


