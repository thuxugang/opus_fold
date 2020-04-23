#!/bin/bash

echo "=====Compiling====="

if [ ! -d "src/build" ]; then
mkdir src/build
else
rm -r src/build
mkdir src/build
fi

cp src/algorithm/* src/build/
cp src/buildprotein/* src/build/
cp src/database/* src/build/
cp src/myclass/* src/build/
cp src/myio/* src/build/
cp src/potential/* src/build/
cp src/opus_fold.cpp src/build/

cd src/build
cmake ..
make
cp ./opus_fold ../../
cp ../opus_fold.ini ../../

echo "=====Cleaning====="
cd ..
rm -r ./build

echo "=====Done====="

