#!/usr/bin/env bash

CALL_DIR=$PWD

BASEDIR=$(dirname "$0")
cd $BASEDIR

# Submodule downloads
git submodule init
git submodule update --recursive

# Prepare SDSL for compilation
cd lib/sdsl
sed -i "s/cmake_minimum_required.*/cmake_minimum_required(VERSION 3.16)/g" CMakeLists.txt
git submodule init
git submodule update
sed -i "s/cmake_minimum_required.*/cmake_minimum_required(VERSION 3.16)/g" external/libdivsufsort/CMakeLists.txt
sed -i "s/cmake_minimum_required.*/cmake_minimum_required(VERSION 3.16)/g" external/googletest/CMakeLists.txt
sed -i "s/cmake_minimum_required.*/cmake_minimum_required(VERSION 3.16)/g" external/googletest/googletest/CMakeLists.txt
cd -

# Prepare KMC for compilation
cp lib/CMakeLists_kmc.txt lib/KMC/kmc_api/CMakeLists.txt

# cmph 2.0.2 dowload
if [ ! -d lib/cmph-2.0.2 ]; then
	curl -O https://altushost-swe.dl.sourceforge.net/project/cmph/v2.0.2/cmph-2.0.2.tar.gz && mv cmph-2.0.2.tar.gz lib/
	tar -xzf lib/cmph-2.0.2.tar.gz -C lib/ && rm lib/cmph-2.0.2.tar.gz
	# Prepare cmph partial compilation
	cp lib/CMakeLists_cmph.txt lib/cmph-2.0.2/src/CMakeLists.txt
fi



# cd lib
# git clone --recursive https://github.com/medvedevgroup/ESSCompress.git
# cd ../
# cd lib/ESSCompress
# ./INSTALL
# cd ../../

# git clone --recursive https://github.com/jermp/sshash.git

cp lib/parse_file.hpp lib/sshash/include/builder/
# Compilation
cmake . && make -j

#export PATH=$BASEDIR/lib/ESSCompress/bin:$BASEDIR/bin:$PATH

#echo "Path of ESS-Compress $(which ESSCompress)"
echo "Path of ESS-Color $(which genMatrix;)"
cd $CALL_DIR
