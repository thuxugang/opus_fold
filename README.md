# OPUS-Fold

Only source codes are uploaded to this repo, the contents in constrain_files and other_potentials are not uploaded due to their sizes. The complete version will comming soon.


## Requirements

For your convenience, the essential packages we used can be found in *requirements* folder.

***gcc/g++ 4.8 is necessary to reproduce our results.***

```
gcc/g++ 4.8
automake
autoconf
libtool
mongodb
mongo-c-driver
mongo-cxx-driver
eigen3
boost
```

## Installation

To help you deploy OPUF-Fold in your environment, we show the docker image making piplines from ubuntu:18.04 docker image as follows: 

1. Download OPUS-Fold and its databases.

   comming soon.

2. Start Docker image and copy the files of OPUS-Fold into it.
   ```
   docker pull ubuntu:18.04
   docker run -it --name fold ubuntu:18.04
   docker start fold
   docker attach fold 
   docker cp ./opus-fold/ fold:/home/opus_fold
   ```
3. Install automake, autoconf and libtool
   ```
   apt-get update
   apt-get install cmake
   apt-get install -y build-essential
   apt-get install git
   apt-get install zip
   apt-get install automake autoconf libtool
   ```
4. Use gcc/g++ 4.8
   ```
   apt-get install -y gcc-4.8
   apt-get install -y g++-4.8
   cd /usr/bin
   rm gcc
   ln -s gcc-4.8 gcc
   rm g++
   ln -s g++-4.8 g++
   ```
5. Install MongoDB 
   ```
   tar xzf mongodb-linux-x86_64-3.0.6.tgz
   mv mongodb-linux-x86_64-3.0.6 mongodb
   cd mongodb
   mkdir data
   mkdir log
   ```
6. Install MongoDB C Driver
   ```
   tar xzf mongo-c-driver-1.13.0.tar.gz
   cd mongo-c-driver-1.13.0
   mkdir cmake-build
   cd cmake-build
   cmake -DENABLE_AUTOMATIC_INIT_AND_CLEANUP=OFF ..
   make
   make install
   ```
7. Install mongocxx Driver
   ```
   tar xzf r3.3.1.tar.gz
   cd mongo-cxx-driver-r3.3.1/build
   cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local ..
   make EP_mnmlstc_core
   make
   make install
   ```
8. Install eigen
   ```
   tar xzf 3.3.7.tar.gz
   cd eigen-eigen-323c052e1731
   mkdir build
   cd build
   cmake ..
   make install
   ```
9. Install Boost
   ```
   apt-get install libboost-all-dev
   ```
10. Import and index OPUS-TA/OPUS-CSF/OPUS-DASF/RotamerLib into MongoDB
   ```
   cd /home/mongodb/bin
   ./mongod -dbpath /home/mongodb/data -logpath /home/mongodb/log/mongodb.log -logappend -fork -port 27017

   ./mongoimport -d csf_db -c csf_5 --file /home/opus_fold/opus_fold_db/csf_cpp_5.dat --type json
   ./mongoimport -d csf_db -c csf_7 --file /home/opus_fold/opus_fold_db/csf_cpp_7.dat --type json
   ./mongoimport -d csf_db -c csf_9 --file /home/opus_fold/opus_fold_db/csf_cpp_9.dat --type json
   ./mongoimport -d csf_db -c csf_11 --file /home/opus_fold/opus_fold_db/csf_cpp_11.dat --type json

   ./mongoimport -d ta_db -c ta_3 --file /home/opus_fold/opus_fold_db/ta_cpp_3.dat --type json
   ./mongoimport -d ta_db -c ta_5 --file /home/opus_fold/opus_fold_db/ta_cpp_5.dat --type json
   ./mongoimport -d ta_db -c ta_7 --file /home/opus_fold/opus_fold_db/ta_cpp_7.dat --type json

   ./mongoimport -d dasf_db -c dasf_5 --file /home/opus_fold/opus_fold_db/dasf_cpp_5.dat --type json
   ./mongoimport -d dasf_db -c dasf_7 --file /home/opus_fold/opus_fold_db/dasf_cpp_7.dat --type json
   ./mongoimport -d dasf_db -c dasf_9 --file /home/opus_fold/opus_fold_db/dasf_cpp_9.dat --type json
   ./mongoimport -d dasf_db -c dasf_11 --file /home/opus_fold/opus_fold_db/dasf_cpp_11.dat --type json
   ./mongoimport -d dasf_db -c rotamer --file /home/opus_fold/opus_fold_db/dasf_rotamer.dat --type json

   ./mongo
   use csf_db
   db.csf_5.ensureIndex({"key":1},{"unique":true})
   db.csf_7.ensureIndex({"key":1},{"unique":true})
   db.csf_9.ensureIndex({"key":1},{"unique":true})
   db.csf_11.ensureIndex({"key":1},{"unique":true})

   use ta_db
   db.ta_3.ensureIndex({"key":1},{"unique":true})
   db.ta_5.ensureIndex({"key":1},{"unique":true})
   db.ta_7.ensureIndex({"key":1},{"unique":true})

   use dasf_db
   db.dasf_5.ensureIndex({"key":1},{"unique":true})
   db.dasf_7.ensureIndex({"key":1},{"unique":true})
   db.dasf_9.ensureIndex({"key":1},{"unique":true})
   db.dasf_11.ensureIndex({"key":1},{"unique":true})
   db.rotamer.ensureIndex({"key":1},{"unique":true})
   ```  
11. Build OPUS_Fold
   ```
   cd opus_fold
   sh ./compile.sh
   ```   
12. Run OPUS_Fold

   Put the constrained contact map files and their corresponding list into `constains_files/contact_maps`, constrained torsional angles files and their corresponding list into `constains_files/torsion_angles` and initial structures files and their corresponding list into `constains_files/init_structures`.***The data in each line of the input file should be separated by a single ' ' or '\t'. '#' should be added at the beginning of each comment line.***
   
   Set the parameters in *opus_fold.ini*. 
   ```
   ./opus_fold
   ```


