[![Build Status](https://travis-ci.org/UMD-AOSC/UMD-SST.svg?branch=develop)](https://travis-ci.org/UMD-AOSC/UMD-SST)

# UMD-SST

Testbed sea surface temperature analysis system.

A variational data assmilation system built on the JCSDA Joint Effort for Data Assimilation Integration (JEDI). This project is just beginning, and further documentation will be created as the project progresses.

## Getting Started
The following steps assume you already have one of the JEDI containers setup and working on your computer. See the [JEDI Documentation](https://jointcenterforsatellitedataassimilation-jedi-docs.readthedocs-hosted.com/en/latest/) for additional details on how to do this.

To build UMD-SST and run the continuous integration tests:
```
git clone https://github.com/UMD-AOSC/UMD-SST.git umdsst.src
mkdir -p umdsst.build
cd umdsst.build
ecbuild ../umdsst.src/bundle
cd umdsst
make -j <nprocs>
ctest
```
