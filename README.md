# Generic Marine JEDI

A generic, marine centric, model-less interface to the JCSDA Joint Effort for Data Assimilation Integration (JEDI).

## Getting Started
The following steps assume you already have one of the JEDI containers setup and working on your computer. See the [JEDI Documentation](https://jointcenterforsatellitedataassimilation-jedi-docs.readthedocs-hosted.com/en/latest/) for additional details on how to do this.

### Building / testing

To build Generic Marine JEDI and run the continuous integration tests:
```
git clone https://github.com/travissluka/generic-marine-jedi.git
cd generic-marine-jedi
mkdir build; cd build
ecbuild ../bundle
cd genericMarine
make -j <nprocs>
ctest
```

### Using

TODO: add word here!