This directory provies an example starting point for creating cycling 3dvar experiments. The files provided will generate an SST analysis over 2 cycles, with the sample observation files provided. During the first cycle, the BUMP correlation operator will be initialized.

## To test
The following assumes that you have already checked out the repository and built it in `generic-marine-jedi/build`. 

Create a directory for your experiment and copy the default cycling configuration and example observations. Within the root `generic-marine-jedi` directory run: 

```
mkdir exp
cd exp
cp ../cycling/* . -r
```

Run the experiment

```
./cycle.sh
```

you should see several directories (`ana`, `inc`,`obs_out`) populated with the output from the 3dvar.

## To use for real
To use for your own experiments, we'll obviously want to obtain your own observations, and modify the `config/*` and `cycle.sh` files accordingly.