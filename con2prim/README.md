This is a simple con2prim code based on the RePrimAnd library.

Build instructions for RePrimAnd are available on its documentation page:

https://wokast.github.io/RePrimAnd/installing.html#building-from-source

Which requires the `meson` build system plus a number of extra libraries
(`GSL`, `Boost`).

`con2prim/minimal.cc` contains a demo code that takes an input file with 3 data
items per line (dens, tau, S, all other fields are ignored) and computes the
primite variable "pressure" from them and stores the results in a file
`output.txt`.

Build instructions for NCSA Delta
=================================

These instructions assume you are logging in to Delta using `ssh` via
`login.delta.ncsa.illinois.edu` and that you have *no* conda environment set up
(i.e. your prompt does not have "(base)" or so in it). If you have conda set up
(ever ran `conda init`) then please edit your `.bashrc` and remove / comment
out the lines that `conda init` added (near the end usually). This also assumes
you only have the default modules loaded so that `module list` returns:

```
ekohaes8: ~$ ssh rhaas@login.delta.ncsa.illinois.edu

NCSA Delta System

Login with NCSA Kerberos + NCSA Duo multi-factor.

[...]

[rhaas@dt-login04 ~]$ module list

Currently Loaded Modules:
  1) gcc/11.4.0      3) cuda/11.8.0         5) slurm-env/0.1
  2) openmpi/4.1.6   4) cue-login-env/1.0   6) default-s11
```

## Building RePrimAnd

1. load a modern gcc `module load gcc/11.4.0`
2. load the gsl module `module load gsl/2.7.1`
3. load the boost module `module load boost/1.83.0`
4. create a working directory `mkdir reprimand` and enter it `cd reprimand`
5. set up a Python virtualenv `python3 -m venv $PWD`
6. activate the env `source bin/activate`
7. install `meson` and `ninja` using `pip install meson ninja`
8. make an installation directory `mkdir reprimand-build`
9. clone a copy of RePrimAnd with my changes in it `git clone -b con2prim https://github.com/rhaas80/RePrimAnd.git` and enter it `cd RePrimAnd`
10. set up the build `meson setup --buildtype=release --prefix=$PWD/../reprimand-build mbuild`
11. build and install `ninja -C mbuild install`

After this you will have built RePrimAnd and installed the library and its include files into `reprimand-build`

## Building the con2prim example

Now let's build the demo code to run con2prim for a SLy equqation of state.

1. enter the con2prim directory `cd con2prim`
2. compile using `make` 

## Running a test

The executable produced is called `minimal` and takes as its single input a
file with `dens` `tau` and `scon` values for which to compute the primitive
values for.

Run it as

```
./minimal data.txt
```

and it should produce output to screen like this

```
(reprimand) [rhaas@dt-login04 con2prim]$ ./minimal data.txt
valid range for rho in cold EOS: [0,0.00323642]
Read 1 elements
took: 2.98661e-05 seconds

```

and write the solution found to `output.txt`:

```
(reprimand) [rhaas@dt-login04 con2prim]$ cat output.txt
0.00113933749999999999 0.000382669128088171991 0.000950243985729741351 0.000911466957661988763 0.0443130663757325613 6.17239504222290823e-05 0.600003560342634801
```

where the values output are `dens` `tau` `scon` `rho` `eps` `press` `vel`.

The C++ code `minimal.cc` shows exactly what is done.

The jupyter notebook `SLy.ipnb` is *only* for illustration and to check eg that
one understands how the piecewise polytropic EOS is constructed.

You can have multiple lines of input in `data.txt` so that you can use it to
process a larger number of values for speed tests.
