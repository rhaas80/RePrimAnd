This is a simple con2prim code based on the RePrimAnd library.

Build instructions for RePrimAnd are available on its documentation page:

https://wokast.github.io/RePrimAnd/installing.html#building-from-source

Which requires the `meson` build system plus a number of extra libraries
(`GSL`, `Boost`).

`con2prim/minimal.cc` contains a demo code that takes an input file with 3 data
items per line (dens, tau, S) and computes the primite variable "pressure" from
them and stores the results in a file `output.txt`.
