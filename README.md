# QTE (Quantum Time Evolution)

This repository contains a set of Max/MSP externals for exploring quantum time evolution.

## Package layout

```
README.md
LICENSE
package-info.json
source/      # C sources and CMake project
externals/   # compiled externals (not tracked)
help/        # help patchers
```

## Building

The externals are built with CMake against the Max SDK.

```
cd source
mkdir build && cd build
cmake ..
make
```

The resulting `.mxo` bundles will be found in the `../externals` folder.

### Using `qte.timedev`

Instantiate `[qte.timedev n tsteps tmin tmax]`. Send a list of `n` eigenvalues
followed by `2*n` floats representing the complex coefficients `(re, im)`. A
`bang` then outputs a list of magnitude/phase pairs for each time step.

