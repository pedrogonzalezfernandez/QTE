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

