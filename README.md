# fl0p : C++ Interpolated flow

[![standard-readme compliant](https://img.shields.io/badge/readme%20style-standard-brightgreen.svg?style=flat-square)](https://github.com/RichardLitt/standard-readme)

A C++ library that provides flow field interpolation of stored data.

This C++ library provides a continuous flow implementation that interpolates fixed position data.
Simple examples are provided and should be enough to explain how to use this library.

This repository contains:

1. The software itself provided as a header only library in the directory [include/sl0](./include/sl0)
2. A few [examples](./examples).

## Table of Contents

- [Background](#background)
- [Install](#install)
- [License](#license)

## Background

This library has been produced during my PhD thesis and as part as the European Research Council project: [C0PEP0D](https://c0pep0d.github.io/)
This library is used as part of [SHELD0N](https://github.com/C0PEP0D/sheld0n), a lagrangian particle advection software.

## Install

### Dependencies

* [**CMake** `v?`](https://cmake.org/download/) or higher must be installed
* a c++17 compliant compiler, such as [**gcc** `v9`](https://gcc.gnu.org/) or higher must be installed. For Ubuntu users: [ppa](https://launchpad.net/%7Ejonathonf/+archive/ubuntu/gcc?field.series_filter=bionic).
* the [**Threading Building Block Library** `v2018`](https://github.com/ibaned/tbb) or higher must be installed ([this version](https://github.com/wjakob/tbb) that enables is installing it using CMake is advised)

Examples:
* [**Eigen**](https://eigen.tuxfamily.org) must be installed
* [**fl0w**](https://github.com/C0PEP0D/fl0w) must be installed
* [**m0sh**](https://github.com/C0PEP0D/m0sh) must be installed

Vtk example:
* [**v0l**](https://github.com/C0PEP0D/v0l) must be installed

The examples assume the following directory tree structure:
```bash
..
 ├── .
 │   │── fl0p
 │   │── fl0w
 │   │── m0sh
 │   └── (v0l)
 └── thirdparty
     └── eigen
```
One should either install these dependencies accordingly, or adapt their path in the **CMakeList.txt** file of the examples.

### Installing

Start by cloning this repository.

```sh
$ git clone https://gitlab.com/rmonthil/c0pep0d.git
```

### Examples

Running an example:

```bash
$ cd examples/stationary
$ mkdir build
$ cd build
$ cmake ..
$ make
$ ./stationary
```

### Updating

A simple pull should be enough.

```sh
$ git pull
```

## Maintainers

Rémi Monthiller - [@rmonthil](https://gitlab.com/rmonthil) - remi.monthiller@gmail.com

## Contributing

Feel free to dive in! [Open an issue](https://github.com/rmonthil/c0pep0d/issues/new) or submit PRs.

## License

[MIT © Centrale Marseille, Rémi MONTHILLER.](./LICENSE)
