# Description

Luca Brambilla: luca13.brambilla@mail.polimi.it

919812 - 10510718

Challenge source code available [here](https://github.com/luca-brambilla/APSC_challenge2)

# Sparse matrix

Implementation of the sparse Matrix

Representations:

- COO: coordinates uncompressed, mapping indices to values

- CSR: compressed sparse row

- CSC: compressed sparse column

Documentation for the template class `Matrix` available 
[here](https://luca-brambilla.github.io/APSC_challenge2/classalgebra_1_1Matrix.html)


# Installation

Download the repository using an SSH key:

```sh
git clone git@github.com:luca-brambilla/APSC_challenge2.git
```

run the code using:

```sh
make
```

# Additional instructions

If you want to read a full matrix you can set a threshold for considering a number as zero, thus not adding it as an element of the matrix.

To do so you can edit the `Makefile` adding a new flag to the compiler:
```
CPPFLAGS += -D ZERO_TOL=1e-8
```

If the variable `ZERO_TOL` is not redefined, the default value in the header source code is `1e-08`.

# Testing

In the `main.cpp` file is possible to test different sections of the cody by changing the values from `false` to `true` in the `if` statements.