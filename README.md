# APSC_challenge2

919812 - 10510718

Implementation of the sparse Matrix

Representations:

- COO: coordinates uncompressed, mapping indices to values

- CSR: compressed sparse row

- CSC: compressed sparse column

Documentation available at:
https://luca-brambilla.github.io/APSC_challenge2/

CHECKLIST:

- [ ] README
- [ ] comments
- [ ] dynamic constructor
- [ ] constructor taking dimensions
- [ ] resize uncompressed
- [x] `compress()`
- [x] `uncompress()`
- [ ] `is_compressed()`
- [ ] call `operator[]`
- [ ] friend matrix-vector `operator*`
- [ ] extension to matrix-vector product with `Matrix` of just one column
- [ ] matrix-matrix `operator*`
- [x] `norm()`
- [ ] extension of `norm()` working with complex
- [x] reader of matrix market format
- [ ] test code with `Chrono` utility
- [ ] working for both compressed and uncompressed
- [ ] everything working for CSR
- [ ] everything working for CSC