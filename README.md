# APSC_challenge2

919812 - 10510718

Implementation of the sparse Matrix

Representations:

- COO: coordinates uncompressed, mapping indices to values

- CSR: compressed sparse row

- CSC: compressed sparse column


CHECKLIST:

- [ ] README
- [ ] comments
- [ ] dynamic constructor
- [ ] `compress()`
- [ ] `uncompress()`
- [ ] `is_compressed()`
- [ ] call `operator[]`
- [ ] friend matrix-vector `operator*`
- [ ] extension to matrix-vector product with `Matrix` of just one column
- [ ] matrix-matrix `operator*`
- [x] `norm()`
- [ ] extension of `norm()` working with complex
- [ ] reader of matrix market format
- [ ] test code with `Chrono` utility