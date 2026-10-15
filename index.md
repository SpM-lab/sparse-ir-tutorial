## Welcome to sparse-ir tutorial

<!--
This online document provides comprehensive tutorials of sparse-ir libraries in various languages:

- Python: [sparse-ir](https://github.com/SpM-lab/sparse-ir)
- Julia: [SparseIR.jl](https://github.com/SpM-lab/SparseIR.jl)
- C & Fortran: [libsparseir](https://github.com/SpM-lab/libsparseir)

The original `libsparseir` is written in C++ and provides C-API, which is robust and accessible from any programming languages.
The Python, Julia and Fortran libraries are wrappers of `libsparseir` and provide user-friendly interfaces.
We recommend to start with Python or Julia libraries, which are easier to install and use.

The version 2 of `sparse-ir` and `SparseIR.jl` have changed the internal implementation to the C++ implementation of `libsparseir`.
Although most of the interfaces remain unchanged, there are some substantial changes, which **break the backward compatibility**.
Please refer to [the migration guide](src/migration.md) for details.
-->

This online document provides comprehensive tutorials of the Python library [sparse-ir](https://github.com/SpM-lab/sparse-ir) library and Julia library [SparseIR.jl](https://github.com/SpM-lab/SparseIR.jl).

**This tutorial is part of the [sparse-ir ecosystem](https://spm-lab.github.io/sparse-ir-doc/).** Start there for an overview of the libraries (including Rust, the C API and the Fortran bindings) and for the [history and comparison of IR, DLR and MiniPole](https://spm-lab.github.io/sparse-ir-doc/src/history_comparison.html). The Rust versions of these tutorials are in the [Rust guide](https://spm-lab.github.io/sparse-ir-rs/).

The tutorials follow the notation defined on the [notation and conventions page](https://spm-lab.github.io/sparse-ir-doc/src/notation.html) of the sparse-ir documentation.
The basic theory is explained there as well: the [intermediate representation](https://spm-lab.github.io/sparse-ir-doc/src/IR_py.html), [sparse sampling](https://spm-lab.github.io/sparse-ir-doc/src/sparse_sampling_py.html), the [discrete Lehmann representation](https://spm-lab.github.io/sparse-ir-doc/src/DLR.html), and [summation over the Matsubara axis](https://spm-lab.github.io/sparse-ir-doc/src/matsubarasum.html).

Please free to make a contact with the authors of each tutorial page if you have any questions!

## Contributors
* Satoshi Terasaki (AtelierArith)
* Markus Wallerberger (TU Wien)
* Samuel Badr (TU Wien)
* Shintaro Hoshino (Saitama Univ.)
* Fumiya Kakizawa (Saitama Univ.)
* Takashi Koretsune (Tohoku Univ.)
* Yuki Nagai (JAEA, Riken)
* Kosuke Nogaki (Kyoto Univ.)
* Takuya Nomoto (Univ. Tokyo)
* Hitoshi Mori (Riken)
* Junya Otsuki (Okayama Univ.)
* Soshun Ozaki (Univ. Tokyo)
* Rihito Sakurai (Saitama Univ.)
* Constanze Vogel (TU Wien)
* Niklas Witt (Univ. Hamburg)
* Kazuyoshi Yoshimi (ISSP)
* Hiroshi Shinaoka (Saitama Univ.)
