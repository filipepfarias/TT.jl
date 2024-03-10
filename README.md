# TT.jl

[![Build Status](https://github.com/filipepfarias/TT.jl/actions/workflows/CI.yml/badge.svg?branch=dev)](https://github.com/filipepfarias/TT.jl/actions/workflows/CI.yml?query=branch%3Adev)
[![codecov](https://codecov.io/gh/filipepfarias/TT.jl/graph/badge.svg?token=K1MPRAZWDH)](https://codecov.io/gh/filipepfarias/TT.jl)

### This package is a Julia version of [TT-Toolbox](https://github.com/oseledets/TT-Toolbox) for MATLAB .

### <span style="color:red">This package is under development!

The following functions are partially implemented and not fully tested. Instabilities may occur.

- [ ] funcrs
  - [x] full_lu
  - [x] reort
  - [x] maxvol2
  - [x] tt_x
  - [x] tt_ones
  - [x] mtimes
  - [x] plus
  - [x] minus
  - [x] tt_mem
  - [x] tt_ranks
  - [x] tt_size
  - [x] norm
  - [x] dot
  - [x] rounded_diagonal_index (my_chop2)

- [ ] multifuncrs
  - [ ] qr (lr)
  - [ ] sum
  - [ ] tt_rand
  - [ ] core_to_vector (core2cell)
  - [ ] vector_to_core (cell2core)

- [ ] tAMEn <5>
  - [ ] tt_matrix <4>
  - [ ] grumble_vector <3>
  - [ ] grumble_matrix <3>
  - [ ] parse_matrix <2>
  - [ ] parse_rhs <2>
  - [ ] amenany_sweepzz <2>
  - [ ] chebdiff <1>
  - [ ] cheb2_interpolant
Complexity to implement and test: <1-5>

- [ ] Documentation (TO-DO)