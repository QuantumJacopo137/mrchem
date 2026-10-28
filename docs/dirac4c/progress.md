# Progress

Milestone: 1 (four-component spinor kernel) implemented and tested, uncommitted in worktree `dirac4c-source-map`. M2+ not started.
M1 (2026-10-28): new `src/qmfunctions/spinor4c.{h,cpp}`, `src/qmoperators/one_electron/dirac_free.{h,cpp}`, tests `tests/qmfunctions/spinor4c.cpp` (tag `spinor4c`) and `tests/qmoperators/dirac_operator.cpp` (tag `dirac_operator`). Both pass at OMP_NUM_THREADS=1 and 2. Full `ctest` (61 tests): same 37 failures as baseline, 0 new. He HF -2.861679945033 reproduced (|dE| = 4.5e-13). Build: 0 warnings. Codex review (4 findings) addressed 2026-10-28, see implementation-notes. Caveats investigated (plan §10): derivative error depends on order and operator; rest-mass round-off and crop amplification measured; proposed decisions D19–D21 await user confirmation. Raw logs in `logs/m1/`.
(Earlier status, milestone 0:) reconnaissance of `DHF_Quentin_based` partially done; no 4c code written.
Verified source SHA: MRChem `feb232a0` (branch `DHF_Quentin_based`); MRCPP pinned `Dheasra/mrcpp@6cc41b9b` (branch), `@ed46025e` (ASC); ⚠ `build/_deps` holds upstream `f73b863d`, reconfigure needed
Branch reconnaissance: `branch-status.md` (2026-10-27)
Tests and raw log paths (2026-10-27, OpenMP-only, OMP_NUM_THREADS=1, `ctest -j 12`):
- `master` `c6a6f174` + MRCPP `f73b863d` (worktree `../mrchem-master-baseline/build-omp`): **59/59 pass** (`logs/master-ctest-omp1.log`)
- `DHF_Quentin_based` `feb232a0` + `Dheasra/mrcpp@6cc41b9b` (`build-omp/`): **22/59 pass, 37 fail** (`logs/ctest-omp1.log`): 5 unit tests abort (operator_composition: `QMPotential::calc` 'Output not empty'; reaction_operator, PB_solver: MRCPP `CompFunction.cpp:989` multiply 'Not implemented'; xc_hessian_lda/pbe), 12 integration aborts, 20 energy deviations (1e-7 … 1.3 Ha). All are regressions w.r.t. master.
Next action: reconfigure build against pinned MRCPP, run baseline OpenMP build + `ctest`; ASC fetched (`ac3755f0`), CompFunction verified (branch-status §5); interview done, decisions D01–D18 in `decisions.md`

## Docs map (this folder)

- `progress.md`: this file, current status and test logs
- `milestones_progress.md`: step-by-step log of each milestone, with tests
- `takeaways.md`: short, plain-language summary of what was built and why (for a Fortran 90 / Python reader)
- `implementation-notes.md`: detailed C++ decisions per milestone (source paths, ownership, alternatives, tests)
- `architecture.md`, `source-map.md`, `branch-status.md`: MRChem and MRCPP structure and verification
- `scope.md`, `theory.md`, `plan.md`, `decisions.md`: physics, plan and decision log
