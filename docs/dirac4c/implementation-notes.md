# Dirac4c implementation notes (teaching log)

Short explanations for a C++ learner, one section per completed milestone. Explain why implementation choices were made and link to real source symbols. Do not put tutorial paragraphs in production code.

## Entry template

### Milestone N: topic
- **Goal:**
- **Code:** `path:Class::method`
- **C++ concept:** define the concept, explain why this usage is safe
- **Why this design:** alternatives and trade-off
- **Physics invariant:** what must mathematically hold
- **Proof:** focused test and expected/observed value

## Milestone 1: four-component spinor kernel (`spinor4c`, `dirac_free`)

### 1. A spinor is an `Orbital` with four complex trees
- **Goal:** one representation, no new class hierarchy.
- **Code:** `src/qmfunctions/spinor4c.cpp: make_spinor, assert_spinor`
- **C++ concept:** `mrcpp::CompFunction` holds `shared_ptr<TreePtr>`, so copying an `Orbital` is a *shallow* copy: two `Orbital`s then alias the same trees and a mutation through one changes the other. Every function therefore takes `const Orbital &` and returns a fresh object built with `mrcpp::deep_copy` (`owned_copy`). Mutation (`fold_c1`, tree pointer swaps) only ever happens on an object the function itself owns.
- **Why this design:** the alternative (a `Spinor4c` wrapper class) would have to re-implement `dot`, `rotate`, `calc_overlap_matrix` and I/O that already honour `Ncomp = 4`. Invariants (`n1[0] = 0`, `d1[0] = 1`, `c1 = 1`, complex) are instead enforced by `assert_spinor`, which every function calls on input and establishes on output.
- **Physics invariant:** `n1[0] = 0` keeps Kramers partners non-orthogonal-by-label; `occ = 1` because both partners are occupied individually (D14).
- **Proof:** `Spinor construction`; every test calls `assert_spinor` through the functions under test.

### 2. Lazy prefactors `c1` and why `fold_c1` runs on every output
- **Code:** `fold_c1`, `apply_beta`, `apply_alpha`, `time_reverse`
- **C++ concept:** `mrcpp::apply_gamma` implements the Pauli matrices without touching trees: it swaps `CompC[i]` pointers and multiplies `c1[i]` by `±1, ±i`. This is O(1), but `norm`, `rescale`, `gradient` ignore `c1`. `fold_c1` multiplies the prefactor into the tree (`FunctionTree::rescale`) and resets `c1 = 1`.
- **Why this design:** `α_k` = `apply_gamma(k+1)` (which yields `diag(σ_k, σ_k)`) followed by swapping the L and S blocks (`swap_large_small`); β is a sign on `c1[Sa], c1[Sb]`. No new tree arithmetic is needed.
- **Physics invariant:** `α_k² = β² = 1`, `{α_j, α_k} = 2δ_jk`, `{α_k, β} = 0`, α_y signs `(-i, +i, -i, +i)`.
- **Proof:** `Dirac matrices beta and alpha_k` compares against Eigen 4x4 matrices pointwise (|Δ| < 1e-10) and checks the algebra on a spinor.

### 3. Time reversal with hard conjugation
- **Code:** `time_reverse`
- **C++ concept:** the soft flag `dagger()` sets is a bit in `CompFunctionData`; `rotate` and `multiply(ket)` abort on it. `FunctionTree::map(conj)` instead conjugates the point values of a tree we own, so the output has no flag. Complex conjugation is a lambda passed as `FMap<T,T>`.
- **Why this design:** Θ(a,b,c,d) = (b*, -a*, d*, -c*) is a permutation with signs (again tree-pointer swaps plus `c1`) on top of conjugation.
- **Physics invariant:** Θ² = -1, ⟨Ψ|ΘΨ⟩ = 0, ‖ΘΨ‖ = ‖Ψ‖, ⟨ΘΨ|ΘΦ⟩ = ⟨Ψ|Φ⟩*.
- **Proof:** `Time reversal`; pointwise `(b*, -a*, d*, -c*)` to 1e-10, Θ² to 1e-6.

### 4. Kramers sets and the guarded Löwdin step
- **Code:** `kramers_expand`, `kramers_rebuild_partners`, `kramers_pairing_error`, `kramers_matrix_error`, `orthonormalize`
- **Why this design:** a Kramers set is an ordinary `OrbitalVector` with layout `[Ψ_1..Ψ_M | ΘΨ_1..ΘΨ_M]`. Its overlap has block form (P Q; -Q* P*), which `kramers_matrix_error` measures. `S = 1` alone does not prove pairing (multiplying one partner by i keeps S and the block form), so `kramers_pairing_error` compares partners with Θ directly; it equals √2 for that corruption. `hermitian_matrix_pow` silently zeroes tiny eigenvalues, so `orthonormalize` rejects sets with λ_min(S) ≤ τ and verifies ‖S-1‖ < prec afterwards.
- **Proof:** `Kramers set`, `Linear-dependence guard`: Löwdin keeps the pairing error at 3.3e-8 *without* rebuilding partners (measured with the old cancellation-limited formula; the current `distance` gives round-off level); the duplicate representative gives eigenvalues {0, 0, 2, 2}; a nearly dependent set (λ_min = 0.1019) is rejected for τ = 2λ_min and accepted for τ = λ_min/2.

### 5. Free Dirac operator
- **Code:** `src/qmoperators/one_electron/dirac_free.cpp: apply_sigma_p, apply_alpha_p, apply_h0`
- **C++ concept:** `RankZeroOperator::operator()(Orbital, alpha=0)` applies p_k to all four components. σ_k is then applied with `apply_gamma(k+1)` on the owned result of each term, and the three terms are summed with `mrcpp::linear_combination`. `CompFunction` declares a move constructor that is never defined (source-map §1.1): `std::vector::push_back(rvalue)` then fails at *link* time, so the term vectors are built from an initializer list of lvalues (copies).
- **Physics invariant:** `σ·p (u,v) = (p_z u + (p_x - i p_y) v, (p_x + i p_y) u - p_z v)`; h0 Hermitian, Θ-invariant, `h0 → c² + T` for kinetically balanced spinors.
- **Proof:** `dirac_operator` tests: σ·p vs analytic Gaussian derivatives (relative 6.4e-5); `‖σ·pΨ‖² = 2⟨T⟩` (both 23.15444314593); ⟨φ|h0ψ⟩ = ⟨h0φ|ψ⟩ to all printed digits; ‖h0Θψ - Θh0ψ‖/‖·‖ = 1.2e-15; kinetic balance |E - T| = 1.1e-4, 1.1e-6, 2.0e-7 for c = 1e2, 1e3, 1e4 (the c = 1e4 value is dot-product round-off, about c²·2.1e-15; see plan §10.2 and the shifted-form test).

### Review fixes (Codex review of M1)
- **Amplitude cutoff (P1).** `mrcpp::linear_combination` skips any term whose *unscaled* tree squared norm or `|c·c1|²` is below `MachineZero = 1e-14` (`CompFunction.cpp:719,737,766`), before the coefficient is applied. `RankZeroOperator::operator()` uses it too. A spinor scaled by 1e-8 (‖Ψ‖² = 7.5e-16) came back from `apply_h0` as exactly zero. Fix: `spinor4c::combine` sums with `FunctionTree::add_inplace` (union grid, no cutoff), and `apply_p_componentwise` scales each component to unit norm before `p[k]` and back after. Test `Homogeneity at small amplitude`: relative error of `h0(sΨ)` vs `s h0(Ψ)` is ~1e-16 for s = 1e-2…1e-12 (‖sΨ‖² down to 7.5e-24); a 1e-9 small block through σ·p has relative error 1.8e-15.
- **Pairing diagnostic (P2).** `kramers_pairing_error` now uses `spinor4c::distance` (difference on the union grid, norm measured directly) and returns +infinity for non-finite data (`std::max(err, NaN)` would keep the old value). Test: partner `(1 + 1e-9 i)ΘΨ` gives 1.000000e-09; a NaN partner gives infinity. The test helpers `diff_norm` delegate to `distance`.
- **Validation (P3).** `orthonormalize` calls `assert_spinor` on every input. `is_spinor` is the non-aborting form, so tests can check rejection of Alpha/Beta labels, occupation 2, 2 components, soft conjugation and unfolded `c1`. Identical spinors labelled Alpha and Beta are reported by `calc_overlap_matrix` as orthogonal (0) while `dot` gives 1: the reason the check is needed. Abort paths themselves (`MSG_ABORT`) remain untested.
- **Guard at τ = prec (P4).** Test `Guard at tau = prec` tunes the perturbation (λ_min ∝ δ²) to 0.5·prec and 2·prec at prec = 1e-3 and 1e-4: rejected below, accepted above with ‖S−1‖ ≤ 1.5e-9 and pairing error ≤ 2e-13. The rotation error is far below prec because it is exact algebra on the same trees, so no amplification problem showed up at τ = prec; refinement below 1e-4 was not run.

### Test pitfalls found
- `Orbital::dagger()` is a shallow copy and sets the conjugation flag on the *shared* data, so it also flips the flag on the original. Tests that need a soft-conjugated spinor deep-copy first and set `data.conj` directly.
- `orbital::load_orbitals(file)` without `n_orbs` probes a missing `_idx_1` and aborts in `alloc()` ("Function must be defined either real or complex"). Pre-existing, not edited; the I/O test passes `n_orbs = 1`.
- Catch tags with brackets (`[spinor4c]`) are run through `ctest -R spinor4c`.

## Ownership and lifetime of the M1 decisions

Plain summary for the decisions in Milestone 1 above. "Owner" means the object that is responsible for freeing the trees.

| Decision | Source | Owner of the trees | When freed | Alternative and why not |
|---|---|---|---|---|
| Spinor is an `Orbital` | `spinor4c.cpp: make_spinor` | the `Orbital` returned to the caller | when the last `Orbital` sharing `func_ptr` is destroyed (`shared_ptr` reference count reaches 0) | new `Spinor4c` class: would re-implement overlap, rotate and I/O |
| Lazy `c1` | `fold_c1`, `apply_beta`, `time_reverse` | the object that owns the trees, after `owned_copy` | with the tree | fold on every operation: one extra tree pass per term |
| `owned_copy` | `spinor4c.cpp: owned_copy` | the new `Orbital` it returns (deep copy) | when the caller drops it | shallow copy: two objects share trees, so changing one changes the other (Python `b = a` on a list) |
| `swap_components` | `spinor4c.cpp: swap_components` | the object passed by reference, which must already be owned | with that object | copy the trees: O(N) work instead of pointer swap |
| Kramers partners | `kramers_expand`, `kramers_rebuild_partners` | the `OrbitalVector` | with the vector | store partners independently: they could drift from $\Theta\Phi_k$ |
| Löwdin rotation | `orthonormalize` | `Phi` is replaced only on success | with `Phi` | in-place rotation: a rejected set would be lost |
| Operator outputs | `dirac_free.cpp: apply_h0` etc. | returned by value | with the caller | output parameter (`Orbital &out`): more error-prone, no gain with move semantics |

**Test evidence for lifetime.** Each function's output was checked to satisfy `assert_spinor`, and the owned copies were checked not to alias the input (see the test pitfalls below). The `Orbital::dagger()` pitfall is the clearest example of a shared-object bug caught in testing.

**LMEW:** no new `LMEW:` comments were added in this pass. The existing comments already explain the non-obvious choices (cutoffs in `linear_combination`, `apply_gamma` behaviour, `make_density`). No runtime code was changed to teach.
