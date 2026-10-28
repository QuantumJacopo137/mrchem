# Takeaways: the four-component spinor, in plain language

**Status:** Milestone 1 (M1) is implemented and tested, but not committed. M2 to M7 are not started. The step-by-step log is in [`milestones_progress.md`](milestones_progress.md). The detailed C++ notes are in [`implementation-notes.md`](implementation-notes.md). The test numbers below are copied from those notes. I did not re-run the tests for this file.

---

## 1. C++ for a Fortran 90 / Python reader

| You know | In this code | What it means |
|---|---|---|
| module | `namespace` | a named group of functions, e.g. `mrchem::spinor4c` |
| `private` inside a module | `namespace { ... }` | helper visible only inside this `.cpp` file |
| module interface | `.h` file | declarations other files may use |
| module body | `.cpp` file | the actual code |
| `intent(in)` | `const T &` | read only, not copied |
| plain argument that may change | `T &` | the caller's object itself; changes are visible outside |
| derived type | `class` / `struct` | an object with data and functions |
| allocatable array | `std::vector<T>` | resizable list |
| `complex(dp)` | `ComplexDouble` | complex number |
| `x%y` | `x.y`, or `x->y` through a pointer | access a member |
| Python `a = b` on a list | shallow copy | two names, one object |
| `copy.deepcopy(a)` | `mrcpp::deep_copy` | new, independent data |
| Python `lambda` | `[](...) { ... }` | a small unnamed function |

**Indices start at 0 in C++.** `for (int a = 0; a < 4; a++)` runs `a = 0, 1, 2, 3`. In this code that means `La = 0, Lb = 1, Sa = 2, Sb = 3`. Fortran's 1-based habit causes off-by-one errors here.

---

## 2. Why pass by reference or by value, and why `const`

Three rules are used throughout the code:

1. **`const T &`** means "read it, do not copy it, and the compiler will refuse any change". This is the same idea as `intent(in)`.
2. **`T &`** means "do not copy it, you may change it". This is the same as an ordinary Fortran argument.
3. **`T` (no `&`)** means "make a copy". For `Orbital` that copy is *shallow*, which matters a lot (see section 3).

Examples from `spinor4c.cpp`:

```cpp
Orbital apply_beta(const Orbital &psi);          // reads psi, returns a new spinor; psi unchanged
void swap_components(Orbital &psi, int a, int b); // changes psi in place
```

Returning an `Orbital` by value (no `&`) is fine. The compiler moves the object out cheaply, so you do not need to think about it.

---

## 3. Milestone 1: the four-component spinor kernel

**What we built.** A spinor is one `Orbital` with four complex components, stored in the order $(L_\alpha, L_\beta, S_\alpha, S_\beta)$. Large components are $L$, small components are $S$. Files: `src/qmfunctions/spinor4c.{h,cpp}` and `src/qmoperators/one_electron/dirac_free.{h,cpp}`.

### Step 1. Reuse `Orbital`; do not write a new class

- **Why:** MRChem already has overlap, rotation and file I/O for `Orbital`. A new class would need all of them rewritten.
- **Fortran/Python parallel:** using an existing module instead of copying its code.
- **Rules:** `assert_spinor` checks each input: 4 components, complex, label 0, occupation 1, no pending conjugation flag. It is like an `assert` at the top of a subroutine.

### Step 2. Lazy factors `c1`: scale now, multiply later

Each component has a number `c1[a]` stored beside its tree. Changing the sign of the small block only changes `c1`, which is cheap. Before `norm`, `rescale` or `gradient`, the code calls `fold_c1`, which multiplies the factor into the tree.

- **Example:** `Sa` has tree $T$ and `c1[Sa] = -1`. Its true value is $-T$. `fold_c1` replaces the tree with $-T$ and resets the factor to 1.
- **Python parallel:** keeping a scale factor next to an array instead of rescaling every element.

### Step 3. Shallow versus deep copy (`owned_copy`)

- **Shallow copy:** `Orbital out = psi;` gives two names for the same trees. Changing one changes the other. Python: `b = a` for a list.
- **Deep copy:** `mrcpp::deep_copy(out, psi);` gives new trees with the same values. Python: `copy.deepcopy`.
- **`owned_copy`:** deep copy, then `fold_c1`. The result is independent and has no pending factors. Every function that changes a spinor starts from one of these.

```text
a = [1, 2]
b = a          # shallow: b[0] = 5 also changes a
c = deepcopy(a)  # deep: c[0] = 5 leaves a alone
```

### Step 4. Kramers pairs and time reversal $\Theta$

A time-reversal partner is built by $\Theta(a,b,c,d) = (b^*, -a^*, d^*, -c^*)$.

**Toy example.** Take $\Psi = (1, 0, 0, 0)$, so only $L_\alpha = 1$.

- $\Theta\Psi = (0, -1, 0, 0)$. The overlap $\langle\Psi|\Theta\Psi\rangle = 0$, so the partners are orthogonal.
- $\Theta(\Theta\Psi) = (-1, 0, 0, 0) = -\Psi$, so $\Theta^2 = -1$.

A Kramers set is stored as $[\Psi_1 \dots \Psi_M \mid \Theta\Psi_1 \dots \Theta\Psi_M]$. Partners are never updated on their own. They are always rebuilt from their representative, so the pairing cannot drift.

**Check used:** `kramers_pairing_error` measures $\|\Phi_{k+M} - \Theta\Phi_k\|$ directly. An overlap of 1 alone does not prove the pairing. Multiplying one partner by $i$ keeps the overlap matrix the same but breaks the pairing.

### Step 5. Guarded Löwdin orthonormalization

Löwdin rotates a set so that its overlap matrix becomes the identity: $\Phi \to \Phi\, S^{-1/2}$.

**Toy example.** Two identical spinors give $S = \begin{pmatrix}1&1\\1&1\end{pmatrix}$, with eigenvalues $0$ and $2$. The smallest eigenvalue is 0, which is below the threshold $\tau$, so `orthonormalize` refuses and returns `false` without changing the input.

- **Why refuse instead of clamping?** `hermitian_matrix_pow` silently sets tiny eigenvalues to zero. That would hide a real linear dependence. The same policy is used elsewhere in MRChem: abort rather than guess.

### Step 6. The free Dirac operator $h_0 = c\,\boldsymbol\alpha\cdot\mathbf p + \beta c^2$

Built from four pieces: $\mathbf p = -i\nabla$, applied to each component, then the Pauli-matrix pattern.

**Toy example.** Let $(L_\alpha, L_\beta) = (1, 2)$ and look at $\sigma_k$ acting on it:

- $\sigma_x$: $(2, 1)$
- $\sigma_y$: $(-2i, i)$
- $\sigma_z$: $(1, -2)$

These are the operations `mrcpp::apply_gamma` performs by swapping tree pointers and multiplying factors by $\pm 1$ or $\pm i$.

**Tests in `dirac_operator`:**

| Check | Result |
|---|---|
| $\sigma\cdot p$ against an analytic derivative | relative error $6.4\cdot10^{-5}$ |
| $\|\sigma\cdot p\,\Psi\|^2 = 2\langle T\rangle$ | $23.15444314593$ |
| $\langle\phi\|h_0\psi\rangle = \langle h_0\phi\|\psi\rangle$ | equal to all printed digits |
| Time reversal commutes with $h_0$ | $1.2\cdot10^{-15}$ |

### Problems found in review and fixed

**Amplitude cutoff.** MRCPP's `linear_combination` skips a term when its squared norm is below $10^{-14}$, before the coefficient is applied. A spinor scaled by $10^{-8}$ has $\|\Psi\|^2 = 7.5\cdot10^{-16}$. It was silently dropped, and `apply_h0` returned exactly zero.

- **Fix:** `combine` adds the terms with `add_inplace`, which has no cutoff. `apply_p_componentwise` scales each component to unit norm before differentiating, then scales back.
- **Test:** $h_0(s\Psi)$ matches $s\,h_0(\Psi)$ to relative error about $10^{-16}$ for $s$ from $10^{-2}$ to $10^{-12}$.

---

## 4. Open points before M2

1. `apply_gamma` was written for 2-component spinors (its own comment says so). The tests check its output against analytic results, but the function has not been reviewed for 4 components.
2. The sign convention of $\sigma_y$ inside MRCPP's `apply_gamma` is listed as unverified in [`branch-status.md`](branch-status.md) §1B. It should be checked against [`theory.md`](theory.md) §1.
3. M1 is not committed. `spinor4c.cpp` holds all the parts, so splitting it by milestone would need hunk-level staging.

---

## 5. Where things are

- `src/qmfunctions/spinor4c.{h,cpp}`: spinor container, Kramers and Löwdin helpers
- `src/qmoperators/one_electron/dirac_free.{h,cpp}`: $\sigma\cdot p$, $\alpha\cdot p$, $h_0$
- `tests/qmfunctions/spinor4c.cpp`, `tests/qmoperators/dirac_operator.cpp`: tests (ctest labels `spinor4c`, `dirac_operator`)
- `docs/dirac4c/`: this file, the milestone log, detailed notes, and the MRCPP and MRChem source map
