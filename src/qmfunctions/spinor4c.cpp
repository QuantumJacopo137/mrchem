/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2023 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
 *
 * This file is part of MRChem.
 *
 * MRChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MRChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MRChem.  If not, see <https://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to MRChem, see:
 * <https://mrchem.readthedocs.io/>
 */

#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

#include <MRCPP/Printer>
#include <MRCPP/trees/FunctionTree.h>
#include <MRCPP/treebuilders/project.h>
#include <MRCPP/utils/spinor_utils.h>

#include "Orbital.h"
#include "orbital_utils.h"
#include "spinor4c.h"

namespace mrchem {
namespace spinor4c {

namespace {

/** @brief Deep copy of a spinor with all lazy prefactors folded into the trees. */
Orbital owned_copy(const Orbital &psi) {
    assert_spinor(psi, false);
    Orbital out;
    mrcpp::deep_copy(out, psi);     // plain copy would share the trees with psi
    fold_c1(out);
    return out;
}

/** @brief Swap tree pointers, prefactors and chunk counters of two components of an owned spinor. */
void swap_components(Orbital &psi, int a, int b) {
    std::swap(psi.CompC[a], psi.CompC[b]);
    std::swap(psi.func_ptr->data.c1[a], psi.func_ptr->data.c1[b]);
    std::swap(psi.func_ptr->data.Nchunks[a], psi.func_ptr->data.Nchunks[b]);
}

} // namespace

/** @brief Create a zero 4c spinor
 *
 * Four complex zero trees, spin() = 0, occ() = 1, c1 = 1 (invariants of assert_spinor).
 * The Orbital(SPIN::Paired, 4) constructor sets occ = 2, which is overwritten here.
 */
Orbital make_spinor() {
    Orbital psi(SPIN::Paired, 4);           // First create an orbital doubly occupied, to take advantage of the constructor
    psi.func_ptr->data.d1[0] = 1.0;         // Then overwrite the occupation with 1, as a spinor requires
    psi.defcomplex();                       // Components are complex-valued
    psi.alloc(4);                           // Four zero trees, one per component
    return psi;
}

/** @brief Project a complex function into one component of a spinor
 *
 * @param prec: projection precision
 * @param psi: spinor created by make_spinor, modified in place (only component `comp` changes)
 * @param comp: component index, see Comp
 * @param f: complex scalar function of the position
 *
 * LMEW: projects straight into the component tree. The CompFunction overload
 * mrcpp::project(out, f, prec, comp) re-allocates and zeroes all other components.
 */
void project(double prec, Orbital &psi, int comp, std::function<ComplexDouble(const mrcpp::Coord<3> &)> f) {
    assert_spinor(psi, false);
    if (comp < 0 or comp > 3) INVALID_ARG_ABORT;
    psi.alloc_comp(comp, true);
    psi.func_ptr->data.c1[comp] = 1.0;      // Projected values are the function itself, no prefactor
    mrcpp::project<3>(prec, *psi.CompC[comp], f);
}

/** @brief Check the spinor invariants without aborting
 *
 * @param psi: spinor to check
 * @param require_unit_c1: also require all prefactors c1[a] == 1 (default)
 * @param why: if given, receives the first violated invariant
 * @returns true if all invariants hold
 */
bool is_spinor(const Orbital &psi, bool require_unit_c1, std::string *why) {
    auto fail = [&](const std::string &msg) {
        if (why != nullptr) *why = msg;
        return false;
    };
    if (psi.Ncomp() != 4) return fail("Spinor must have 4 components");
    if (!psi.iscomplex()) return fail("Spinor must be complex");
    if (psi.func_ptr->data.n1[0] != 0) return fail("Spinor label n1[0] must be 0");
    if (psi.func_ptr->data.d1[0] != 1.0) return fail("Spinor occupation must be 1");
    if (psi.func_ptr->data.conj) return fail("Spinor must not be soft-conjugated");
    for (int a = 0; a < 4; a++) {
        if (psi.CompC[a] == nullptr) return fail("Spinor component is not allocated");
        if (require_unit_c1 and psi.func_ptr->data.c1[a] != ComplexDouble(1.0, 0.0)) return fail("Spinor prefactor c1 not folded");
    }
    return true;
}

/** @brief Abort if a spinor violates the invariants of this namespace, see is_spinor */
void assert_spinor(const Orbital &psi, bool require_unit_c1) {
    std::string why;
    if (!is_spinor(psi, require_unit_c1, &why)) MSG_ABORT(why);
}

/** @brief Fold the lazy prefactors c1[a] into the trees and reset them to 1
 *
 * In place: the caller must own the trees (no aliasing copies).
 * LMEW: norm(), rescale(), gradient() and orthogonalize() ignore c1, so every output is folded.
 */
void fold_c1(Orbital &psi) {
    for (int a = 0; a < psi.Ncomp(); a++) {
        ComplexDouble &c = psi.func_ptr->data.c1[a];
        if (c != ComplexDouble(1.0, 0.0)) {
            psi.CompC[a]->rescale(c);
            c = 1.0;
        }
    }
}

/** @brief beta Psi = (L_a, L_b, -S_a, -S_b) */
Orbital apply_beta(const Orbital &psi) {
    Orbital out = owned_copy(psi);
    out.func_ptr->data.c1[Sa] = -1.0;       // Sign flip of the small block is just a prefactor, no tree arithmetic
    out.func_ptr->data.c1[Sb] = -1.0;
    fold_c1(out);
    return out;
}

/** @brief (L, S) -> (S, L): exchange of the large and small blocks */
Orbital swap_large_small(const Orbital &psi) {
    Orbital out = owned_copy(psi);
    swap_components(out, La, Sa);
    swap_components(out, Lb, Sb);
    return out;
}

/** @brief alpha_k Psi = (sigma_k S, sigma_k L), k = 0,1,2 for x,y,z
 *
 * LMEW: mrcpp::apply_gamma(k+1) gives diag(sigma_k, sigma_k) lazily (tree swaps and c1);
 * swapping the L and S blocks turns it into the off-diagonal alpha_k.
 */
Orbital apply_alpha(const Orbital &psi, int k) {
    if (k < 0 or k > 2) INVALID_ARG_ABORT;
    Orbital out = owned_copy(psi);
    mrcpp::apply_gamma(out, k + 1);         // Works on the trees in place, hence the owned copy
    fold_c1(out);
    return swap_large_small(out);
}

/** @brief Time reversal Theta (a, b, c, d) = (b*, -a*, d*, -c*), Theta^2 = -1
 *
 * LMEW: hard conjugation of the tree values via FunctionTree::map. The soft flag from dagger()
 * is not used because rotate() and multiply() abort on soft-conjugated functions.
 */
Orbital time_reverse(const Orbital &psi) {
    Orbital out = owned_copy(psi);
    for (int a = 0; a < 4; a++) out.CompC[a]->map([](ComplexDouble z) { return std::conj(z); });
    swap_components(out, La, Lb);           // Theta pairs the a and b components of each block...
    swap_components(out, Sa, Sb);           // ...and the minus signs go into the prefactors
    out.func_ptr->data.c1[Lb] = -1.0;
    out.func_ptr->data.c1[Sb] = -1.0;
    fold_c1(out);
    return out;
}

/** @brief out = sum_i coefs[i] terms[i], tree by tree on the union grid
 *
 * @param coefs: complex coefficients
 * @param terms: spinors (all invariants of assert_spinor), not modified
 * @param prec: relative crop of the result, negative means no truncation
 *
 * LMEW: mrcpp::linear_combination skips every term whose *unscaled* tree squared norm is below
 * MachineZero (CompFunction.cpp:767) before the coefficient is applied, so spinors of small
 * amplitude (or small components) are silently dropped. FunctionTree::add_inplace has no such cutoff.
 */
Orbital combine(const std::vector<ComplexDouble> &coefs, const std::vector<Orbital> &terms, double prec) {
    if (coefs.size() != terms.size()) INVALID_ARG_ABORT;
    Orbital out = make_spinor();            // Starts from zero trees, each term is added on top
    for (size_t i = 0; i < terms.size(); i++) {
        assert_spinor(terms[i]);
        for (int a = 0; a < 4; a++) out.CompC[a]->add_inplace(coefs[i], *terms[i].CompC[a]);
    }
    out.crop(prec, false);
    return out;
}

/** @brief || a - b || formed on the union grid
 *
 * The difference is built first and its norm measured directly, so resolution is not limited
 * by the cancellation in ||a||^2 + ||b||^2 - 2 Re<a|b> (about 1e-8). Returns +infinity if the
 * result is not finite, so NaN data can never look like a small distance.
 */
double distance(const Orbital &a, const Orbital &b) {
    Orbital d = combine({1.0, -1.0}, {a, b});
    double n2 = d.getSquareNorm();
    if (!std::isfinite(n2)) return std::numeric_limits<double>::infinity();
    return std::sqrt(n2);
}

/** @brief [Psi_1..Psi_M] -> [Psi_1..Psi_M | Theta Psi_1..Theta Psi_M] (deep copies) */
OrbitalVector kramers_expand(const OrbitalVector &reps) {
    OrbitalVector out;
    const int M = reps.size();
    for (int k = 0; k < M; k++) out.push_back(owned_copy(reps[k]));
    for (int k = 0; k < M; k++) out.push_back(time_reverse(reps[k]));
    return out;
}

/** @brief Phi[k+M] = Theta Phi[k] for a vector of 2M spinors; the partners are regenerated, never updated independently */
void kramers_rebuild_partners(OrbitalVector &Phi) {
    if (Phi.size() % 2 != 0) MSG_ABORT("Kramers set needs an even number of spinors");
    const int M = Phi.size() / 2;
    for (int k = 0; k < M; k++) {
        Orbital partner = time_reverse(Phi[k]);
        partner.func_ptr->rank = k + M;     // Rank follows the slot the partner now occupies
        Phi[k + M] = partner;
    }
}

/** @brief max_k || Phi[k+M] - Theta Phi[k] || over the M pairs
 *
 * Uses distance(), so the resolution is round-off level and non-finite data gives +infinity.
 */
double kramers_pairing_error(const OrbitalVector &Phi) {
    if (Phi.size() % 2 != 0) MSG_ABORT("Kramers set needs an even number of spinors");
    const int M = Phi.size() / 2;
    double err = 0.0;
    for (int k = 0; k < M; k++) {
        double d = distance(Phi[k + M], time_reverse(Phi[k]));
        if (!(d <= err)) err = d; // also propagates infinity: std::max would drop a NaN
    }
    return err;
}

/** @brief Deviation of a 2M x 2M matrix over [reps | partners] from the block form (P Q; -Q* P*)
 *
 * Returns max(|A_22 - P*|, |A_21 + Q*|) in max-norm, with A = (P Q; A_21 A_22).
 * Holds for the overlap matrix and for any time-reversal even operator.
 */
double kramers_matrix_error(const ComplexMatrix &A) {
    if (A.rows() != A.cols() or A.rows() % 2 != 0) MSG_ABORT("Kramers matrix must be square with even dimension");
    const int M = A.rows() / 2;
    ComplexMatrix P = A.topLeftCorner(M, M);
    ComplexMatrix Q = A.topRightCorner(M, M);
    double e1 = (A.bottomRightCorner(M, M) - P.conjugate()).cwiseAbs().maxCoeff();
    double e2 = (A.bottomLeftCorner(M, M) + Q.conjugate()).cwiseAbs().maxCoeff();
    return std::max(e1, e2);
}

/** @brief Smallest eigenvalue of a Hermitian overlap matrix */
double min_overlap_eigenvalue(const ComplexMatrix &S) {
    Eigen::SelfAdjointEigenSolver<ComplexMatrix> es(S, Eigen::EigenvaluesOnly);
    return es.eigenvalues().minCoeff();
}

/** @brief Loewdin orthonormalization behind a linear-dependence guard
 *
 * @param prec: rotation precision, also the tolerance for the final check ||S - 1||_max < prec
 * @param tau: guard threshold, the set is rejected unless lambda_min(S) > tau (M1 default: tau = prec)
 * @param Phi: set to orthonormalize; replaced on success, untouched if rejected
 * @return false if lambda_min(S) <= tau (nothing changed), true on success
 *
 * LMEW: hermitian_matrix_pow silently maps |lambda| < MachineZero to 0, so a rank-deficient set
 * would otherwise come back "orthonormalized". A failed post-check aborts.
 * No pair-preserving step is included: for a Kramers set check kramers_pairing_error afterwards.
 */
bool orthonormalize(double prec, double tau, OrbitalVector &Phi) {
    for (int i = 0; i < Phi.size(); i++) assert_spinor(Phi[i]); // labels n1 != 0 would hide overlaps in calc_overlap_matrix
    ComplexMatrix S = orbital::calc_overlap_matrix(Phi);
    double lmin = min_overlap_eigenvalue(S);
    if (lmin <= tau) {
        MSG_WARN("Overlap eigenvalue " << lmin << " below threshold " << tau);
        return false;
    }
    ComplexMatrix U = orbital::calc_lowdin_matrix(Phi);     // Löwdin: U = S^(-1/2), so the rotated set has overlap 1
    OrbitalVector Psi = orbital::rotate(Phi, U, prec);
    ComplexMatrix S_out = orbital::calc_overlap_matrix(Psi);
    double err = (S_out - ComplexMatrix::Identity(S_out.rows(), S_out.cols())).cwiseAbs().maxCoeff();  // Re-check, since small eigenvalues may have been zeroed silently
    if (err >= prec) MSG_ABORT("Orthonormalization failed: ||S - 1|| = " << err);
    Phi = Psi;
    return true;
}

} // namespace spinor4c
} // namespace mrchem
