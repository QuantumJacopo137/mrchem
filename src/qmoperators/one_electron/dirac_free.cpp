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

#include <cmath>
#include <vector>

#include <MRCPP/utils/CompFunction.h>
#include <MRCPP/utils/spinor_utils.h>

#include "dirac_free.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/spinor4c.h"

namespace mrchem {
namespace dirac {

namespace {

/** @brief p_k applied to every component, independent of the amplitude of the components
 *
 * LMEW: RankZeroOperator::operator() assembles its result with mrcpp::linear_combination, which skips a
 * tree whose unscaled squared norm is below MachineZero (CompFunction.cpp:767). Each component is
 * therefore scaled to unit norm before the derivative and scaled back afterwards (derivatives of the
 * unit-norm function are O(1), far above the cutoff). All-zero components stay zero.
 */
Orbital apply_p_componentwise(MomentumOperator &p, int k, const Orbital &psi) {
    double s[4];
    Orbital unit;
    mrcpp::deep_copy(unit, psi);            // Copy, so psi is not rescaled
    for (int a = 0; a < 4; a++) {
        s[a] = std::sqrt(psi.CompC[a]->getSquareNorm());
        if (s[a] > 0.0) unit.CompC[a]->rescale(ComplexDouble(1.0 / s[a], 0.0));
    }
    Orbital out = p[k](unit);
    spinor4c::fold_c1(out);                 // p_k sets lazy prefactors, fold before rescaling back
    for (int a = 0; a < 4; a++) {
        if (s[a] > 0.0) out.CompC[a]->rescale(ComplexDouble(s[a], 0.0));  // Restore the original amplitude
    }
    return out;
}

} // namespace

/** @brief diag(sigma.p, sigma.p) Psi = (sigma.p L, sigma.p S)
 *
 * @param p: MomentumOperator, set up by the caller
 * @param psi: 4c spinor (spinor4c invariants), not modified
 * @returns new spinor, with sigma.p (u, v) = (p_z u + (p_x - i p_y) v, (p_x + i p_y) u - p_z v)
 *
 * LMEW: p_k acts on every component (alpha = 0 is the identity in RankZeroOperator),
 * so sigma_k is applied afterwards with apply_gamma(k+1) on the owned result of each term.
 * The three terms are summed with spinor4c::combine (no amplitude cutoff).
 */
Orbital apply_sigma_p(MomentumOperator &p, const Orbital &psi) {
    spinor4c::assert_spinor(psi);
    std::vector<Orbital> terms;
    terms.reserve(3);
    for (int k = 0; k < 3; k++) {
        Orbital t = apply_p_componentwise(p, k, psi);   // p_k on every component
        mrcpp::apply_gamma(t, k + 1);                   // Multiply by sigma_k (tree swaps, lazy prefactors)
        spinor4c::fold_c1(t);
        terms.push_back(t);
    }
    return spinor4c::combine({1.0, 1.0, 1.0}, terms);   // sum over k of sigma_k p_k
}

/** @brief alpha.p Psi = (sigma.p S, sigma.p L): the L and S blocks of apply_sigma_p are exchanged */
Orbital apply_alpha_p(MomentumOperator &p, const Orbital &psi) {
    return spinor4c::swap_large_small(apply_sigma_p(p, psi));   // Off-diagonal blocks: sigma.p -> alpha.p
}

/** @brief h0 Psi = c alpha.p Psi + c^2 beta Psi = (c^2 L + c sigma.p S, c sigma.p L - c^2 S)
 *
 * @param prec: precision of the final linear combination
 * @param p: MomentumOperator, set up by the caller
 * @param c: speed of light in atomic units
 * @param psi: 4c spinor, not modified
 */
Orbital apply_h0(double prec, MomentumOperator &p, double c, const Orbital &psi) {
    std::vector<Orbital> terms{apply_alpha_p(p, psi), spinor4c::apply_beta(psi)};   // Free Dirac: kinetic (alpha.p) and mass (beta) terms
    return spinor4c::combine({ComplexDouble(c, 0.0), ComplexDouble(c * c, 0.0)}, terms, prec);
}

} // namespace dirac
} // namespace mrchem
