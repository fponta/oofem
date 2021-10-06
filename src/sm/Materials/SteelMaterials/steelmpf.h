/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef steelmpf_h
#define steelmpf_h

#include "sm/Materials/structuralmaterial.h"

#define _IFT_SteelMPF_Name "steelMPF"
#define _IFT_SteelMPF_fyP "fyp"
#define _IFT_SteelMPF_fyN "fyn"
#define _IFT_SteelMPF_E0 "e0"
#define _IFT_SteelMPF_bP "bp"
#define _IFT_SteelMPF_bN "bn"
#define _IFT_SteelMPF_R0 "r0"
#define _IFT_SteelMPF_a1 "a1"
#define _IFT_SteelMPF_a2 "a2"
#define _IFT_SteelMPF_a3P "a3p"
#define _IFT_SteelMPF_a4P "a4p"
#define _IFT_SteelMPF_a3N "a3n"
#define _IFT_SteelMPF_a4N "a4n"

namespace oofem {

enum class SteelMPF_Properties : int {
    SMPF_FYP = 300,
    SMPF_FYN = 301,
    SMPF_E0  = 302,
    SMPF_BP  = 303,
    SMPF_BN  = 304,
    SMPF_R0  = 305,
    SMPF_A1  = 306,
    SMPF_A2  = 307,
    SMPF_A3P = 308,
    SMPF_A4P = 309,
    SMPF_A3N = 310,
    SMPF_A4N = 311,
};

class GaussPoint;
class Domain;

/**
 * Oofem port of Kolozvari's implementation of the steel material model described in
 * Filippou, F. C., Popov, E. P., Bertero, V. V. (1983). "Effects of Bond Deterioration on Hysteretic Behavior of Reinforced Concrete Joints". Report EERC 83-19, Earthquake Engineering Research Center, University of California, Berkeley. 
 * uses a 2 step memory to mitigate stress overshooting.
 * @author Francesco Pontarin
 */
class SteelMPF : public StructuralMaterial
{
public:
    // Constructor
    SteelMPF( int n, Domain *d ) :
        StructuralMaterial( n, d ) {}

    MaterialStatus *CreateStatus( GaussPoint *gp ) const override;

    FloatArrayF<1> giveRealStressVector_1d( const FloatArrayF<1> &reducedE, GaussPoint *gp, TimeStep *tStep ) const override;

    FloatMatrixF<1, 1> give1dStressStiffMtrx( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep ) const override;

    int giveIPValue( FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep ) override;

    double give( int, GaussPoint *gp ) const override;

    bool hasMaterialModeCapability( MaterialMode mode ) const override;

    bool isCharacteristicMtrxSymmetric( MatResponseMode rMode ) const override { return true; }
    const char *giveClassName() const override { return "SteelMPF"; }
    const char *giveInputRecordName() const override { return _IFT_SteelMPF_Name; }
    void initializeFrom( InputRecord &ir ) override;

protected:
    // Helper functions to calculate stresses and strains according to the formulation
    /**
     * helper function to calculate @f$ \varepsilon^{*} @f$
     */
    static double computeEpsilonStar( const double &eps, const double &eps0P, const double &epsrP );
    /**
     * helper function to calculate the stress according to the material law
     */
    static double computeSigma( const double &b, const double &epsStar, const double &R, const double &sig0, const double &sigr );
    /**
     * helper function to calculate @f$ E^{*} @f$
     */
    static double computeEStar( const double &b, const double &eps, const double &R );
    /**
     * helper function to calculate the stress shift
     */
    static double computeShift( const double &sigy, const double &epsMax, const double &epsy, const double &a3, const double &a4 );

    // Material properties
    double fyP = 0.0; // : positive yield stress
    double fyN = 0.0; // : negative yield stress
    double E0  = 0.0; // : initial stiffness
    double bP  = 0.0; // : hardening ratio (Esh/E0)
    double bN  = 0.0; // : hardening ratio (Esh/E0)
    double R0  = 0.0; // : elastic-plastic transition
    double a1  = 0.0; // : coefficient a1 for changing the transition curvature
    double a2  = 0.0; // : coefficient a2 for changing the transition curvature
    double a3P = 0.0; // : coefficient a3 for hardening in tension
    double a4P = 0.0; // : coefficient a4 for hardening in tension
    double a3N = 0.0; // : coefficient a3 for hardening in compression
    double a4N = 0.0; // : coefficient a4 for hardening in compression

    double a1R0      = 0.0; // : a1 * R0
    double epsYieldP = 0.0; // : yield strain in tension
    double epsYieldN = 0.0; // : yield strain in compression
};

} // namespace oofem

#endif //steelmpf_h