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

#include "steelmpf.h"
#include "classfactory.h"
#include "inputrecord.h"
#include "contextioerr.h"
#include "gausspoint.h"

namespace oofem {
REGISTER_Material( SteelMPF );

void SteelMPF::initializeFrom( InputRecord &ir )
{
    StructuralMaterial::initializeFrom( ir );

    IR_GIVE_FIELD( ir, fyP, _IFT_SteelMPF_fyP );
    fyN = fyP;
    IR_GIVE_OPTIONAL_FIELD( ir, fyN, _IFT_SteelMPF_fyN );
    IR_GIVE_FIELD( ir, E0, _IFT_SteelMPF_E0 );
    IR_GIVE_FIELD( ir, bP, _IFT_SteelMPF_bP );
    bN = bP;
    IR_GIVE_OPTIONAL_FIELD( ir, bN, _IFT_SteelMPF_bN );
    IR_GIVE_FIELD( ir, R0, _IFT_SteelMPF_R0 );
    IR_GIVE_FIELD( ir, a1, _IFT_SteelMPF_a1 );
    IR_GIVE_FIELD( ir, a2, _IFT_SteelMPF_a2 );
    IR_GIVE_FIELD( ir, a3P, _IFT_SteelMPF_a3P );
    IR_GIVE_FIELD( ir, a4P, _IFT_SteelMPF_a4P );
    a3N = a3P;
    IR_GIVE_OPTIONAL_FIELD( ir, a3N, _IFT_SteelMPF_a3N );
    a4N = a4P;
    IR_GIVE_OPTIONAL_FIELD( ir, a4N, _IFT_SteelMPF_a4N );

    // curvature parameter
    a1R0 = a1 * R0;
    // yield strains
    epsYieldP = fyP / E0;
    epsYieldN = fyN / E0;
}


MaterialStatus *
SteelMPF::CreateStatus( GaussPoint *gp ) const
{
    return new SteelMPFMaterialStatus( gp, fyP, fyN, E0, R0 );
}


FloatMatrixF<1, 1>
SteelMPF::give1dStressStiffMtrx( MatResponseMode mode,
    GaussPoint *gp,
    TimeStep *tStep ) const
{
    if ( this->isActivated( tStep ) ) {
        // static_cast as I'm the one overriding CreateStatus and the type HAS to be what the one needed
        const auto ms = static_cast<SteelMPFMaterialStatus *>( this->giveStatus( gp ) );
        if ( !ms->hasTangent() ) { ///@todo Make this hack fit more nicely into OOFEM in general;
            // Evaluating the function once, so that the tangent can be obtained.
            this->giveRealStressVector_1d( zeros<1>(), gp, tStep );
        }

        return ms->giveTangentTemp();
    }
    return { 0. };
}


FloatArrayF<1>
SteelMPF::giveRealStressVector_1d( const FloatArrayF<1> &reducedE, GaussPoint *gp, TimeStep *tStep ) const
{
    SteelMPFMaterialStatus *ms = static_cast<SteelMPFMaterialStatus *>( this->giveStatus( gp ) );
    ms->copyCommittedToTemp();

    double eps           = reducedE.at( 1 );
    double Et            = 0;
    double sig           = 0;
    const double &epsOld = ms->giveStrainVector().at( 1 );
    const double &sigOld = ms->giveStressVector().at( 1 );

    SteelMPF_LoadStatus &inc = ms->giveIncTemp();
    int &outP                = ms->giveOutPTemp();
    int &outN                = ms->giveOutNTemp();

    double &eps0P        = ms->giveEps0PTemp();
    double &sig0P        = ms->giveSig0PTemp();
    double &eps0PTwoPrev = ms->giveEps0PTwoPrevTemp();
    double &sig0PTwoPrev = ms->giveSig0PTwoPrevTemp();

    double &eps0N        = ms->giveEps0NTemp();
    double &sig0N        = ms->giveSig0NTemp();
    double &eps0NTwoPrev = ms->giveEps0NTwoPrevTemp();
    double &sig0NTwoPrev = ms->giveSig0NTwoPrevTemp();

    double &epsrP        = ms->giveEpsrPTemp();
    double &sigrP        = ms->giveSigrPTemp();
    double &epsrPTwoPrev = ms->giveEpsrPTwoPrevTemp();
    double &sigrPTwoPrev = ms->giveSigrPTwoPrevTemp();

    double &epsrN        = ms->giveEpsrNTemp();
    double &sigrN        = ms->giveSigrNTemp();
    double &epsrNTwoPrev = ms->giveEpsrNTwoPrevTemp();
    double &sigrNTwoPrev = ms->giveSigrNTwoPrevTemp();

    double &epsrPMax = ms->giveEpsrPMaxTemp();
    double &epsrNMax = ms->giveEpsrNMaxTemp();
    double &RP          = ms->giveRPTemp();
    double &RN          = ms->giveRNTemp();
    int &nLoop          = ms->giveNLoopTemp();

    double &RPTwoPrev = ms->giveRPTwoPrevTemp();
    double &RNTwoPrev = ms->giveRNTwoPrevTemp();

    // Determine initial loading condition
    if ( inc == SteelMPF_LoadStatus::SMP_LS_Unloaded ) {
        eps0P        = epsYieldP;
        sig0P        = fyP;
        eps0PTwoPrev = eps0P;
        sig0PTwoPrev = sig0P;

        eps0N        = -epsYieldN;
        sig0N        = -fyN;
        eps0NTwoPrev = eps0N;
        sig0NTwoPrev = sig0N;

        RP = R0;
        RN = R0;

        RPTwoPrev = R0;
        RNTwoPrev = R0;

        nLoop = eps == 0.0 ? 0 : 1;

        if ( eps > 0 ) {
            inc = SteelMPF_LoadStatus::SMP_LS_PositiveLoading;

            eps0P        = epsYieldP;
            sig0P        = fyP;
            eps0PTwoPrev = eps0P;
            sig0PTwoPrev = sig0P;

            double epsStarP = computeEpsilonStar( eps, eps0P, epsrP );
            sig             = computeSigma( bP, epsStarP, RP, sig0P, sigrP );
            Et              = ( ( sig0P - sigrP ) / ( eps0P - epsrP ) ) * computeEStar( bP, epsStarP, RP );
        } else {
            inc = SteelMPF_LoadStatus::SMP_LS_NegativeLoading;

            eps0N        = -epsYieldN;
            sig0N        = -fyN;
            eps0NTwoPrev = eps0N;
            sig0NTwoPrev = sig0N;

            double epsStarN = computeEpsilonStar( eps, eps0N, epsrN );
            sig             = computeSigma( bN, epsStarN, RN, sig0N, sigrN );
            Et              = ( ( sig0N - sigrN ) / ( eps0N - epsrN ) ) * computeEStar( bN, epsStarN, RN );
        }
    }

    else { // not unloaded
        // committed state variables have been copied already to the temporary state.
        if ( inc == SteelMPF_LoadStatus::SMP_LS_PositiveLoading ) {
            if ( eps < epsOld ) {
                RNTwoPrev    = RN;
                epsrNTwoPrev = epsrN;
                sigrNTwoPrev = sigrN;
                eps0NTwoPrev = eps0N;
                sig0NTwoPrev = sig0N;

                epsrN = epsOld;
                sigrN = sigOld;
                nLoop++;

                epsrNMax = max( epsrN, epsrNMax );

                double sigShift;
                if ( abs( epsrNMax ) < epsYieldP ) {
                    sigShift = 0.0;
                } else {
                    sigShift = max( 0.0, computeShift( fyN, epsrNMax, epsYieldP, a3P, a4P ) );
                }

                double sigy = -fyN - sigShift;
                eps0N       = ( sigy * ( 1.0 - bN ) + E0 * epsrN - sigrN ) / ( E0 * ( 1.0 - bN ) );
                sig0N       = sigrN + E0 * ( eps0N - epsrN );

                double em;
                if ( nLoop == 1 ) {
                    em = eps0N;
                } else if ( nLoop == 2 ) {
                    em = -epsYieldN;
                } else {
                    em = epsrP;
                }

                double ksi = abs( ( em - eps0N ) / epsYieldN );
                RN         = min( RN, R0 - a1R0 * ksi / ( a2 + ksi ) );

                double epsStarN = computeEpsilonStar( eps, eps0N, epsrN );
                sig             = computeSigma( bN, epsStarN, RN, sig0N, sigrN );
                Et              = ( ( sig0N - sigrN ) / ( eps0N - epsrN ) ) * computeEStar( bN, epsStarN, RN );

                double epsStarNLimit = computeEpsilonStar( eps, eps0NTwoPrev, epsrNTwoPrev );
                double sigLimit      = computeSigma( bN, epsStarNLimit, RNTwoPrev, sig0NTwoPrev, sigrNTwoPrev );
                double EtLimit       = ( ( sig0NTwoPrev - sigrNTwoPrev ) / ( eps0NTwoPrev - epsrNTwoPrev ) ) * computeEStar( bN, epsStarNLimit, RNTwoPrev );

                double RControl   = RNTwoPrev;
                double estControl = computeEpsilonStar( epsrN, eps0NTwoPrev, epsrNTwoPrev );
                double sigControl = computeSigma( bN, estControl, RControl, sig0NTwoPrev, sigrNTwoPrev );

                outN = sigrN < sigControl ? 1 : 0;

                if ( epsrN < epsrNTwoPrev && outN == 0 && sig < sigLimit ) {
                    sig = sigLimit;
                    Et  = EtLimit;
                    ms->copyNTwoPrevToCurrentTemp();
                }
            } else {
                double epsStarP = computeEpsilonStar( eps, eps0P, epsrP );
                sig             = computeSigma( bP, epsStarP, RP, sig0P, sigrP );
                Et              = ( ( sig0P - sigrP ) / ( eps0P - epsrP ) ) * computeEStar( bP, epsStarP, RP );

                double epsStarPLimit = computeEpsilonStar( eps, eps0PTwoPrev, epsrPTwoPrev );
                double sigLimit      = computeSigma( bP, epsStarPLimit, RPTwoPrev, sig0PTwoPrev, sigrPTwoPrev );
                double EtLimit       = ( ( sig0PTwoPrev - sigrPTwoPrev ) / ( eps0PTwoPrev - epsrPTwoPrev ) ) * computeEStar( bP, epsStarPLimit, RPTwoPrev );

                if ( epsrP > epsrPTwoPrev && outP == 0 && sig > sigLimit ) {
                    sig = sigLimit;
                    Et  = EtLimit;
                    ms->copyPTwoPrevToCurrentTemp();
                }
            }
        }

        if ( inc == SteelMPF_LoadStatus::SMP_LS_NegativeLoading ) {
            if ( eps > epsOld ) {
                RPTwoPrev    = RP;
                epsrPTwoPrev = epsrP;
                sigrPTwoPrev = sigrP;
                eps0PTwoPrev = eps0P;
                sig0PTwoPrev = sig0P;

                epsrP = epsOld;
                sigrP = sigOld;
                nLoop++;

                epsrPMax = min( epsrP, epsrPMax );

                double sigShift;
                if ( abs( epsrPMax ) < epsYieldN ) {
                    sigShift = 0.0;
                } else {
                    sigShift = max( 0.0, computeShift( fyP, epsrPMax, epsYieldN, a3N, a4N ) );
                }

                double sigy = fyP + sigShift;
                eps0P       = ( sigy * ( 1.0 - bP ) + E0 * epsrP - sigrP ) / ( E0 * ( 1.0 - bP ) );
                sig0P       = sigrP + E0 * ( eps0P - epsrP );

                double em;
                if ( nLoop == 1 ) {
                    em = eps0P;
                } else if ( nLoop == 2 ) {
                    em = epsYieldP;
                } else {
                    em = epsrN;
                }

                double ksi = abs( ( em - eps0P ) / epsYieldP );
                RP         = min( RP, R0 - a1R0 * ksi / ( a2 + ksi ) );

                double epsStarP = computeEpsilonStar( eps, eps0P, epsrP );
                sig             = computeSigma( bP, epsStarP, RP, sig0P, sigrP );
                Et              = ( ( sig0P - sigrP ) / ( eps0P - epsrP ) ) * computeEStar( bP, epsStarP, RP );

                double epsStarPLimit = computeEpsilonStar( eps, eps0PTwoPrev, epsrPTwoPrev );
                double sigLimit      = computeSigma( bP, epsStarPLimit, RPTwoPrev, sig0PTwoPrev, sigrPTwoPrev );
                double EtLimit       = ( ( sig0PTwoPrev - sigrPTwoPrev ) / ( eps0PTwoPrev - epsrPTwoPrev ) ) * computeEStar( bP, epsStarPLimit, RPTwoPrev );

                double RControl       = RPTwoPrev;
                double epsStarControl = computeEpsilonStar( epsrP, eps0PTwoPrev, epsrPTwoPrev );
                double sigControl     = computeSigma( bP, epsStarControl, RControl, sig0PTwoPrev, sigrPTwoPrev );

                outP = sigrP > sigControl ? 1 : 0;

                if ( epsrP > epsrPTwoPrev && outP == 0 && sig > sigLimit ) {
                    sig = sigLimit;
                    Et  = EtLimit;
                    ms->copyPTwoPrevToCurrentTemp();
                }
            } else {
                double epsStarN = computeEpsilonStar( eps, eps0N, epsrN );
                sig             = computeSigma( bN, epsStarN, RN, sig0N, sigrN );
                Et              = ( ( sig0N - sigrN ) / ( eps0N - epsrN ) ) * computeEStar( bN, epsStarN, RN );

                double epsStarNLimit = computeEpsilonStar( eps, eps0NTwoPrev, epsrNTwoPrev );
                double sigLimit      = computeSigma( bN, epsStarNLimit, RNTwoPrev, sig0NTwoPrev, sigrNTwoPrev );
                double EtLimit       = ( ( sig0NTwoPrev - sigrNTwoPrev ) / ( eps0NTwoPrev - epsrNTwoPrev ) ) * computeEStar( bN, epsStarNLimit, RNTwoPrev );

                if ( epsrN < epsrNTwoPrev && outN == 0 && sig < sigLimit ) {
                    sig = sigLimit;
                    Et  = EtLimit;
                    ms->copyNTwoPrevToCurrentTemp();
                }
            }
        }
        if ( eps > epsOld ) {
            inc = SteelMPF_LoadStatus::SMP_LS_PositiveLoading;
        } else if ( eps < epsOld ) {
            inc = SteelMPF_LoadStatus::SMP_LS_NegativeLoading;
        }
    }
    ms->letTempStrainVectorBe( { eps } );
    ms->letTempStressVectorBe( { sig } );
    ms->setTangentTemp( Et );
    return sig;
}

int SteelMPF::giveIPValue( FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep )
{
    const auto ms = static_cast<SteelMPFMaterialStatus *>( this->giveStatus( gp ) );
    if ( type == IST_Undefined ) {
        // The undefined value is used to just dump the entire state vector.
        answer = ms->giveStateVector();
        return 1;
    }
    if ( type == IST_StressTensor ) {
        answer = ms->giveStressVector();
        return 1;
    }
    return StructuralMaterial::giveIPValue( answer, gp, type, tStep );
}

double SteelMPF::computeShift( const double &sigy, const double &epsmax, const double &epsy, const double &a3, const double &a4 )
{
    return sigy * a3 * ( abs( epsmax ) / epsy - a4 );
}


double
SteelMPF::computeEpsilonStar( const double &eps, const double &eps0, const double &epsr )
{
    return ( eps - epsr ) / ( eps0 - epsr );
}


double
SteelMPF::computeSigma( const double &b, const double &epsStar, const double &R, const double &sig0, const double &sigr )
{
    const double sigStar = b * epsStar + ( ( 1.0 - b ) / pow( ( 1.0 + pow( epsStar, R ) ), ( 1.0 / R ) ) ) * epsStar;
    return sigr + sigStar * ( sig0 - sigr );
}


double
SteelMPF::computeEStar( const double &b, const double &eps, const double &R )
{
    return ( b + ( ( 1.0 - b ) / pow( ( 1 + pow( eps, R ) ), ( 1.0 / R ) ) ) * ( 1.0 - pow( eps, R ) / ( 1 + pow( eps, R ) ) ) );
}


bool SteelMPF::hasMaterialModeCapability( MaterialMode mode ) const
{
    return mode == _1dMat;
}


double
SteelMPF::give( int aProperty, GaussPoint *gp ) const
{
    switch ( static_cast<SteelMPF_Properties>( aProperty ) ) {
    case SteelMPF_Properties::SMPF_FYP:
        return this->fyP;
    case SteelMPF_Properties::SMPF_FYN:
        return this->fyN;
    case SteelMPF_Properties::SMPF_E0:
        return this->E0;
    case SteelMPF_Properties::SMPF_BP:
        return this->bP;
    case SteelMPF_Properties::SMPF_BN:
        return this->bN;
    case SteelMPF_Properties::SMPF_R0:
        return this->R0;
    case SteelMPF_Properties::SMPF_A1:
        return this->a1;
    case SteelMPF_Properties::SMPF_A2:
        return this->a2;
    case SteelMPF_Properties::SMPF_A3P:
        return this->a3P;
    case SteelMPF_Properties::SMPF_A4P:
        return this->a4P;
    case SteelMPF_Properties::SMPF_A3N:
        return this->a3N;
    case SteelMPF_Properties::SMPF_A4N:
        return this->a4N;

    default:
        return Material::give( aProperty, gp );
    }
}


#define READ_STREAM_FUNC( x )      \
    if ( !stream.read( x ) ) {     \
        THROW_CIOERR( CIO_IOERR ); \
    }
#define WRITE_STREAM_FUNC( x )     \
    if ( !stream.write( x ) ) {    \
        THROW_CIOERR( CIO_IOERR ); \
    }

SteelMPFMaterialStatus ::SteelMPFMaterialStatus( GaussPoint *g, const double &fyP, const double &fyN, const double &E0, const double &R0 ) :
    StructuralMaterialStatus( g )
{
    // temporary state variables
    EtTemp = E0;

    incTemp = SteelMPF_LoadStatus::SMP_LS_Unloaded;

    RPTwoPrevTemp = 0.0;
    RNTwoPrevTemp = 0.0;

    outPTemp = 0;
    outNTemp = 0;

    epsrPTemp = 0.0;
    sigrPTemp = 0.0;

    epsrNTemp = 0.0;
    sigrNTemp = 0.0;

    epsrPMaxTemp = 0.0;
    epsrNMaxTemp = 0.0;

    eps0PTemp = 0.0;
    sig0PTemp = 0.0;

    eps0NTemp = 0.0;
    sig0NTemp = 0.0;

    epsrNTwoPrevTemp = 0.0;
    sigrNTwoPrevTemp = 0.0;
    eps0NTwoPrevTemp = 0.0;
    sig0NTwoPrevTemp = 0.0;

    epsrPTwoPrevTemp = 0.0;
    sigrPTwoPrevTemp = 0.0;
    eps0PTwoPrevTemp = 0.0;
    sig0PTwoPrevTemp = 0.0;

    RPTemp = 0.0;
    RNTemp = 0.0;

    nLoopTemp = 0;

    // committed state variables
    Et = E0;

    inc = SteelMPF_LoadStatus::SMP_LS_Unloaded;

    RPTwoPrev = 0.0;
    RNTwoPrev = 0.0;

    outP = 0;
    outN = 0;

    epsrP = 0.0;
    sigrP = 0.0;

    epsrN = 0.0;
    sigrN = 0.0;

    epsrPMax = 0.0;
    epsrNMax = 0.0;

    eps0P = 0.0;
    sig0P = 0.0;

    eps0N = 0.0;
    sig0N = 0.0;

    epsrNTwoPrev = 0.0;
    sigrNTwoPrev = 0.0;
    eps0NTwoPrev = 0.0;
    sig0NTwoPrev = 0.0;

    epsrPTwoPrev = 0.0;
    sigrPTwoPrev = 0.0;
    eps0PTwoPrev = 0.0;
    sig0PTwoPrev = 0.0;

    RP = R0;
    RN = R0;

    nLoop = 0;
}

void SteelMPFMaterialStatus ::initTempStatus()
{
    StructuralMaterialStatus::initTempStatus();

    EtTemp = Et;

    incTemp          = SteelMPF_LoadStatus::SMP_LS_Unloaded;
    RPTwoPrevTemp    = RPTwoPrev;
    RNTwoPrevTemp    = RNTwoPrev;
    outPTemp         = 1;
    outNTemp         = 1;
    epsrPTemp        = 0.0;
    sigrPTemp        = 0.0;
    epsrNTemp        = 0.0;
    sigrNTemp        = 0.0;
    epsrPMaxTemp  = 0.0;
    epsrNMaxTemp  = 0.0;
    epsrNTwoPrevTemp = 0.0;
    sigrNTwoPrevTemp = 0.0;
    epsrPTwoPrevTemp = 0.0;
    sigrPTwoPrevTemp = 0.0;
    nLoopTemp        = 0;
    eps0PTemp        = eps0P;
    sig0PTemp        = sig0P;
    eps0NTemp        = eps0N;
    sig0NTemp        = sig0N;
    eps0NTwoPrevTemp = eps0NTwoPrev;
    sig0NTwoPrevTemp = sig0NTwoPrev;
    eps0PTwoPrevTemp = eps0PTwoPrev;
    sig0PTwoPrevTemp = sig0PTwoPrev;
    RPTemp           = RP;
    RNTemp           = RN;
}

void SteelMPFMaterialStatus::copyCommittedToTemp() noexcept
{
    EtTemp = Et;

    incTemp          = inc;
    RPTwoPrevTemp    = RPTwoPrev;
    RNTwoPrevTemp    = RNTwoPrev;
    outPTemp         = outP;
    outNTemp         = outN;
    epsrPTemp        = epsrP;
    sigrPTemp        = sigrP;
    epsrNTemp        = epsrN;
    sigrNTemp        = sigrN;
    epsrPMaxTemp  = epsrPMax;
    epsrNMaxTemp  = epsrNMax;
    eps0PTemp        = eps0P;
    sig0PTemp        = sig0P;
    eps0NTemp        = eps0N;
    sig0NTemp        = sig0N;
    epsrNTwoPrevTemp = epsrNTwoPrev;
    sigrNTwoPrevTemp = sigrNTwoPrev;
    eps0NTwoPrevTemp = eps0NTwoPrev;
    sig0NTwoPrevTemp = sig0NTwoPrev;
    epsrPTwoPrevTemp = epsrPTwoPrev;
    sigrPTwoPrevTemp = sigrPTwoPrev;
    eps0PTwoPrevTemp = eps0PTwoPrev;
    sig0PTwoPrevTemp = sig0PTwoPrev;
    RPTemp           = RP;
    RNTemp           = RN;
    nLoopTemp        = nLoop;
}

void SteelMPFMaterialStatus ::updateYourself( TimeStep *tStep ) noexcept
{
    StructuralMaterialStatus::updateYourself( tStep );

    Et = EtTemp;

    inc          = incTemp;
    RPTwoPrev    = RPTwoPrevTemp;
    RNTwoPrev    = RNTwoPrevTemp;
    outP         = outPTemp;
    outN         = outNTemp;
    epsrP        = epsrPTemp;
    sigrP        = sigrPTemp;
    epsrN        = epsrNTemp;
    sigrN        = sigrNTemp;
    epsrPMax  = epsrPMaxTemp;
    epsrNMax  = epsrNMaxTemp;
    eps0P        = eps0PTemp;
    sig0P        = sig0PTemp;
    eps0N        = eps0NTemp;
    sig0N        = sig0NTemp;
    epsrNTwoPrev = epsrNTwoPrevTemp;
    sigrNTwoPrev = sigrNTwoPrevTemp;
    eps0NTwoPrev = eps0NTwoPrevTemp;
    sig0NTwoPrev = sig0NTwoPrevTemp;
    epsrPTwoPrev = epsrPTwoPrevTemp;
    sigrPTwoPrev = sigrPTwoPrevTemp;
    eps0PTwoPrev = eps0PTwoPrevTemp;
    sig0PTwoPrev = sig0PTwoPrevTemp;
    RP           = RPTemp;
    RN           = RNTemp;
    nLoop        = nLoopTemp;
}

void SteelMPFMaterialStatus::copyPTwoPrevToCurrentTemp() noexcept
{
    epsrPTemp = epsrPTwoPrev;
    sigrPTemp = sigrPTwoPrev;
    eps0PTemp = eps0PTwoPrev;
    sig0PTemp = sig0PTwoPrev;
    RPTemp    = RPTwoPrev;
}

void SteelMPFMaterialStatus::copyNTwoPrevToCurrentTemp() noexcept
{
    epsrNTemp = epsrNTwoPrev;
    sigrNTemp = sigrNTwoPrev;
    eps0NTemp = eps0NTwoPrev;
    sig0NTemp = sig0NTwoPrev;
    RNTemp    = RNTwoPrev;
}

FloatArray
SteelMPFMaterialStatus::giveStateVector() const
{
    return FloatArray{
        static_cast<double>( inc ),
        static_cast<double>( outP ),
        static_cast<double>( outN ),
        static_cast<double>( nLoop ),
        eps0P,
        sig0P,
        eps0PTwoPrev,
        sig0PTwoPrev,
        eps0N,
        sig0N,
        eps0NTwoPrev,
        sig0NTwoPrev,
        epsrP,
        sigrP,
        epsrPTwoPrev,
        sigrPTwoPrev,
        epsrN,
        sigrN,
        epsrNTwoPrev,
        sigrNTwoPrev,
        epsrPMax,
        epsrNMax,
        RP,
        RN,
        RPTwoPrev,
        RNTwoPrev,
    };
}

void SteelMPFMaterialStatus ::saveContext( DataStream &stream, ContextMode mode )
{
    StructuralMaterialStatus::saveContext( stream, mode );

    WRITE_STREAM_FUNC( static_cast<int>( inc ) )
    WRITE_STREAM_FUNC( outP )
    WRITE_STREAM_FUNC( outN )
    WRITE_STREAM_FUNC( nLoop )
    WRITE_STREAM_FUNC( eps0P )
    WRITE_STREAM_FUNC( sig0P )
    WRITE_STREAM_FUNC( eps0PTwoPrev )
    WRITE_STREAM_FUNC( sig0PTwoPrev )
    WRITE_STREAM_FUNC( eps0N )
    WRITE_STREAM_FUNC( sig0N )
    WRITE_STREAM_FUNC( eps0NTwoPrev )
    WRITE_STREAM_FUNC( sig0NTwoPrev )
    WRITE_STREAM_FUNC( epsrP )
    WRITE_STREAM_FUNC( sigrP )
    WRITE_STREAM_FUNC( epsrPTwoPrev )
    WRITE_STREAM_FUNC( sigrPTwoPrev )
    WRITE_STREAM_FUNC( epsrN )
    WRITE_STREAM_FUNC( sigrN )
    WRITE_STREAM_FUNC( epsrNTwoPrev )
    WRITE_STREAM_FUNC( sigrNTwoPrev )
    WRITE_STREAM_FUNC( epsrPMax )
    WRITE_STREAM_FUNC( epsrNMax )
    WRITE_STREAM_FUNC( RP )
    WRITE_STREAM_FUNC( RN )
    WRITE_STREAM_FUNC( RPTwoPrev )
    WRITE_STREAM_FUNC( RNTwoPrev )
}

void SteelMPFMaterialStatus::restoreContext( DataStream &stream, ContextMode mode )
{
    StructuralMaterialStatus::restoreContext( stream, mode );

    int x;
    READ_STREAM_FUNC( x )
    inc = static_cast<SteelMPF_LoadStatus>( x );
    READ_STREAM_FUNC( outP )
    READ_STREAM_FUNC( outN )
    READ_STREAM_FUNC( nLoop )
    READ_STREAM_FUNC( eps0P )
    READ_STREAM_FUNC( sig0P )
    READ_STREAM_FUNC( eps0PTwoPrev )
    READ_STREAM_FUNC( sig0PTwoPrev )
    READ_STREAM_FUNC( eps0N )
    READ_STREAM_FUNC( sig0N )
    READ_STREAM_FUNC( eps0NTwoPrev )
    READ_STREAM_FUNC( sig0NTwoPrev )
    READ_STREAM_FUNC( epsrP )
    READ_STREAM_FUNC( sigrP )
    READ_STREAM_FUNC( epsrPTwoPrev )
    READ_STREAM_FUNC( sigrPTwoPrev )
    READ_STREAM_FUNC( epsrN )
    READ_STREAM_FUNC( sigrN )
    READ_STREAM_FUNC( epsrNTwoPrev )
    READ_STREAM_FUNC( sigrNTwoPrev )
    READ_STREAM_FUNC( epsrPMax )
    READ_STREAM_FUNC( epsrNMax )
    READ_STREAM_FUNC( RP )
    READ_STREAM_FUNC( RN )
    READ_STREAM_FUNC( RPTwoPrev )
    READ_STREAM_FUNC( RNTwoPrev )
}

} // namespace oofem