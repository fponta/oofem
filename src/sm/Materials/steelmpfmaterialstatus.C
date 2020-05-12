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

#include "sm/Materials/steelmpfmaterialstatus.h"
#include "sm/Materials/structuralms.h"
#include "sm/Materials/structuralmaterial.h"
#include "contextioerr.h"
#include "gausspoint.h"

namespace oofem {

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

    epsrPMaxMaxTemp = 0.0;
    epsrNMaxMaxTemp = 0.0;

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

    epsrPMaxMax = 0.0;
    epsrNMaxMax = 0.0;

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
    epsrPMaxMaxTemp  = 0.0;
    epsrNMaxMaxTemp  = 0.0;
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
    epsrPMaxMaxTemp  = epsrPMaxMax;
    epsrNMaxMaxTemp  = epsrNMaxMax;
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
    epsrPMaxMax  = epsrPMaxMaxTemp;
    epsrNMaxMax  = epsrNMaxMaxTemp;
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
        epsrPMaxMax,
        epsrNMaxMax,
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
    WRITE_STREAM_FUNC( epsrPMaxMax )
    WRITE_STREAM_FUNC( epsrNMaxMax )
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
    READ_STREAM_FUNC( epsrPMaxMax )
    READ_STREAM_FUNC( epsrNMaxMax )
    READ_STREAM_FUNC( RP )
    READ_STREAM_FUNC( RN )
    READ_STREAM_FUNC( RPTwoPrev )
    READ_STREAM_FUNC( RNTwoPrev )
}
} // namespace oofem
