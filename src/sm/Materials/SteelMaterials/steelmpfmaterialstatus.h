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

#ifndef steelmpfmaterialstatus_h
#define steelmpfmaterialstatus_h

#include "matstatus.h"
#include "floatarray.h"
#include "sm/Materials/structuralms.h"
#include "gausspoint.h"
#include "timestep.h"

namespace oofem {
class GaussPoint;
class Domain;

enum class SteelMPF_LoadStatus {
    SMP_LS_Unloaded        = 0,
    SMP_LS_PositiveLoading = 1,
    SMP_LS_NegativeLoading = 2,
};

/**
 * This class implements a structural material status information for SteelMP material.
 * Filippou, F. C., Popov, E. P., Bertero, V. V. (1983). "Effects of Bond Deterioration on Hysteretic Behavior of Reinforced Concrete Joints". Report EERC 83-19, Earthquake Engineering Research Center, University of California, Berkeley. 
 * @author Francesco Pontarin
 */
class SteelMPFMaterialStatus : public StructuralMaterialStatus
{
public:
    SteelMPFMaterialStatus( GaussPoint *g, const double &fyP, const double &fyN,
        const double &E0, const double &R0 );

    void initTempStatus() override;
    void updateYourself( TimeStep *tStep ) noexcept override;

    void copyCommittedToTemp() noexcept;
    void copyPTwoPrevToCurrentTemp() noexcept;
    void copyNTwoPrevToCurrentTemp() noexcept;

    FloatArray giveStateVector() const;
    const double &giveTangentTemp() { return this->EtTemp; }
    void setTangentTemp( const double &E )
    {
        this->EtTemp   = E;
        hasTangentFlag = true;
    }
    bool hasTangent() const { return hasTangentFlag; }

    SteelMPF_LoadStatus &giveIncTemp() { return this->incTemp; }
    const SteelMPF_LoadStatus &giveIncCommitted() { return this->inc; }
    double &giveRPTwoPrevTemp() { return this->RPTwoPrevTemp; }
    double &giveRNTwoPrevTemp() { return this->RNTwoPrevTemp; }
    int &giveOutPTemp() { return this->outPTemp; }
    int &giveOutNTemp() { return this->outNTemp; }
    double &giveEpsrPTemp() { return this->epsrPTemp; }
    double &giveSigrPTemp() { return this->sigrPTemp; }
    double &giveEpsrNTemp() { return this->epsrNTemp; }
    double &giveSigrNTemp() { return this->sigrNTemp; }
    double &giveEpsrPMaxMaxTemp() { return this->epsrPMaxMaxTemp; }
    double &giveEpsrNMaxMaxTemp() { return this->epsrNMaxMaxTemp; }
    double &giveEps0PTemp() { return this->eps0PTemp; }
    double &giveSig0PTemp() { return this->sig0PTemp; }
    double &giveEps0NTemp() { return this->eps0NTemp; }
    double &giveSig0NTemp() { return this->sig0NTemp; }
    double &giveEpsrNTwoPrevTemp() { return this->epsrNTwoPrevTemp; }
    double &giveSigrNTwoPrevTemp() { return this->sigrNTwoPrevTemp; }
    double &giveEps0NTwoPrevTemp() { return this->eps0NTwoPrevTemp; }
    double &giveSig0NTwoPrevTemp() { return this->sig0NTwoPrevTemp; }
    double &giveEpsrPTwoPrevTemp() { return this->epsrPTwoPrevTemp; }
    double &giveSigrPTwoPrevTemp() { return this->sigrPTwoPrevTemp; }
    double &giveEps0PTwoPrevTemp() { return this->eps0PTwoPrevTemp; }
    double &giveSig0PTwoPrevTemp() { return this->sig0PTwoPrevTemp; }
    double &giveRPTemp() { return this->RPTemp; }
    double &giveRNTemp() { return this->RNTemp; }
    int &giveNLoopTemp() { return this->nLoopTemp; }

    void saveContext( DataStream &stream, ContextMode mode ) override;
    void restoreContext( DataStream &stream, ContextMode mode ) override;

    const char *giveClassName() const noexcept override { return "SteelMPMaterialStatus"; }

protected:
    bool hasTangentFlag = false;

    // TRIAL State Variables
    double EtTemp;

    // TRIAL History Variables
    SteelMPF_LoadStatus incTemp;

    int outPTemp;
    int outNTemp;

    int nLoopTemp;

    double eps0PTemp;
    double sig0PTemp;
    double eps0PTwoPrevTemp;
    double sig0PTwoPrevTemp;

    double eps0NTemp;
    double sig0NTemp;
    double eps0NTwoPrevTemp;
    double sig0NTwoPrevTemp;

    double epsrPTemp;
    double sigrPTemp;
    double epsrPTwoPrevTemp;
    double sigrPTwoPrevTemp;

    double epsrNTemp;
    double sigrNTemp;
    double epsrNTwoPrevTemp;
    double sigrNTwoPrevTemp;

    double epsrPMaxMaxTemp;
    double epsrNMaxMaxTemp;

    double RPTemp;
    double RNTemp;

    double RPTwoPrevTemp;
    double RNTwoPrevTemp;

    // CONVERGED State Variables
    double Et;

    // CONVERGED History Variables
    SteelMPF_LoadStatus inc;

    int outP;
    int outN;

    int nLoop;

    double eps0P;
    double sig0P;
    double eps0PTwoPrev;
    double sig0PTwoPrev;

    double eps0N;
    double sig0N;
    double eps0NTwoPrev;
    double sig0NTwoPrev;

    double epsrP;
    double sigrP;
    double epsrPTwoPrev;
    double sigrPTwoPrev;

    double epsrN;
    double sigrN;
    double epsrNTwoPrev;
    double sigrNTwoPrev;

    double epsrPMaxMax;
    double epsrNMaxMax;

    double RP;
    double RN;

    double RPTwoPrev;
    double RNTwoPrev;
};

} // namespace oofem

#endif //steelmpfmaterialstatus_h