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
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef linkslip_h
#define linkslip_h

#include "structuralinterfacematerial.h"
#include "structuralinterfacematerialstatus.h"

///@name Input fields for LatticeSlip
//@{
#define _IFT_LinkSlip_Name "linkslip"
#define _IFT_LinkSlip_talpha "talpha"
#define _IFT_LinkSlip_kn "kn"
#define _IFT_LinkSlip_kl "kl"
#define _IFT_LinkSlip_type "type"
#define _IFT_LinkSlip_t0 "t0"
#define _IFT_LinkSlip_s1 "s1"
#define _IFT_LinkSlip_s2 "s2"
#define _IFT_LinkSlip_s3 "s3"
#define _IFT_LinkSlip_tf "tf"
#define _IFT_LinkSlip_alpha "alpha"
//@}

namespace oofem {



class LinkSlipStatus : public StructuralInterfaceMaterialStatus
{
protected:

  double kappa=0;

  double tempKappa=0;

public:

    /// Constructor
    LinkSlipStatus(GaussPoint *g);
    /// Destructor
    ~LinkSlipStatus() {}
    
    double  giveKappa(){ return this->kappa; }

    void  letTempKappaBe(const double &v)
    { this->tempKappa = v; }

    void   printOutputAt(FILE *file, TimeStep *tStep) const override;

    const char *giveClassName() const override { return "LinkSlipStatus"; }
    
    virtual void initTempStatus();

    virtual void updateYourself(TimeStep *); // update after new equilibrium state reached

    void saveContext(DataStream &stream, ContextMode mode);

    void restoreContext(DataStream &stream, ContextMode mode);
};

/**
 * This class implements a slip model for concrete for lattice elements.
 */
class LinkSlip : public StructuralInterfaceMaterial
    //
{
protected:

    ///Normal modulus
    double kNormal;

    ///Lateral modulus
    double kLateral;

    int type;
    
    ///Strength for slip component
    double tauMax,tauFinal;
    
    double s1,s2,s3;

    double alpha;

public:

    LinkSlip(int n, Domain *d);

    /// Destructor
    virtual ~LinkSlip();

    virtual bool hasAnalyticalTangentStiffness() const override { return true; }
    
    virtual const char *giveInputRecordName() const { return _IFT_LinkSlip_Name; }
    virtual const char *giveClassName() const { return "LinkSlip"; }

    virtual IRResultType initializeFrom(InputRecord *ir);

    //  virtual void computeStressIndependentStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN, ValueModeType mode);

    virtual void  giveThermalDilatationVector(FloatArray &answer,  GaussPoint *gp,  TimeStep *tStep);


    virtual bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) { return false; }
    
    double evaluateBondStress(const double kappa) const;

    FloatArrayF<3> giveEngTraction_3d(const FloatArrayF<3> &jump, GaussPoint *gp, TimeStep *tStep) const override;
	
    FloatMatrixF<3,3> give3dStiffnessMatrix_Eng(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const override;
    
    virtual Interface *giveInterface(InterfaceType);
                                  

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;


protected:

    virtual int giveIPValue(FloatArray &answer,
                            GaussPoint *gp,
                            InternalStateType type,
                            TimeStep *atTime);
};
} // end namespace oofem

#endif
