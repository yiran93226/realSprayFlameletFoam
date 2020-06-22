/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "AbramzonSirignano.H"
#include "specie.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::tmp<Foam::scalarField> Foam::AbramzonSirignano<CloudType>::calcXc
(
    const label celli
) const
{
    scalarField Xc(this->owner().thermo().carrier().Y().size());

    forAll(Xc, i)
    {
        Xc[i] =
            this->owner().thermo().carrier().Y()[i][celli]
           /this->owner().thermo().carrier().W(i);
    }

    return Xc/sum(Xc);
}


template<class CloudType>
Foam::scalar Foam::AbramzonSirignano<CloudType>::Sh
(
    const scalar Re,
    const scalar Sc
) const
{
    return 2.0 + 0.6*Foam::sqrt(Re)*cbrt(Sc);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::AbramzonSirignano<CloudType>::AbramzonSirignano
(
    const dictionary& dict,
    CloudType& owner
)
:
    PhaseChangeModel<CloudType>(dict, owner, typeName),
    liquids_(owner.thermo().liquids()),
    activeLiquids_(this->coeffDict().lookup("activeLiquids")),
    liqToCarrierMap_(activeLiquids_.size(), -1),
    liqToLiqMap_(activeLiquids_.size(), -1)
{
    if (activeLiquids_.size() == 0)
    {
        WarningInFunction
            << "Evaporation model selected, but no active liquids defined"
            << nl << endl;
    }
    else
    {
        Info<< "Participating liquid species:" << endl;

        // Determine mapping between liquid and carrier phase species
        forAll(activeLiquids_, i)
        {
            Info<< "    " << activeLiquids_[i] << endl;
            liqToCarrierMap_[i] =
                owner.composition().carrierId(activeLiquids_[i]);
        }

        // Determine mapping between model active liquids and global liquids
        const label idLiquid = owner.composition().idLiquid();
        forAll(activeLiquids_, i)
        {
            liqToLiqMap_[i] =
                owner.composition().localId(idLiquid, activeLiquids_[i]);
        }
    }
}


template<class CloudType>
Foam::AbramzonSirignano<CloudType>::AbramzonSirignano
(
    const AbramzonSirignano<CloudType>& pcm
)
:
    PhaseChangeModel<CloudType>(pcm),
    liquids_(pcm.owner().thermo().liquids()),
    activeLiquids_(pcm.activeLiquids_),
    liqToCarrierMap_(pcm.liqToCarrierMap_),
    liqToLiqMap_(pcm.liqToLiqMap_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::AbramzonSirignano<CloudType>::~AbramzonSirignano()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::AbramzonSirignano<CloudType>::calculate
(
    const scalar dt,
    const label celli,
    const scalar Re,
    const scalar Pr,
    const scalar d,
    const scalar nu,
    const scalar T,
    const scalar Ts,
    const scalar pc,
    const scalar Tc,
    const scalarField& X,
    const scalar mass,
    scalarField& dMassPC,
    scalar& mtc
) const
{
    // immediately evaporate mass that has reached critical condition
    if ((liquids_.Tc(X) - T) < small)
    {
        if (debug)
        {
            WarningInFunction
                << "Parcel reached critical conditions: "
                << "evaporating all available mass" << endl;
        }

        forAll(activeLiquids_, i)
        {
            const label lid = liqToLiqMap_[i];
            dMassPC[lid] = great;
        }

        return;
    }

    // construct carrier phase species volume fractions for cell, celli
    const scalarField Xc(calcXc(celli));

    // carrier thermo properties
    scalar muc = 0.0;
    scalar Wc = 0.0;
    scalar Cpc = 0.0;
    forAll(this->owner().thermo().carrier().Y(), i)
    {
        scalar Yc = this->owner().thermo().carrier().Y()[i][celli];
        muc += Yc*this->owner().thermo().carrier().mu(i, pc, Ts);
        Wc += Xc[i]*this->owner().thermo().carrier().W(i);
        Cpc += Yc*this->owner().thermo().carrier().Cp(i, pc, Ts);
    }

    const scalar md = mass;

    const scalar rhod = md*6.0/(pi*d*d*d);

    // time constant [s]
    const scalar taud = rhod*sqr(d + rootVSmall)/(18.0*muc);

    // calculate liquid heat capacity
    scalar Cpd = 0.0;

    // liquid molar weight
    scalar Wd = 0.0;

    // latent heat of vaporization
    //scalar Lv = 0.0;

    forAll(activeLiquids_, i)
    {
        const label lid = liqToLiqMap_[i];
        Wd += X[lid]*liquids_.properties()[lid].W();
    }

    forAll(activeLiquids_, i)
    {
        const label lid = liqToLiqMap_[i];
        Cpd += liquids_.properties()[lid].Cp(pc, T)*X[lid]*liquids_.properties()[lid].W()/Wd;
        //Lv += liquids_.properties()[lid].hl(pc, T)*X[lid]*liquids_.properties()[lid].W()/Wd;
    }

    //scalar mdot = 0.0;
    mtc = 0.0;

    // calculate mass transfer of each specie in liquid
    forAll(activeLiquids_, i)
    {
        const label gid = liqToCarrierMap_[i];
        const label lid = liqToLiqMap_[i];

        // vapour diffusivity [m2/s]
        const scalar Dab = liquids_.properties()[lid].D(pc, Ts);

        // saturation pressure for species i [pa]
        // - carrier phase pressure assumed equal to the liquid vapour pressure
        //   close to the surface
        // NOTE: if pSat > pc then particle is superheated
        // calculated evaporation rate will be greater than that of a particle
        // at boiling point, but this is not a boiling model
        const scalar pSat = liquids_.properties()[lid].pv(pc, T);

        // Schmidt number
        const scalar Sc = nu / (Dab + rootVSmall);

        // Sherwood number
        const scalar Sh = this->Sh(Re, Sc);

        // Nusselt number
        const scalar Nu = 2.0 + 0.6 * Foam::sqrt(Re) * cbrt(Pr);

        // Lewis number
        const scalar Le = Sc / max(Pr, small);

        // surface molar fraction - Raoult's Law
        const scalar Xs = X[lid] * pSat / pc;

        // Knudsen layer thickness [m]
        //const scalar Lk = muc * sqrt(2.0*pi*T*RR/liquids_.properties()[lid].W()) / (Sc*pc);

        // vapor mole fraction at the droplet surface
        //scalar Xsneq = Xs - 2.0*Lk*beta/d;
        scalar Xseq = Xs;

        // vapor mass fraction at the droplet surface
        //scalar Ysneq = Xsneq*liquids_.properties()[lid].W()/Wc;
        //Ysneq = (Ysneq < 0.99999 ? Ysneq : 0.99999);

        scalar Yseq = Xseq / ( Xseq + (1 - Xseq)*(Wc / liquids_.properties()[lid].W()) );
        // scalar Yseq = Xseq *liquids_.properties()[lid].W()/Wc;

        Yseq = (Yseq < 0.99999 ? Yseq : 0.99999);

        // vapor mass fraction in bulk gas
        scalar Yv = this->owner().thermo().carrier().Y()[gid][celli];

        // Spalding transfer number for mass
        // scalar BM = (Ysneq - Yv) / max(small, 1.0 - Ysneq);
        scalar BM = mag( (Yseq - Yv) / max(small, 1.0 - Yseq) );
        scalar FM = Foam::pow(1 + BM, 0.7) / BM * log(1 + BM);
        scalar ShModify = 2 + (Sh - 2) / FM;

        //test 
        // Info<< "BM = " << BM << ", Sh* = " << ShModify << endl;

        //=================calculate BT iteratively================
        // initial heat transfer number 
        // scalar BT = max(small, Cpd * (Tc - Ts) / Lv);
        scalar BT = max(small, liquids_.properties()[lid].Cpg(pc, Ts) * (Tc - T) / liquids_.properties()[lid].hl(pc, T));

        //test 
        // Info<< "BT_0 = " << liquids_.properties()[lid].Cpg(pc, Ts) * (Tc - T) / liquids_.properties()[lid].hl(pc, T) << endl;

        scalar NuModify = 0.0;

        for (label j = 0; j < 50; j++)
        {
            scalar FT = Foam::pow(1 + BT, 0.7) / BT * log(1 + BT);
            NuModify = 2 + (Nu - 2) / max(small, FT);

            //???????????????????????????
            scalar phi = (liquids_.properties()[lid].Cpg(pc, Ts)  / Cpc * ShModify / max(small, NuModify) ) / Le;

            scalar BT_new = max( small, Foam::pow(1 + BM, phi) - 1 );

            if (mag(BT - BT_new) / BT_new < 1e-4)
            {
                //test 
                // Info<< "Breaking time circle number j = " << j << endl;
                // Info<< "BT_new = " << BT << ", FT = " << FT << ", phi = " << phi << ", Nu* = " << NuModify << endl;
                break;
            }
                   
            BT = BT_new;

            //test
            if (j == 49){
                Info<< "Convergence of BT is not achieved..." << endl;            
                // Info<< "BT = " << BT << ", FT = " << FT << ", phi = " << phi << ", Nu* = " << NuModify << endl;
            }                
        }
        
        //dMassPC[lid] = Sh/(3.0*Sc) * (md/taud) * log(1.0 + BM) * dt;
        dMassPC[lid] = ShModify / (3.0 * Sc) * (md / taud) * log(1.0 + BM) * dt;
        scalar mdot =  dMassPC[lid] / dt;
        mtc += mdot/(md*BT) * (3*Pr/NuModify)* (Cpc / Cpd);

        // test
        if (dMassPC[lid] <= 0)
        {
            Info<< "evporation mass for lid number [" << lid << "] = " << dMassPC[lid] << " < 0, thus something goes wrong... " << endl;
        }

        // rate of heat transfer within the liquid droplet [J/s]
        //Ql += mdot*( Cpc*(Tc - Ts)/BT - liquids_.properties()[lid].hl(pc, Ts));
    }
    // test
    // Info<< "mtc = " << mtc << endl;
    
}


template<class CloudType>
Foam::scalar Foam::AbramzonSirignano<CloudType>::dh
(
    const label idc,
    const label idl,
    const scalar p,
    const scalar T
) const
{
    scalar dh = 0;

    typedef PhaseChangeModel<CloudType> parent;
    switch (parent::enthalpyTransfer_)
    {
        case (parent::etLatentHeat):
        {
            dh = liquids_.properties()[idl].hl(p, T);
            break;
        }
        case (parent::etEnthalpyDifference):
        {
            scalar hc = this->owner().composition().carrier().Ha(idc, p, T);
            scalar hp = liquids_.properties()[idl].h(p, T);

            dh = hc - hp;
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unknown enthalpyTransfer type" << abort(FatalError);
        }
    }

    return dh;
}


template<class CloudType>
Foam::scalar Foam::AbramzonSirignano<CloudType>::Tvap
(
    const scalarField& X
) const
{
    return liquids_.Tpt(X);
}


template<class CloudType>
Foam::scalar Foam::AbramzonSirignano<CloudType>::TMax
(
    const scalar p,
    const scalarField& X
) const
{
    return liquids_.pvInvert(p, X);
}


// ************************************************************************* //
