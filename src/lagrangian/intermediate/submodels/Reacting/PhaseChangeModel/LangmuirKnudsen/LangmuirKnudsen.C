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

#include "LangmuirKnudsen.H"
#include "specie.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::tmp<Foam::scalarField> Foam::LangmuirKnudsen<CloudType>::calcXc
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
Foam::scalar Foam::LangmuirKnudsen<CloudType>::Sh
(
    const scalar Re,
    const scalar Sc
) const
{
    return 2.0 + 0.552*Foam::sqrt(Re)*cbrt(Sc);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::LangmuirKnudsen<CloudType>::LangmuirKnudsen
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
Foam::LangmuirKnudsen<CloudType>::LangmuirKnudsen
(
    const LangmuirKnudsen<CloudType>& pcm
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
Foam::LangmuirKnudsen<CloudType>::~LangmuirKnudsen()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::LangmuirKnudsen<CloudType>::calculate
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

    forAll(activeLiquids_, i)
    {
        const label lid = liqToLiqMap_[i];
        Wd += X[lid]*liquids_.properties()[lid].W();
    }

    forAll(activeLiquids_, i)
    {
        const label lid = liqToLiqMap_[i];
        Cpd += liquids_.properties()[lid].Cp(pc, T)*X[lid]*liquids_.properties()[lid].W()/Wd;
    }

    scalar beta = 0.0;

    scalar mdot = 0.0;

    // iterative method
    for (label j=0; j<50; j++)
    {
        scalar mdotDash = mdot;

        // non-dimensional evaporation parameter
        beta = 1.5*Pr*taud * (mdot/md);

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
            const scalar Sc = nu/(Dab + rootVSmall);

            // Sherwood number
            const scalar Sh = this->Sh(Re, Sc);

            // surface molar fraction - Raoult's Law
            const scalar Xs = X[lid]*pSat/pc;

            // Knudsen layer thickness [m]
            const scalar Lk = muc * sqrt(2.0*pi*T*RR/liquids_.properties()[lid].W()) / (Sc*pc);

            // vapor mole fraction at the droplet surface
            scalar Xsneq = Xs - 2.0*Lk*beta/d;

            // vapor mass fraction at the droplet surface
            scalar Ysneq = Xsneq*liquids_.properties()[lid].W()/Wc;
            Ysneq = (Ysneq < 0.99999 ? Ysneq : 0.99999);

            // vapor mass fraction in bulk gas
            scalar Yv = this->owner().thermo().carrier().Y()[gid][celli];

            // Spalding transfer number for mass
            scalar BM = (Ysneq - Yv) / max(small, 1.0 - Ysneq);

            dMassPC[lid] = Sh/(3.0*Sc) * (md/taud) * log(1.0 + BM) * dt;
        }

        mdot = sum(dMassPC) / dt;

        if (mag(mdot - mdotDash)/max(small, mdotDash) < 1e-4) break;
    }

    beta = max(small, beta);

    mtc = (( beta / (exp(beta) - 1.0) ) / taud) * (Cpc / Cpd);
}


template<class CloudType>
Foam::scalar Foam::LangmuirKnudsen<CloudType>::dh
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
Foam::scalar Foam::LangmuirKnudsen<CloudType>::Tvap
(
    const scalarField& X
) const
{
    return liquids_.Tpt(X);
}


template<class CloudType>
Foam::scalar Foam::LangmuirKnudsen<CloudType>::TMax
(
    const scalar p,
    const scalarField& X
) const
{
    return liquids_.pvInvert(p, X);
}


// ************************************************************************* //
