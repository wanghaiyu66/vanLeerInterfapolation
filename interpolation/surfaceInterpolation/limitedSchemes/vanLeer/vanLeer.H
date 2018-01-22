/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Class
    Foam::vanLeerLimiter

Description
    Class with limiter function which returns the limiter for the
    vanLeer differencing scheme based on r obtained from the LimiterFunc
    class.

    Used in conjunction with the template class LimitedScheme.

SourceFiles
    vanLeer.C

\*---------------------------------------------------------------------------*/

#ifndef vanLeer_H
#define vanLeer_H

#include "vector.H"
#include "limitedSurfaceInterpolationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class vanLeerLimiter Declaration
\*---------------------------------------------------------------------------*/

template<class LimiterFunc>
class vanLeerLimiter
:
    public LimiterFunc
{

public:

    vanLeerLimiter(Istream&)
    {}

    scalar limiter
    (
        const scalar,
        const scalar faceFlux,
        const typename LimiterFunc::phiType& phiP,
        const typename LimiterFunc::phiType& phiN,
        const typename LimiterFunc::gradPhiType& gradcP,
        const typename LimiterFunc::gradPhiType& gradcN,
        const vector& d
    ) const
    {
        scalar r = LimiterFunc::r
        (
            faceFlux, phiP, phiN, gradcP, gradcN, d
        );

        return (r + mag(r))/(1 + mag(r));
    }
};

/*
template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
vanLeerInterpolate(const GeometricField<Type, fvPatchField, volMesh>& vf)
{
    const word name("base");
    //const word schemeName("Gauss vanLeer");
    return limitedSurfaceInterpolationScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().interpolationScheme(name)
    )().interpolate(vf);

}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
vanLeerInterpolate(const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf)
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tinterp =
        vanLeerInterpolate(tvf());
    tvf.clear();
    return tinterp;
}
*/

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
vanLeerInterpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& faceFlux
)
{
    const word name("base");
    //const word schemeName("Gauss vanLeer");
    return limitedSurfaceInterpolationScheme<Type>::New
    (
        vf.mesh(),
        faceFlux,
        vf.mesh().interpolationScheme(name)
    )().interpolate(vf);
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
vanLeerInterpolate
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf,
    const tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >& tsf
)
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tinterp =
        vanLeerInterpolate(tvf(), tsf());
    tvf.clear();
    tsf.clear();
    return tinterp;
}





// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
