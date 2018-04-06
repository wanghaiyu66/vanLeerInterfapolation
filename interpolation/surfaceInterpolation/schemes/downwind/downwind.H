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
    Foam::downwind

Description
    Downwind differencing scheme class.

SourceFiles
    downwind.C

\*---------------------------------------------------------------------------*/

#ifndef downwind_H
#define downwind_H

#include "surfaceInterpolationScheme.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class downwind Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class downwind
:
    public surfaceInterpolationScheme<Type>
{
    // Private data

        const surfaceScalarField& faceFlux_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const downwind&);


public:

    //- Runtime type information
    TypeName("downwind");


    // Constructors

        //- Construct from faceFlux
        downwind
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            faceFlux_(faceFlux)
        {}

        //- Construct from Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        downwind
        (
            const fvMesh& mesh,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            faceFlux_
            (
                mesh.lookupObject<surfaceScalarField>
                (
                    word(is)
                )
            )
        {}

        //- Construct from faceFlux and Istream
        downwind
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream&
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            faceFlux_(faceFlux)
        {}


    // Member Functions

        //- Return the interpolation weighting factors
        virtual tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const
        {
            return neg(faceFlux_);
        }
};


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
downwindInterpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& faceFlux
)
{
    const word name("base1");
    return surfaceInterpolationScheme<Type>::New
    (
        vf.mesh(),
        faceFlux,
        vf.mesh().interpolationScheme(name)
    )().interpolate(vf);
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
downwindInterpolate
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf,
    const tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >& tsf
)
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tinterp =
        downwindInterpolate(tvf(), tsf());
    tvf.clear();
    tsf.clear();
    return tinterp;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
