Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field V (Electrical Potential)\n" << endl;
volScalarField V
(
    IOobject
    (
        "V",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Electric field
volVectorField E
(
    IOobject
    (
        "E",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,    // There is nothing to read
        IOobject::AUTO_WRITE
    ),
    mesh,
    //-fvc::grad(V)
    dimensionedVector(Foam::dimensionSet(1, 1, -3, 0, 0, -1, 0), One) // [m*kg/s3/A1] = [V/m]
);


// Electric Field strength
volVectorField EdotGradE
(
    IOobject
    (
        "EdotGradE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    //2*E & fvc::grad(E)
    dimensionedVector(Foam::dimensionSet(2, 1, -6, 0, 0, -2, 0), Zero) // [kg2m/s6/A2] = [V2/m3]
);


#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, simple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

/*

If you want to add an entry to the transportProperties dictionary do so like
following

// Schmidt number
dimensionedScalar Sc("Sc", dimless, laminarTransport);

// Turbulent Schmidt number
dimensionedScalar Prt("Sct", dimless, laminarTransport);

*/

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

#include "createMRF.H"
#include "createFvOptions.H"
