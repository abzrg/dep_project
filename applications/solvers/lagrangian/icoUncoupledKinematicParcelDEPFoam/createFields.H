#include "readGravitationalAcceleration.H"

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

#include "createPhi.H"

singlePhaseTrahsportModel laminarTransport(U, phi);

dimensionedScalar rhoInfValue
(
    "rhoInf",
    dimDensity,
    laminarTransport
);

volScalarField rhoInf
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    rhoInfValue
);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

volScalarField mu
(
    IOobject
    (
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    laminarTransport.nu()*rhoInfValue
);

// Initialize the 'basicKinematicCollidingCloud' (includes particle-particle
// interactions):
const word kinematicCloudName
(
    args.getOrDefault<word>("cloud", "kinematicCloud")
);

Info<< "Constructing kinematicCloud " << kinematicCloudName << endl;
basicKinematicCollidingCloud kinematicCloud
(
    kinematicCloudName,
    rhoInf,
    U,
    mu,
    g
);

IOobject Eheader
(
    "E",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
);

autoPtr<volVectorField> EPtr;

if (Eheader.typeHeaderOk<volVectorField>(true))
{
    Info<< "\nReading field E\n" << endl;

    EPtr.reset(new volVectorField (Eheader, mesh));
}

IOobject EdotGradEheader
(
    "EdotGradE",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
);

autoPtr<volVectorField> EdotGradEPtr;

if (EdotGradEheader.typeHeaderOk<volVectorField>(true))
{
    Info<< "Reading field EdotGradE" << endl;

    EdotGradEPtr.reset(new volVectorField(EdotGradEheader, mesh));
}

#include "createNonInertialFrameFields.H"
