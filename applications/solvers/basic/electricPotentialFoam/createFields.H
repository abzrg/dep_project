Info<< "Reading field V\n" << endl;

// Electric potential
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
