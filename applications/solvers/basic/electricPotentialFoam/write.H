    if (runTime.writeTime())
    {
        // Eelctric field
        volVectorField E
        (
            IOobject
            (
                "E",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            -fvc::grad(V)
        );

        // Electric field strength
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
            2*E & fvc::grad(E)
        );

        runTime.write();
    }
