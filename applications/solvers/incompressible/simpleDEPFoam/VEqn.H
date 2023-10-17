    // Electrical potential laplacian equation

    fvScalarMatrix VEqn
    (
        fvm::laplacian(V)
    );

    VEqn.solve();


    // Electric field strength

    E = -fvc::grad(V);
    EdotGradE = 2*(E & fvc::grad(E));
