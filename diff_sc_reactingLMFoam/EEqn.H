// Solving energy equation

{
    volScalarField& he = thermo.he();

//line added for the turbulent Prt
    scalar prt = 1.0;
    extraDict.readIfPresent("turbPrandtl", prt, true);

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + mvConvection->fvmDiv(phi, he)
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), U),
                p,
                "div(phiv,p)"
            )
          : 0*dpdt // low mach number assumption
        )
      - fvm::laplacian((turbulence->alpha())+(turbulence->alphat()/prt), he)
     ==
        reaction->Sh()
      + fvOptions(rho, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();

    Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;
}
