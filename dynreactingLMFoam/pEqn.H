rho = thermo.rho();

volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));

volVectorField HbyA("HbyA", U);
HbyA = rAU*UEqn.H();

surfaceScalarField phiHbyA
(
    "phiHbyA",
    (
        (fvc::interpolate(rho*HbyA) & mesh.Sf())
      + rhorAUf*fvc::ddtCorr(rho, U, rhoUf)
    )
);

fvc::makeRelative(phiHbyA, rho, U);
fvOptions.makeRelative(fvc::interpolate(rho), phiHbyA);

while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pdEqn
    (
        fvc::ddt(psi) * pReff
      + fvc::div(phiHbyA)
      - fvm::laplacian(rho*rAU, pd)
     ==
        fvOptions(psi, pd, rho.name())
    );

    fvOptions.constrain(pdEqn);

    pdEqn.solve(mesh.solver(pd.select(pimple.finalInnerIter())));

    if (pimple.finalNonOrthogonalIter())
    {
        phi = phiHbyA + pdEqn.flux();
    }
}


#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

U = HbyA - rAU*fvc::grad(pd);
U.correctBoundaryConditions();
fvOptions.correct(U);

{
    rhoUf = fvc::interpolate(rho*U);
    surfaceVectorField n(mesh.Sf()/mesh.magSf());
    rhoUf += n*(fvc::absolute(phi, rho, U)/mesh.magSf() - (n & rhoUf));
}

