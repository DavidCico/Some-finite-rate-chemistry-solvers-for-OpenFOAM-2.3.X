# Some-solvers-for-OpenFOAM-2.3.X
This repository contains different solvers for OpenFOAM 2.3.X

## Files description

<ul>
    <li><p align="justify">"<em>diff_sc_reactingLMFoam</em>"</p></li>
    <li><p align="justify">"<em>dynreactingLMFoam</em>"</p></li>
    <li><p align="justify">"<em>reactingBoxturbFoam</em>"</p></li>
</ul>


## Prerequisites

<p align="justify">This library is developed for <strong>OpenFOAM 2.3.X</strong>, and requires the latter version of the Open Source CFD Toolbox. For more information on how to install this version, the reader is referred to <a href="https://sites.google.com/site/foamguides/installation/installing-openfoam-2-3-x">this web page</a>.</p>

## Instructions on program

Compile the library using:

    wmake libso

This will create the library <em>libInflowTurbulence.so</em> in your $FOAM_USER_LIBBIN folder.

Add this to the case controlDict:

    libs ("libInflowTurbulence.so");

And specify for U at inlet with settings fitting your case

    inlet
    {
      overlap         0.5;    // How much the vortons can overlap each other
        L               0.01;    // Integral length scale
        eta             1e-04;  // Kolmogorov length
        Cl              6.783;  // Constant used for E(k) 6.783 in paper
        Ceta            0.4;  // Constant used for E(k) 0.4 in paper
        type            inflowGenerator<homogeneousTurbulence>;
        fluctuationScale (1 1 1); // Fluctuation scale in each direction
        referenceField  uniform (10 0 0); // Reference value for inlet
        R               uniform (1 0 0 1 0 1); // Reynolds stress tensor
        value           uniform (10 0 0); 
    }

