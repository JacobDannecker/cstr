
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

Info<< "Reading field CY\n" <<endl;
volScalarField CY 
(
IOobject
(
"CY",
runTime.timeName(),
mesh,
IOobject::MUST_READ,
IOobject::AUTO_WRITE
),
mesh
);

Info<< "Reading field S\n" <<endl;
volScalarField S
(
IOobject
(
"S",
runTime.timeName(),
mesh,
IOobject::MUST_READ,
IOobject::AUTO_WRITE
),
mesh
);

Info<< "Reading field nut\n" << endl;
volScalarField nut 
(
    IOobject
    (
        "nut",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField D
(
    IOobject
    (
        "D",                
        runTime.timeName(), 
        mesh,               
        IOobject::MUST_READ,  
        IOobject::AUTO_WRITE 
    ),
    mesh
);

Info<< "Reading field C0\n" <<endl;
volScalarField C0
(
IOobject
(
"C0",
runTime.timeName(),
mesh,
IOobject::NO_READ,
IOobject::AUTO_WRITE
),
mesh
);

Info<< "Reading field C1\n" <<endl;
volScalarField C1
(
IOobject
(
"C1",
runTime.timeName(),
mesh,
IOobject::NO_READ,
IOobject::AUTO_WRITE
),
mesh
);

Info<< "Reading field DT\n" <<endl;
volScalarField DT
(
IOobject
(
"DT",
runTime.timeName(),
mesh,
IOobject::NO_READ,
IOobject::AUTO_WRITE
),
mesh
);
Info<< "Reading transportProperties\n" << endl;
IOdictionary transportProperties
(
IOobject
(
"transportProperties",
runTime.constant(),
mesh,
IOobject::MUST_READ,
IOobject::NO_WRITE
)
);


dimensionedScalar D0
(
transportProperties.lookup("D0")
);
dimensionedScalar kappa0 
(
transportProperties.lookup("kappa0")
);

dimensionedScalar kappa1 
(
transportProperties.lookup("kappa1")
);

dimensionedScalar sigma0 
(
transportProperties.lookup("sigma0")
);

dimensionedScalar sigma1 
(
transportProperties.lookup("sigma1")
);

dimensionedScalar mu0 
(
transportProperties.lookup("mu0")
);

dimensionedScalar mu1 
(
transportProperties.lookup("mu1")
);

dimensionedScalar cy0 
(
transportProperties.lookup("cy0")
);

dimensionedScalar cye 
(
transportProperties.lookup("cye")
);

#include "createPhi.H"

#include "createFvModels.H"
#include "createFvConstraints.H"
