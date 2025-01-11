/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

Application
    scalarTransportFoam

Description
    Solves the steady or transient transport equation for a passive scalar.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

    #include "CourantNo.H"

    while (simple.loop(runTime))
    {
        Info<< "Time = " << runTime.userTimeName() << nl << endl;

        fvModels.correct();
		DT = nut;
		D = D0 + DT;
		C0 = (kappa0 /(sqrt(3.14 * pow(sigma0, 2)))) * exp(-pow((CY - mu0), 2) / (2*pow(sigma0, 2)));
        	C1 = - kappa1 * pow((mu1 *  mag(symm(fvc::grad(U)))), sigma1 ) *  (1 - exp(-pow((CY/cy0), cye)));
		S = C0 + C1;
		Info << "MEAN S: " << sum(S * mesh.V())/sum(mesh.V()) << endl;
		Info << "MEAN STRAIN: " << sum(mag(symm(fvc::grad(U))) * mesh.V())/sum(mesh.V()) << endl;
        while (simple.correctNonOrthogonal())
        {
        

            fvScalarMatrix CYEqn
            (
                fvm::ddt(CY)
              + fvm::div(phi, CY)
              - fvm::laplacian(D, CY)
             ==
                S
            );

            CYEqn.relax();
            fvConstraints.constrain(CYEqn);
            CYEqn.solve();
            fvConstraints.constrain(CY);
        }

        forAll(CY, celli)
		{
		    if (CY[celli] < 0.0)
		    {
			CY[celli] = 0.0;
		    }
		}

	dimensionSet customDim(-1, 1, 0, 0, 0, 0, 0);
	dimensionedScalar cy_mean = sum(CY * mesh.V()) / sum(mesh.V());
	Info << sum(mesh.V());
	Info << cy_mean << endl;

	std::ofstream outFile("scalar_output.txt", std::ios::app);
	if (outFile.is_open())
	{
    	    // Write the time and the scalar value (with unit)
    	    outFile << runTime.timeName() << " "  // Time of the simulation step
            << cy_mean.value() << "\n";  // Value of the scalar

    	    outFile.close();  // Always close the file after writing
	}
	else
	{
    	    // Handle file open error
    	   std::cerr << "Error: Unable to open the file for writing!" << std::endl;
	}

        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
