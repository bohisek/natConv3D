/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2015 FlowKit Sarl
 * Route d'Oron 2
 * 1010 Lausanne, Switzerland
 * E-mail contact: contact@flowkit.com
 *
 * The most recent release of Palabos can be downloaded at 
 * <http://www.palabos.org/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/* Main author: Orestis Malaspinas
 */

/** \file
 * A fluid constrained between a hot bottom wall (no-slip for the velocity) and a cold
 * top wall (no-slip for the velocity). The lateral walls are periodic. Under the
 * influence of gravity, convection rolls are formed. Thermal effects are modelled
 * by means of a Boussinesq approximation: the fluid is incompressible, and the influence
 * of the temperature is visible only through a body-force term, representing buoyancy
 * effects. The temperature field obeys an advection-diffusion equation.
 *
 * The simulation is first created in a fully symmetric manner. The symmetry is therefore
 * not spontaneously broken; while the temperature drops linearly between the hot and
 * and cold wall, the convection rolls fail to appear at this point. In a second stage, a
 * random noise is added to trigger the instability.
 *
 * This application is technically a bit more advanced than the other ones, because it
 * illustrates the concept of data processors. In the present case, they are used to
 * create the initial condition, and to trigger the instability.
 **/

#include "palabos3D.h"
#include "palabos3D.hh"

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>

using namespace plb;
using namespace std;

typedef double T;

#define NSDESCRIPTOR descriptors::ForcedD3Q19Descriptor
#define ADESCRIPTOR descriptors::AdvectionDiffusionD3Q7Descriptor

#define ADYNAMICS AdvectionDiffusionBGKdynamics
#define NSDYNAMICS GuoExternalForceBGKdynamics

/// Fixing temperature field.
template<typename T, template<typename NSU> class nsDescriptor,
                     template<typename ADU> class adDescriptor>
struct FixTempBoolMaskProcessor3D :
    public BoxProcessingFunctional3D_LS<T,adDescriptor,bool> 
{
    FixTempBoolMaskProcessor3D(RayleighBenardFlowParam<T,nsDescriptor,adDescriptor> parameters_)
        : parameters(parameters_)
    { }
    virtual void process(Box3D domain, BlockLattice3D<T,adDescriptor>& adLattice, ScalarField3D<bool>& boolMask)
    {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            	for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                	if (boolMask.get(iX,iY,iZ)) {
                        	T temperature = (T) 5.0;				// Here, set temperature of fibers
                		adLattice.get(iX,iY,iZ).defineDensity(temperature);
                		Array<T,adDescriptor<T>::d> jEq(0.0, 0.0, 0.0);
                		iniCellAtEquilibrium(adLattice.get(iX,iY,iZ), temperature, jEq);
                	}
                }         
            }
        }
    }
    virtual FixTempBoolMaskProcessor3D<T,nsDescriptor,adDescriptor>* clone() const
    {
        return new FixTempBoolMaskProcessor3D<T,nsDescriptor,adDescriptor>(*this);
    }
    
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::dynamicVariables;
    }
    
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulk;
    }
    
private :
    RayleighBenardFlowParam<T,nsDescriptor,adDescriptor> parameters;
};

//-----------------------------------------------------------------------------------------------------------

/// Initialization of the temperature field.
template<typename T, template<typename NSU> class nsDescriptor,
                     template<typename ADU> class adDescriptor>
struct IniTemperatureRayleighBenardProcessor3D :
    public BoxProcessingFunctional3D_LS<T,adDescriptor,bool> 
{
    IniTemperatureRayleighBenardProcessor3D(RayleighBenardFlowParam<T,nsDescriptor,adDescriptor> parameters_)
        : parameters(parameters_)
    { }
    virtual void process(Box3D domain, BlockLattice3D<T,adDescriptor>& adLattice, ScalarField3D<bool>& boolMask)
    {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            	for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                	               
                	T temperature = (T) 20.0;
                                 
                	if (boolMask.get(iX,iY,iZ)) {
                	temperature = (T) 5.0;
                	}
                
                	Array<T,adDescriptor<T>::d> jEq(0.0, 0.0, 0.0);
                	adLattice.get(iX,iY,iZ).defineDensity(temperature);
                	iniCellAtEquilibrium(adLattice.get(iX,iY,iZ), temperature, jEq);
                }
            }
        }
    }
    virtual IniTemperatureRayleighBenardProcessor3D<T,nsDescriptor,adDescriptor>* clone() const
    {
        return new IniTemperatureRayleighBenardProcessor3D<T,nsDescriptor,adDescriptor>(*this);
    }
    
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;
    }
    
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulkAndEnvelope;
    }
    
private :
    RayleighBenardFlowParam<T,nsDescriptor,adDescriptor> parameters;
};

//----------------------------------------------------------------------------------------------------------

void rayleighBenardSetup (
        MultiBlockLattice3D<T, NSDESCRIPTOR>& nsLattice,
        MultiBlockLattice3D<T, ADESCRIPTOR>& adLattice,
        OnLatticeBoundaryCondition3D<T,NSDESCRIPTOR>& nsBoundaryCondition,
        OnLatticeAdvectionDiffusionBoundaryCondition3D<T,ADESCRIPTOR>& adBoundaryCondition,
        MultiScalarField3D<bool>& boolMask,
        RayleighBenardFlowParam<T,NSDESCRIPTOR,ADESCRIPTOR> &parameters )
{
    nsBoundaryCondition.setVelocityConditionOnBlockBoundaries(nsLattice);  // default dirichlet (no slip)
    
    adBoundaryCondition.setTemperatureConditionOnBlockBoundaries(adLattice);
    
    initializeAtEquilibrium(nsLattice, nsLattice.getBoundingBox(), (T)1., Array<T,3>((T)0.,(T)0.,(T)0.) );
    
    applyProcessingFunctional (
            new IniTemperatureRayleighBenardProcessor3D<T,NSDESCRIPTOR,ADESCRIPTOR>(parameters), 
            adLattice.getBoundingBox(), adLattice, boolMask);
            
    setBoundaryDensity(adLattice,adLattice.getBoundingBox(), 20.0);        
            
    nsLattice.initialize();
    adLattice.initialize();
}

//-----------------------------------------------------------------------------------------------------------

void writeVTK(MultiBlockLattice3D<T,NSDESCRIPTOR>& nsLattice,
              MultiBlockLattice3D<T,ADESCRIPTOR>& adLattice,
              RayleighBenardFlowParam<T,NSDESCRIPTOR,ADESCRIPTOR> const& parameters, plint iter)
{
    T dx = parameters.getDeltaX();
    T dt = parameters.getDeltaT();

    plint nx = parameters.getNx();
    plint ny = 600;
    plint nz = parameters.getNz();

    Box3D xz_vtkDomain(0, nx-1, ny, ny, 0, nz-1);  
    
    VtkImageOutput3D<T> vtkOut1(createFileName("vtk_xz", iter, 6), dx);
    vtkOut1.writeData<float>(*computeDensity(adLattice, xz_vtkDomain), "temperature", (T)1);
    vtkOut1.writeData<float>(*computeVelocityNorm(nsLattice, xz_vtkDomain), "velocityNorm", dx/dt);
    vtkOut1.writeData<3,float>(*computeVelocity(nsLattice, xz_vtkDomain), "velocity", dx/dt);
    //vtkOut.writeData<float>(*computeDensity(nsLattice), "pressure", (T)1/3);


    nx = 600;
    ny = parameters.getNy();
    //nz = parameters.getNz();

    Box3D yz_vtkDomain(nx, nx, 0, ny-1, 0, nz-1);  
    
    VtkImageOutput3D<T> vtkOut2(createFileName("vtk_yz", iter, 6), dx);
    vtkOut2.writeData<float>(*computeDensity(adLattice, yz_vtkDomain), "temperature", (T)1);
    vtkOut2.writeData<float>(*computeVelocityNorm(nsLattice, yz_vtkDomain), "velocityNorm", dx/dt);
    vtkOut2.writeData<3,float>(*computeVelocity(nsLattice, yz_vtkDomain), "velocity", dx/dt);
}

//-----------------------------------------------------------------------------------------------------------

void writeVTK(MultiScalarField3D<bool>& boolMask,
              RayleighBenardFlowParam<T,NSDESCRIPTOR,ADESCRIPTOR> const& parameters)
{
    T dx = parameters.getDeltaX();     
    VtkImageOutput3D<bool> vtkOut("vtk_voxels", dx);
    vtkOut.writeData<bool>(boolMask, "voxels", (T)1);
}

//-----------------------------------------------------------------------------------------------------------

int main(int argc, char *argv[])
{
    plbInit(&argc, &argv);
    
    global::timer("simTime").start();
    
    T Ra=0.;
    try {
        global::argv(1).read(Ra);
    }
    catch(PlbIOException& exception) {
        pcout << exception.what() << endl;
        pcout << "The structure of the input parameters should be : "
              << (string)global::argv(0) << " Ra" << endl;;
        // Exit the program, because wrong input data is a fatal error.
        exit(1);
    }

    plint iniT = 0;   // initial time

    plb_ifstream ifile("lastTimestep.dat");
    if (ifile.is_open())
    {
        ifile >> iniT;
        global::mpi().bCast(&iniT,1);
    }

    const T lx  = 1.;        
    const T ly  = 1.;    
    const T lz  = 1.;
    const T uMax  = 0.018;    // look at "../MATLAB/for_Palabos/LBMunits.m"
    const T Pr = 0.7442;      
    
    const T hotTemperature = 20.0;
    const T coldTemperature = 5.0;
    const plint resolution = 1199;

    global::directories().setOutputDir("./tmp/");
    global::IOpolicy().activateParallelIO(false);
    
    RayleighBenardFlowParam<T,NSDESCRIPTOR,ADESCRIPTOR> parameters (
            Ra, 
            Pr, 
            uMax,
            coldTemperature,
            hotTemperature, 
            resolution, 
            lx, 
            ly,
            lz );
                                        
    writeLogFile(parameters,"palabos.log");
    
    plint nx = parameters.getNx();
    plint ny = parameters.getNy();
    plint nz = parameters.getNz();
    
    T nsOmega = parameters.getSolventOmega();
    T adOmega = parameters.getTemperatureOmega();
    
    MultiBlockLattice3D<T, NSDESCRIPTOR> nsLattice (
            nx,ny,nz,new NSDYNAMICS<T, NSDESCRIPTOR>(nsOmega) );
       
    MultiBlockLattice3D<T, ADESCRIPTOR> adLattice (
            nx,ny,nz,new ADYNAMICS<T, ADESCRIPTOR>(adOmega) );
    
            
    OnLatticeBoundaryCondition3D<T,NSDESCRIPTOR>*
        nsBoundaryCondition = createLocalBoundaryCondition3D<T,NSDESCRIPTOR>();
        
    OnLatticeAdvectionDiffusionBoundaryCondition3D<T,ADESCRIPTOR>*
        adBoundaryCondition = createLocalAdvectionDiffusionBoundaryCondition3D<T,ADESCRIPTOR>();
    
    nsLattice.toggleInternalStatistics(false);
    adLattice.toggleInternalStatistics(false);
    
    pcout << "reading the boolMask ..." << endl;
    // ----- bool mask ------
    MultiScalarField3D<bool> boolMask(parameters.getNx(), parameters.getNy(), parameters.getNz());
    plb_ifstream boolFile("voxels.dat");
    boolFile >> boolMask;
    defineDynamics(nsLattice, boolMask, new BounceBack<T,NSDESCRIPTOR>, true);

    pcout << "finished." << endl;
    
    rayleighBenardSetup(nsLattice, adLattice,*nsBoundaryCondition, *adBoundaryCondition, boolMask, parameters);
    
    plint processorLevel = 1;
    integrateProcessingFunctional(
            new FixTempBoolMaskProcessor3D<T,NSDESCRIPTOR,ADESCRIPTOR>(parameters), 
            adLattice.getBoundingBox(), adLattice, boolMask, processorLevel);
    
    Array<T,NSDESCRIPTOR<T>::d> forceOrientation(T(),T(),(T)1.);
    
    integrateProcessingFunctional (
            new BoussinesqThermalProcessor3D<T,NSDESCRIPTOR,ADESCRIPTOR> (
                parameters.getLatticeGravity(), parameters.getAverageTemperature(),
                parameters.getDeltaTemperature(),forceOrientation ),
            nsLattice.getBoundingBox(),
            nsLattice, adLattice, processorLevel );
    
    T tIni = global::timer("simTime").stop();
    pcout << "time elapsed for rayleighBenardSetup:" << tIni << endl;
    global::timer("simTime").start();
    
    plint iT = 0;
    plint maxT = 100000;
    plint saveIter = 1000;
    T tCheck = 0.;
    T tMax = 16.0*3600.;


    if (iniT > 0)
	{
	    loadBinaryBlock(nsLattice, "checkpointTemp.dat");
    	    loadBinaryBlock(adLattice, "checkpointVel.dat");
	}
    else
	{
	    writeVTK(boolMask, parameters);
	}    

    // Main loop over time iterations.
    for (iT = iniT; iT <= maxT+iniT; ++iT) 
    {  
        //pcout << "iter: " << iT << endl;

        if (iT % saveIter == 0)
        {
            pcout << "At time " << iT * parameters.getDeltaT() << " : Writing VTK." << endl;
            writeVTK(nsLattice, adLattice, parameters, iT);
        } 

        tCheck = global::timer("simTime").stop();
        global::timer("simTime").start();

        if (tCheck>tMax) break;
 
        // Lattice Boltzmann iteration step.
        adLattice.collideAndStream();
        nsLattice.collideAndStream();
    }

    pcout << "writing CheckpointTemp.dat ... " << endl;
    saveBinaryBlock(nsLattice, "checkpointTemp.dat");
    pcout << "writing CheckpointVel.dat ... " << endl;
    saveBinaryBlock(adLattice, "checkpointVel.dat");
    
    plb_ofstream ofile("lastTimestep.dat");
    ofile << iT-1 << endl;
           
    T tEnd = global::timer("simTime").stop();
    
    T totalTime = tEnd-tIni;
    pcout << "N=" << resolution << endl;
    pcout << "number of processors: " << global::mpi().getSize() << endl;
    pcout << "simulation time: " << totalTime << endl;
    pcout << "total time: " << tEnd << endl;
    pcout << "total iterations: " << iT << endl;
}
