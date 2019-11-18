#include <cstdio>
#include <iostream>

#include "gmsh.h"
#include "mesh.hpp"
#include "solver.hpp"

// Main of the program.

int main(int argc, char **argv)
{
    std::string mshPath; // path to the msh file
    int numUnknown; // Number of unknown of the problem.
    std::vector<std::string> fluxName; // All name of the fluxes in the unknown order.
    std::vector<std::string> numFluxName; // All name of the numerical fluxes.
    std::vector<std::vector<double>> parameters; // Various parameters.
    double timeStep; // timeStep of the simulation
    double timeMax; // maximum time of the simulation.
    double t; // Current simulation time
    std::string gaussType; // Type of Gauss integration.
    int viewTag; // tag of the view.
    std::vector<std::string> modelNames;


    std::vector<int> elementType;
    std::vector<std::vector<std::size_t>> bin, bin1;
    gmsh::vectorpair dimTags;

    std::size_t i, j, k, l, m;

    parametersLoading(argv[1],
                      mshPath,
                      numUnknown,
                      fluxName,
                      numFluxName,
                      parameters,
                      timeStep,
                      timeMax,
                      gaussType);


    gmsh::initialize(argc, argv);
    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::option::setNumber("Mesh.SaveAll", 1);
    gmsh::open(mshPath);


    std::vector<Unknown *> u(numUnknown); // Value at the current time
    std::vector<std::vector<double>> allNodesValues(numUnknown); // Vector that contains the node values of all the unknwons. Very useful for coupled systems
    std::vector<std::vector<std::pair<double,double>>> allGaussValues(numUnknown); // smae principle for Gauss.

    // Gets all elements in the mesh.
    gmsh::model::mesh::getElements(elementType,
                                   bin,
                                   bin1,
                                   gmsh::model::getDimension());


    gmsh::model::getEntities(dimTags,
                             gmsh::model::getDimension());

    if(elementType.size() > 1)
    {
        gmsh::logger::write("No hybrid mesh can be taken into account.");
    }

    std::cout << elementType.size() << std::endl;

    Element mainElements(gaussType, elementType[0], dimTags[0].second);

    mainElements.frontierAndNeighbouring(gaussType);

    for(i = 0; i < u.size(); ++i)
    {
        u[i] = new Unknown (mainElements.getFrontierElements()->getTotalGaussPointNumber(),
                            mainElements.getTotalNumNodes(),
                            fluxName[i],
                            numFluxName[i],
                            parameters[i]);

        
        u[i]->getBoundaryConditions(mainElements.getFrontierElements()->getNodeTags(),
                                    mainElements.getFrontierElements()->getGaussPointNumber(),
                                    mainElements.getFrontierElements()->getNumNodes(),
                                    mainElements.getFrontierElements()->getType());
        

        
    }

    for(i = 0; i < numUnknown; ++i)
    {
        allNodesValues[i].resize(mainElements.getTotalNumNodes(),0);
        allGaussValues[i].resize(mainElements.getFrontierElements()->getTotalGaussPointNumber(), std::make_pair (0,0));
    }

    gmsh::model::list(modelNames);
    viewTag = gmsh::view::add("View");

    for(t = 0; t < timeMax; t += timeStep)
    {
        for(i = 0; i < u.size(); ++i)
        {
            // Acquisition of the Gauss points values.
            u[i]->getGaussPointValues(mainElements.getCorrespondance(),
                                      mainElements.getFrontierElements()->getNumNodes(),
                                      mainElements.getFrontierElements()->getBasisFunctions(),
                                      mainElements.getFrontierElements()->getGaussPointNumber());
                                      
            // Acquisition of the fluxes.
            u[i]->getFluxes(mainElements.getFrontierElements()->normals(), 
                            mainElements.getNeighbours(),
                            mainElements.getFrontierElements()->getGaussPointNumber(),
                            allNodesValues,
                            allGaussValues);

            // Multiplication of the fluxes and the stiffness matrices
            std::vector<double> S(mainElements.getTotalNumNodes(), 0);

            for(j = 0; j < S.size()/(mainElements.getNumNodes()); ++j)
            {
                for(k = 0; k < mainElements.getNumNodes(); ++k)
                {
                    for(l = 0; l < mainElements.getNumNodes(); ++l)
                    {
                        S[j * mainElements.getNumNodes() + k] += (mainElements.getStiffnessMatrixX(j * mainElements.getNumNodes() * mainElements.getNumNodes() 
                                                                                                  + k * mainElements.getNumNodes() 
                                                                                                  + l)
                                                               * u[i]->fluxX(j * mainElements.getNumNodes() + l)
                                                               + mainElements.getStiffnessMatrixY(j * mainElements.getNumNodes() * mainElements.getNumNodes() 
                                                                                                  + k * mainElements.getNumNodes() 
                                                                                                  + l)
                                                               * u[i]->fluxY(j * mainElements.getNumNodes() + l)
                                                               + mainElements.getStiffnessMatrixZ(j * mainElements.getNumNodes() * mainElements.getNumNodes() 
                                                                                                  + k * mainElements.getNumNodes() 
                                                                                                  + l)
                                                               * u[i]->fluxZ(j * mainElements.getNumNodes() + l));
                    }
                    
                }

            }
            

            // Computation of the matrix F's, containing the numerical fluxes.
            std::vector<double> F(mainElements.getTotalNumNodes(), 0);
            std::vector<double> frontierBasisFunctions = mainElements.getFrontierElements()->getBasisFunctions();
            std::vector<double> frontierGaussWeight = mainElements.getFrontierElements()->getGaussWeight();
            std::vector<double> frontierJacobians = mainElements.getFrontierElements()->getJacobians();
            std::vector<std::pair<int,int>> nodeCorrespondance = mainElements.getCorrespondance();

            for(j = 0; j < mainElements.getFrontierElements()->getTotalNumNodes()/mainElements.getFrontierElements()->getNumNodes(); ++j)
            {
                for(k = 0; k < mainElements.getFrontierElements()->getNumNodes(); ++k)
                {
                    for(l = 0; l < mainElements.getFrontierElements()->getGaussPointNumber(); ++l)
                    {
                        double tmp = frontierBasisFunctions[l * mainElements.getFrontierElements()->getNumNodes() + k]
                                   * frontierGaussWeight[l]
                                   * frontierJacobians[j * mainElements.getFrontierElements()->getGaussPointNumber() + l]
                                   *  ( u[i]->numFluxX(j * mainElements.getFrontierElements()->getGaussPointNumber() + l)
                                      + u[i]->numFluxY(j * mainElements.getFrontierElements()->getGaussPointNumber() + l)
                                      + u[i]->numFluxZ(j * mainElements.getFrontierElements()->getGaussPointNumber() + l));

                        F[nodeCorrespondance[j * mainElements.getFrontierElements()->getNumNodes() + k].first] -= tmp;
                            
                        if(nodeCorrespondance[j * mainElements.getFrontierElements()->getNumNodes() + k].second >= 0)
                        {
                            F[nodeCorrespondance[j * mainElements.getFrontierElements()->getNumNodes() + k].second] += tmp;
                        }
                    }

                }

            }
            

            u[i]->computeNextStep(S, 
                                  F, 
                                  mainElements.getMassMatrixInverse(), 
                                  mainElements.getNumNodes(), 
                                  timeStep);

            

            u[i]->incrementTime(timeStep);
            

        }


        for(i = 0; i < numUnknown; ++i)
        {
            allGaussValues[i] = u[i]->gaussPointValues();
            allNodesValues[i] = u[i]->nodeValues();
        }
        

        std::vector<std::vector<double>> showNodeValues(mainElements.getTotalNumNodes());

        for(i = 0; i < showNodeValues.size(); ++i)
        {
            showNodeValues[i].resize(numUnknown);

            for(j = 0; j < showNodeValues[i].size(); ++j)
            {
                showNodeValues[i][j] = allNodesValues[j][i];
            }
        }
        

        gmsh::view::addModelData(viewTag, 
                                 int(t/timeStep), 
                                 modelNames[0], 
                                 "NodeData" , 
                                 mainElements.getNodeTags(), 
                                 showNodeValues, 
                                 t, 
                                 numUnknown);
                                 

        gmsh::view::write(viewTag, "testView.msh");
        
        
    }

    for(i = 0; i < u.size(); ++i)
    {
        if( u[i] != NULL)
        {
            delete u[i];
        }
    }

    gmsh::finalize();

    return 0;
}