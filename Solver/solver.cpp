#include <cstdio>
#include <iostream>

#include "gmsh.h"
#include "mesh.hpp"
#include "solver.hpp"

std::vector<std::vector<double>> computeRightHandSide(std::vector<Unknown *> & u, 
                                                      Element & mainElements,
                                                      std::vector<std::vector<double>> & allNodesValues,
                                                      std::vector<std::vector<std::pair<double, double>>> & allGaussValues)
{

    std::size_t i, j, k, l;

    std::vector<std::vector<double>> rightHandSide(u.size());

    for(i = 0; i < u.size(); ++i)
    {
        rightHandSide[i].resize(mainElements.getTotalNumNodes(), 0);
    }

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

        std::vector<double> F(mainElements.getTotalNumNodes(), 0);

        std::vector<double> M = mainElements.getMassMatrixInverse();

        std::vector<double> frontierBasisFunctions = mainElements.getFrontierElements()->getBasisFunctions();

        std::vector<double> frontierGaussWeight = mainElements.getFrontierElements()->getGaussWeight();

        std::vector<double> frontierJacobians = mainElements.getFrontierElements()->getJacobians();

        std::vector<std::pair<int,int>> nodeCorrespondance = mainElements.getCorrespondance();

        int frontierGpNumber = mainElements.getFrontierElements()->getGaussPointNumber();

        int frontierNumNodes = mainElements.getFrontierElements()->getNumNodes();

        int mainNumNodes = mainElements.getNumNodes();

        for(j = 0; j < S.size()/mainNumNodes; ++j)
        {

            for(k = 0; k < mainNumNodes; ++k)
            {

                for(l = 0; l < mainNumNodes; ++l)
                {

                    int stiffIndex = j * mainNumNodes * mainNumNodes + k * mainNumNodes + l;

                    int fluxIndex = j * mainNumNodes + l;

                    S[j * mainNumNodes + k] += mainElements.getStiffnessMatrixX(stiffIndex)
                                                * u[i]->fluxX(fluxIndex)
                                                + mainElements.getStiffnessMatrixY(stiffIndex)
                                                * u[i]->fluxY(fluxIndex)
                                                + mainElements.getStiffnessMatrixZ(stiffIndex)
                                                * u[i]->fluxZ(fluxIndex);
                                                
                }
                
            }

        }

        for(j = 0; j < mainElements.getFrontierElements()->getTotalNumNodes()/mainElements.getFrontierElements()->getNumNodes(); ++j)
        {

            for(k = 0; k < frontierNumNodes; ++k)
            {
                double tmp = 0;

                int nodeFrontIndex = j * frontierNumNodes + k;

                for(l = 0; l < frontierGpNumber; ++l)
                {

                    int gpFrontIndex = j * frontierGpNumber + l;

                    tmp += frontierBasisFunctions[l * frontierNumNodes + k]
                        * frontierGaussWeight[l]
                        * frontierJacobians[gpFrontIndex]
                        * ( u[i]->numFluxX(gpFrontIndex)
                            + u[i]->numFluxY(gpFrontIndex)
                            + u[i]->numFluxZ(gpFrontIndex));

                }

                F[nodeCorrespondance[nodeFrontIndex].first] -= tmp;
                        
                if(nodeCorrespondance[nodeFrontIndex].second >= 0)
                {

                    F[nodeCorrespondance[nodeFrontIndex].second] += tmp;

                }

            }

        }

        for(j = 0; j < M.size()/(mainNumNodes * mainNumNodes); ++j)
        {

            for(k = 0; k < mainNumNodes; ++k)
            {

                for(l = 0; l < mainNumNodes; ++l)
                {
                    rightHandSide[i][j * mainNumNodes + k] += M[j * mainNumNodes * mainNumNodes + k * mainNumNodes + l] 
                                                            * (S[j * mainNumNodes + l] + F[j * mainNumNodes + l]);
                }

            }

        }

    }

    return rightHandSide;

}

void solver(int numUnknown,
            Element & mainElements,
            std::vector<std::string> fluxName,
            std::vector<std::string> numFluxName,
            std::vector<std::vector<double>> parameters, 
            double timeStep,
            double timeMax, 
            std::vector<std::string> & viewNames,
            std::string timeMethod)
{

    double t;

    std::size_t i, j;

    std::vector<Unknown *> u(numUnknown);

    std::vector<std::vector<double>> showResults(mainElements.getTotalNumNodes());

    std::vector<std::vector<double>> allNodesValues(u.size());

    std::vector<std::vector<std::pair<double,double>>> allGaussValues(u.size());

    std::vector<std::string> modelNames; // Get the model name

    std::vector<int> viewTags(u.size());

    gmsh::model::list(modelNames);

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

    for(i = 0; i < showResults.size(); ++i)
    {

        showResults[i].resize(1);
        
    }

    for(i = 0; i < viewTags.size(); ++i)
    {

        viewTags[i] = gmsh::view::add(viewNames[i]);

    }


    for(i = 0; i < u.size(); ++i)
    {

        allNodesValues[i].resize(mainElements.getTotalNumNodes(),0);

        allGaussValues[i].resize(mainElements.getFrontierElements()->getTotalGaussPointNumber(), std::make_pair (0,0));
        
    }

    for(t = 0; t < timeMax; t += timeStep)
    {


        if(timeMethod.find("Euler") != std::string::npos)
        {
            
            std::vector<std::vector<double>> rightHandSide = computeRightHandSide(u,
                                                                                  mainElements,
                                                                                  allNodesValues,
                                                                                  allGaussValues);

            for(i = 0; i < u.size(); ++i)
            {
                for(j = 0; j < rightHandSide[i].size(); ++j)
                {

                    rightHandSide[i][j] *= timeStep;

                }

                u[i]->update(rightHandSide[i]);

            }

        }

        else if(timeMethod.find("Runge Kutta 4") != std::string::npos)
        {
            

        }


        for(i = 0; i < allGaussValues.size(); ++i)
        {

            allGaussValues[i] = u[i]->gaussPointValues();

            allNodesValues[i] = u[i]->nodeValues();

        }

        for(i = 0; i < u.size(); ++i)
        {

            for(j = 0; j < showResults.size(); ++j)
            {

                showResults[j][0] = u[i]->nodeValue(j);
                
            }

            gmsh::view::addModelData(viewTags[i], 
                                    int(t/timeStep), 
                                    modelNames[0], 
                                    "NodeData" , 
                                    mainElements.getNodeTags(), 
                                    showResults, 
                                    t);

            gmsh::view::write(viewTags[i], viewNames[i]);

        }

        for(i = 0; i < u.size(); ++i)
        {

            u[i]->incrementTime(timeStep);
            
        }

    }

    for(i = 0; i < u.size(); ++i)
    {
        if(u[i])
        {
            delete u[i];
        }
    }

}