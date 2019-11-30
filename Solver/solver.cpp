#include <cstdio>
#include <iostream>

#include "gmsh.h"
#include "mesh.hpp"
#include "solver.hpp"
#include "parameters.hpp"

std::vector<std::vector<double>> computeRightHandSide(std::vector<Unknown *> & u, 
                                                      Element & mainElements)
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
                        mainElements.getFrontierElements()->getGaussPointNumber());

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

void fullUpdate(std::vector<Unknown *> & u,
                std::vector<std::vector<double>> values,
                const double t,
                const double factor = 1)
{

    std::size_t i, j;

    for(i = 0; i < u.size(); ++i)
    {
        for(j = 0; j < values[i].size(); ++j)
        {

            values[i][j] *= factor;

        }

        u[i]->update(values[i]);

        u[i]->setTime(t);

    }

}

void vecMultiply(std::vector<std::vector<double>> & vec, const double multConst)
{

    std::size_t i, j;

    for(i = 0; i < vec.size(); ++i)
    {

        for(j = 0; j < vec[i].size(); ++j)
        {

            vec[i][j] *= multConst;

        }

    }

}

void solver(Parameters & simulInfos, Element & mainElements)
{

    double t;

    int count = 0;

    std::size_t i, j;

    std::vector<Unknown *> u(simulInfos.numUnknown());

    std::vector<std::vector<double>> showResults(mainElements.getTotalNumNodes());

    std::vector<std::vector<double>> rightHandSide(u.size());

    std::vector<std::string> modelNames; // Get the model name

    std::vector<int> viewTags(u.size());

    gmsh::model::list(modelNames);

    for(i = 0; i < u.size(); ++i)
    {

        u[i] = new Unknown (mainElements.getFrontierElements()->getTotalGaussPointNumber(),
                            mainElements.getTotalNumNodes(),
                            simulInfos.fluxName(i),
                            simulInfos.numFluxName(i),
                            simulInfos.parameters(i));

        
        u[i]->getBoundaryConditions(mainElements.getFrontierElements()->getNodeTags(),
                                    mainElements.getFrontierElements()->getGaussPointNumber(),
                                    mainElements.getFrontierElements()->getNumNodes(),
                                    mainElements.getFrontierElements()->getType());

        rightHandSide[i].resize(mainElements.getTotalNumNodes());
        
        
    }

    for(i = 0; i < showResults.size(); ++i)
    {

        showResults[i].resize(1);
        
    }

    for(i = 0; i < viewTags.size(); ++i)
    {

        viewTags[i] = gmsh::view::add(simulInfos.viewNames(i));

    }

    for(t = 0; t < simulInfos.timeMax(); t += simulInfos.timeStep())
    {

        if(simulInfos.timeMethod().find("Euler") != std::string::npos)
        {

            rightHandSide = computeRightHandSide(u, mainElements);

        }

        else if(simulInfos.timeMethod().find("Runge Kutta 4") != std::string::npos)
        {

            std::vector<Unknown *> tmp = u;

            std::vector<std::vector<double>> k1 = computeRightHandSide(tmp, mainElements);

            fullUpdate(tmp, k1, t + simulInfos.timeStep()/2, simulInfos.timeStep()/2);

            std::vector<std::vector<double>> k2 = computeRightHandSide(tmp, mainElements);

            tmp = u;

            fullUpdate(tmp, k2, t + simulInfos.timeStep()/2, simulInfos.timeStep()/2);

            std::vector<std::vector<double>> k3 = computeRightHandSide(tmp, mainElements);

            tmp = u;

            fullUpdate(tmp, k3, t + simulInfos.timeStep(), simulInfos.timeStep());

            std::vector<std::vector<double>> k4 = computeRightHandSide(tmp, mainElements);

            for(i = 0; i < rightHandSide.size(); ++i)
            {

                for(j = 0; j < rightHandSide[i].size(); ++j)
                {

                    rightHandSide[i][j] = 1/6 * (k1[i][j] + 2 * k2[i][j] + 2 * k3[i][j] + k4[i][j]);

                }

            }
            
        }

        fullUpdate(u, rightHandSide, t + simulInfos.timeStep(), simulInfos.timeStep());

        if(!(count % simulInfos.timeRegistration()) && count != 0)
        {

            for(i = 0; i < u.size(); ++i)
            {

                for(j = 0; j < showResults.size(); ++j)
                {

                    showResults[j][0] = u[i]->nodeValue(j);
                    
                }

                gmsh::view::addModelData(viewTags[i], 
                                        int(t/simulInfos.timeStep()), 
                                        modelNames[0], 
                                        "NodeData" , 
                                        mainElements.getNodeTags(), 
                                        showResults, 
                                        t);

                gmsh::view::write(viewTags[i], simulInfos.viewNames(i));

            }

            count = 0;

        }

        ++count;

    }

    for(i = 0; i < u.size(); ++i)
    {

        if(u[i])
        {

            delete u[i];

        }

    }

}