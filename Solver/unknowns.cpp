#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "solver.hpp"
#include <cmath>
#include <vector>
#include <fstream>

Unknown::Unknown(double totalGaussPointsNumber, 
                 double totalNumNodes, 
                 std::string fluxName, 
                 std::string numFluxName, 
                 const std::vector<double> & parameters)
{
    std::size_t i;
    m_nodeValue.resize(totalNumNodes, 0.);

    m_fluxX.resize(totalNumNodes, 0.);
    m_fluxY.resize(totalNumNodes, 0.);
    m_fluxZ.resize(totalNumNodes, 0.);

    m_gaussValue.resize(totalGaussPointsNumber, std::make_pair(0,0));

    m_numFluxX.resize(totalGaussPointsNumber, 0.);
    m_numFluxY.resize(totalGaussPointsNumber, 0.);
    m_numFluxZ.resize(totalGaussPointsNumber, 0.);

    m_fluxName = fluxName;
    m_numFluxName = numFluxName;
    m_time = 0;
    m_parameters = parameters;

}

void Unknown::getFluxes(const std::vector<double> & normals,
                        const std::vector<std::pair<int, int>> & neighbours,
                        int gaussPointsNumbers)
{
    std::size_t i;

    /***********************************************************
     * 
     *          ADVECTION ZONE
     * 
     * ********************************************************/

    std::fill(m_fluxX.begin(), m_fluxX.end(), 0);

    std::fill(m_fluxY.begin(), m_fluxY.end(), 0);

    std::fill(m_fluxZ.begin(), m_fluxZ.end(), 0);

    std::fill(m_numFluxX.begin(), m_numFluxX.end(), 0);

    std::fill(m_numFluxY.begin(), m_numFluxY.end(), 0);

    std::fill(m_numFluxZ.begin(), m_numFluxZ.end(), 0);
    
    if(m_fluxName.find("Advection") != std::string::npos)
    {

        for(i = 0; i < m_fluxX.size(); ++i)
        {

            m_fluxX[i] = m_parameters[0] * m_nodeValue[i];

            m_fluxY[i] = m_parameters[1] * m_nodeValue[i];

            m_fluxZ[i] = m_parameters[2] * m_nodeValue[i];

        }

    }

    if(m_numFluxName.find("AdvectionUpwind") != std::string::npos)
    {

        std::vector<double> scalarProducts(m_gaussValue.size(), 0); // retains the  scalars products between the velocity vector and teh normal

        for(i = 0; i < normals.size(); ++i)
        {

            scalarProducts[i/3] += m_parameters[i % 3] * normals[i];

        }

        for(i = 0; i < m_gaussValue.size(); ++i)
        {

            if(scalarProducts[i] >= 0 && neighbours[i/gaussPointsNumbers].second >= 0)
            {

                m_numFluxX[i] = m_parameters[0] 
                              * normals[i * 3] 
                              * m_gaussValue[i].first;
                              
                m_numFluxY[i] = m_parameters[1] 
                              * normals[i * 3 + 1] 
                              * m_gaussValue[i].first;
                               
                m_numFluxZ[i] = m_parameters[2] 
                              * normals[i * 3 + 2] 
                              * m_gaussValue[i].first;
            }

            else if(scalarProducts[i] < 0 && neighbours[i/gaussPointsNumbers].second >= 0)
            {

                m_numFluxX[i] = m_parameters[0] 
                              * normals[i * 3] 
                              * m_gaussValue[i].second;
                              
                m_numFluxY[i] = m_parameters[1] 
                              * normals[i * 3 + 1] 
                              * m_gaussValue[i].second;
                              
                m_numFluxZ[i] = m_parameters[2] 
                              * normals[i * 3 + 2] 
                              * m_gaussValue[i].second;
                              
            }
            
            else if(m_boundaryType[i].find("Sinusoidal") != std::string::npos)
            {

                m_numFluxX[i] = m_parameters[0] 
                              * normals[i * 3] 
                              * m_parameters[3] 
                              * sin(m_parameters[4] * m_time);
                              
                m_numFluxY[i] = m_parameters[1] 
                              * normals[i * 3 + 1] 
                              * m_parameters[3] 
                              * sin(m_parameters[4] * m_time);
                              
                m_numFluxZ[i] = m_parameters[2] 
                              * normals[i * 3 + 2] 
                              * m_parameters[3] 
                              * sin(m_parameters[4] * m_time);
                              
            }

            else if(m_boundaryType[i].find("Opening") != std::string::npos)
            {

                m_numFluxX[i] = m_parameters[0] 
                              * normals[i * 3] 
                              * m_gaussValue[i].first;
                              
                m_numFluxY[i] = m_parameters[1] 
                              * normals[i * 3 + 1] 
                              * m_gaussValue[i].first;
                              
                m_numFluxZ[i] = m_parameters[2] 
                              * normals[i * 3 + 2] 
                              * m_gaussValue[i].first;
                              
            }

            else if(m_boundaryType[i].find("Wall") != std::string::npos)
            {

                m_numFluxX[i] = 0.;

                m_numFluxY[i] = 0.;

                m_numFluxZ[i] = 0.;

            } 

        }
    
    }

}

void Unknown::getBoundaryConditions(const std::vector<std::size_t> & frontierNodes, 
                                    int numGpPerFrontier,
                                    int numNodesPerFrontier,
                                    int elementType)
{

    std::size_t i, j, k; // Loops variables

    std::vector<std::pair<int,int>> physicalDimTags; //Vector pairs containing the (dimension,tag) of the physical groups.

    std::vector<std::size_t> physicalNodes; // Vector containing the list of physical nodes.

    m_boundaryType.resize(m_gaussValue.size(), "NONE");

    gmsh::model::getPhysicalGroups(physicalDimTags,
                                   gmsh::model::getDimension() - 1);

    for(i = 0; i < physicalDimTags.size(); ++i)
    {

        std::string name; // Contains the name of the physical group currently treated.

        std::vector<double> bin; // Bin variable, stores useless quantities.

        int nodeCount = 0; // Count the number of nodes of an element corresponding to the treated physical group.

        gmsh::model::getPhysicalName(physicalDimTags[i].first,
                                     physicalDimTags[i].second,
                                     name);

        gmsh::model::mesh::getNodesForPhysicalGroup(physicalDimTags[i].first,
                                                    physicalDimTags[i].second,
                                                    physicalNodes,
                                                    bin);

        for(j = 0; j < frontierNodes.size(); ++j)
        {
            if(!(j % numNodesPerFrontier) && j)
            {
                nodeCount = 0;
            }

            for(k = 0; k < physicalNodes.size(); ++k)
            {
                if(frontierNodes[j] == physicalNodes[k])
                {

                    ++nodeCount;

                    break;
                }
            }

            if(nodeCount == numNodesPerFrontier)
            {

                for(k = 0; k < numGpPerFrontier; ++k)
                {

                    m_boundaryType[j/numNodesPerFrontier * numGpPerFrontier + k] = name;

                }

            }

        }
        
    }

}

void Unknown::getGaussPointValues(const std::vector<std::pair<int, int>> & nodeCorrespondance, 
                                  int elementNumNodes, 
                                  const std::vector<double> & basisfunctions, 
                                  int numGpPerFrontier)
{
    std::size_t i, j;

    std::fill(m_gaussValue.begin(), m_gaussValue.end(), std::make_pair(0,0));

    for(i = 0; i < m_gaussValue.size(); ++i)
    {

        for(j = 0; j < elementNumNodes; ++j)
        {

            int frontierNodeIndex = i/numGpPerFrontier * elementNumNodes + j; // Contains the index of the frontier node currently treated

            int basisFunctionIndex = (i % numGpPerFrontier) * elementNumNodes + j; // Contains the index of the basis function that is treated

            m_gaussValue[i].first += m_nodeValue[nodeCorrespondance[frontierNodeIndex].first] 
                                   * basisfunctions[basisFunctionIndex];
                                    
            if(nodeCorrespondance[frontierNodeIndex].second >= 0)
            {

                m_gaussValue[i].second += m_nodeValue[nodeCorrespondance[frontierNodeIndex].second] 
                                        * basisfunctions[basisFunctionIndex];
                                        
            }

        }

    }


}