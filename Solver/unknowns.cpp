#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include "solver.hpp"
#include <cmath>
#include <vector>
#include <fstream>

void parametersLoading(const std::string paramPath,
                       std::string & mshPath,
                       int & numUnknown,
                       std::vector<std::string> & fluxName,
                       std::vector<std::string> & numFluxName,
                       std::vector<std::vector<double>> & parameters)
{

    std::size_t i, j;
    int numParams;
    std::ifstream parametersFile;

    parametersFile.open(paramPath);

    // Getting the name of the msh file.
    getline(parametersFile, mshPath);

    parametersFile >> numUnknown;
    parametersFile.get();

    fluxName.resize(numUnknown);
    numFluxName.resize(numUnknown);

    for(i = 0; i < numUnknown; ++i)
    {
        getline(parametersFile, fluxName[i]);
    }

    for(i = 0; i < numUnknown; ++i)
    {
        getline(parametersFile, numFluxName[i]);
    }

    parametersFile >> numParams;
    parametersFile.get();

    parameters.resize(numUnknown);

    for(i = 0; i < numUnknown; ++i)
    {
        parameters[i].resize(numParams);
        for(j = 0; j < numParams; ++j)
        {
            parametersFile >> parameters[i][j];
            parametersFile.get();
        }

    }

    parametersFile.close();

}

void Unknown::getFluxes(const std::vector<double> normals)
{
    std::size_t i;

    /***********************************************************
     * 
     *          ADVECTION ZONE
     * 
     * ********************************************************/
    
    if(m_fluxName.find("Advection") != std::string::npos)
    {

        m_fluxX.resize(m_nodeValue.size(), 0);
        m_fluxY.resize(m_nodeValue.size(), 0);
        m_fluxZ.resize(m_nodeValue.size(), 0);

        for(i = 0; i < m_fluxX.size(); ++i)
        {
            m_fluxX[i] = m_parameters[0] * m_nodeValue[i];
        }

        for(i = 0; i < m_fluxY.size() && gmsh::model::getDimension() >= 2; ++i)
        {
            m_fluxY[i] = m_parameters[1] * m_nodeValue[i];
        }

        for(i = 0; i < m_fluxY.size() && gmsh::model::getDimension() == 3; ++i)
        {
            m_fluxZ[i] = m_parameters[2] * m_nodeValue[i];
        }

    }

    // Upwind numerical flux for 1D advection

    if(m_numFluxName.find("AdvectionUpwind") != std::string::npos)
    {
        std::vector<double> scalarProducts(m_gaussValue.size(), 0); // retains the  scalars products between the velocity vector and teh normal

        m_numFluxX.resize(m_nodeValue.size(), 0);
        m_numFluxY.resize(m_nodeValue.size(), 0);
        m_numFluxZ.resize(m_nodeValue.size(), 0);

        if(gmsh::model::getDimension() == 1)
        {

        }

        else
        {
            for(i = 0; i < normals.size(); ++i)
            {
                scalarProducts[i/3] += m_parameters[i % 3] * normals[i];
            }

            for(i = 0; i < m_numFluxX.size(); ++i)
            {
                if(scalarProducts[i] > 0)
                {
                    m_numFluxX[i] = m_parameters[0] * m_gaussValue[i].first;
                    m_numFluxY[i] = m_parameters[1] * m_gaussValue[i].first;
                    m_numFluxZ[i] = m_parameters[2] * m_gaussValue[i].first;
                }

                else
                {
                    m_numFluxX[i] = m_parameters[0] * m_gaussValue[i].second;
                    m_numFluxY[i] = m_parameters[1] * m_gaussValue[i].second;
                    m_numFluxZ[i] = m_parameters[2] * m_gaussValue[i].second;
                }
            }
        }
    }

}

void Unknown::getBoundaryConditions(int tag, std::vector<std::pair<std::size_t,std::size_t>> neighbours, int numGpPerFrontier)
{
    std::size_t i, j, k, l;
    std::vector<int> physicalTag;
    m_boundaryType.resize(m_gaussValue.size(), "NONE");

    gmsh::model::getPhysicalGroupsForEntity(gmsh::model::getDimension() - 1, 
                                            tag, 
                                            physicalTag);

    for(i = 0; i < physicalTag.size(); ++i)
    {
        std::string name;
        std::vector<double> bin;
        std::vector<int> bin2;
        std::vector<std::vector<std::size_t>> elementTags, bin1;

        gmsh::model::getPhysicalName(gmsh::model::getDimension() - 1,
                                     physicalTag[i],
                                     name);

        gmsh::model::mesh::getElements(bin2,
                                       elementTags,
                                       bin1,
                                       gmsh::model::getDimension() - 1,
                                       physicalTag[i]);

        for(j = 0; j < neighbours.size(); ++j)
        {
            for(k = 0; k < elementTags[0].size() && neighbours[j].second == 0 && elementTags[0][k] != neighbours[j].first; ++k)
            {
                
                for(l = 0; l < numGpPerFrontier && elementTags[0][k] == neighbours[j].first; ++l)
                {
                    m_boundaryType[j * numGpPerFrontier + l] = name;
                } 

            }
        }
        
    }

}

void Unknown::getGaussPointValues(std::vector<std::pair<int, int>> nodeCorrespondance, 
                                  int elementNumNodes, 
                                  std::vector<double> basisfunctions, 
                                  int numGpPerFrontier)
{
    std::size_t i, j;

    std::fill(m_gaussValue.begin(), m_gaussValue.end(), std::make_pair(0,0));

    for(i = 0; i < m_gaussValue.size(); ++i)
    {
        for(j = 0; j < elementNumNodes; ++j)
        {
            m_gaussValue[i].first += m_nodeValue[nodeCorrespondance[i/numGpPerFrontier * elementNumNodes + j].first] 
                                   * basisfunctions[j * numGpPerFrontier + i % numGpPerFrontier];
                                    
            if(nodeCorrespondance[i/numGpPerFrontier * elementNumNodes + j].second >= 0)
            {
                m_gaussValue[i].second += m_nodeValue[nodeCorrespondance[i/numGpPerFrontier * elementNumNodes + j].second] 
                                        * basisfunctions[j * numGpPerFrontier + i % numGpPerFrontier];
            }

        }
    }
}