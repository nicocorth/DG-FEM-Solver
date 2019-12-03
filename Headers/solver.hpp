#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <cstdio>
#include <vector>
#include <string>
#include <iostream>

#include "mesh.hpp"
#include "parameters.hpp"

// Class containing the unknowns of the problem.
class Unknown{

    protected:
    
    std::vector<double> m_nodeValue; // Nodal values of the unknown.
    std::vector<std::vector<double>> m_nodeProperties; // Properties at the given node. It is of the form "N1prop1, N1prop2, ..., N2prop1, ... NNpropP"
    std::vector<std::pair<double,double>> m_gaussValue; // Gauss point values of the unknown.
    std::vector<std::pair<double,double>> m_gaussProperties; // Same principle as m_nodeProperties pour les pts de gauss.
    std::vector<double> m_fluxX; // Nodal values of the physical fluxes along X.
    std::vector<double> m_fluxY; // Nodal values of the physical fluxes along Y.
    std::vector<double> m_fluxZ; // Nodal values of the physical fluxes along Z.
    std::vector<double> m_numFluxX; // Numerical flux value at the gauss points along X.
    std::vector<double> m_numFluxY; // Numerical flux value at the gauss points along Y.
    std::vector<double> m_numFluxZ; // Numerical flux value at the gauss points along Z.
    std::vector<double> m_parameters; // Contains the parameters for the computation of the different fluxes.
    std::string m_fluxName; // Name of the physical flux.
    std::string m_numFluxName; // Name of the numerical flux.
    double m_time; // Time at which the unknown is evaluated.
    std::vector<std::string> m_boundaryType; // Contains the boundary type of each node.

    public:

    /* This is the Unknown constructor. totalGaussPointsNumber represents the total number of gauss points on
    the frontier elements. totalNumNodes is the total number of nodes on all main elements. fluxNname and
    numFluxName contains the name of the physical and numerical fluxes respectively. parameters gather all
    required parameters of the simulation.*/

    Unknown(double totalGaussPointsNumber, 
            double totalNumNodes, 
            std::string fluxName, 
            std::string numFluxName, 
            const std::vector<double> & parameters);

    /* This function compute the fluxes of the Unknown. normals are the normal to the gauss point of the
    frontier elements. neighbours contains the index of the main elements which are on each side of each frontier
    element. gaussPointsNumber is the number of gauss points per frontier element. allNodeValues si a 2D vector
    whose length is the total number of unknowns of the system. Each entry contains a vector with the value
    of an unknwon at each node. allGaussValues is similar to allNodeValues but it stores the values at the 
    Gauss points. */

    void getFluxes(const std::vector<double> & normals,
                   const std::vector<std::pair<int, int>> & neighbours,
                   int gaussPointsNumbers);

    void getGaussPointValues(const std::vector<std::pair<int, int>> & nodeCorrespondance, 
                             int elementNumNodes, 
                             const std::vector<double> & basisfunctions, 
                             int numGpPerFrontier);

    void getBoundaryConditions(const std::vector<std::size_t> & frontierNodes, 
                               int numGpPerFrontier,
                               int numNodesPerFrontier,
                               int elementType);


    void getNodeProperties();

    void getGaussProperties();

    void update(const std::vector<double> & updater)
    {
        std::size_t i;

        for(i = 0; i < updater.size(); ++i)
        {

            m_nodeValue[i] += updater[i];

        }
        
    }

    std::vector<double> parameters()
    {
        return m_parameters;
    }

    std::vector<std::pair<double,double>> gaussPointValues()
    {
        return m_gaussValue;
    }

    std::vector<double> nodeValues()
    {
        return m_nodeValue;
    }

    std::string fluxName()
    {
        return m_fluxName;
    }

    std::string numFluxName()
    {
        return m_numFluxName;
    }

    double nodeValue(std::size_t i)
    {
        return m_nodeValue[i];
    }

    double fluxX(std::size_t i)
    {
        return m_fluxX[i];
    }

    double fluxY(std::size_t i)
    {
        return m_fluxY[i];
    }

    double fluxZ(std::size_t i)
    {
        return m_fluxZ[i];
    }

    double numFluxX(std::size_t i)
    {
        return m_numFluxX[i];
    }

    double numFluxY(std::size_t i)
    {
        return m_numFluxY[i];
    }

    double numFluxZ(std::size_t i)
    {
        return m_numFluxZ[i];
    }

    void setTime(double t)
    {
        m_time = t;
    }

};

void solver(Parameters & simulInfos, Element & mainElements);


#endif