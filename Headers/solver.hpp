#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <cstdio>
#include <vector>
#include <string>
#include <iostream>

// Class containing the unknowns of the problem.
class Unknown{

    protected:
    
    std::vector<double> m_nodeValue; // Nodal values of the unknown.
    std::vector<std::pair<double,double>> m_gaussValue; // Gauss point values of the unknown.
    std::vector<double> m_fluxX; // Nodal values of the physiqcal fluxes along X.
    std::vector<double> m_fluxY; // Nodal values of the physiqcal fluxes along Y.
    std::vector<double> m_fluxZ; // Nodal values of the physiqcal fluxes along Z.
    std::vector<double> m_numFluxX; // Numerical flux value at the gauss points along X.
    std::vector<double> m_numFluxY; // Numerical flux value at the gauss points along Y.
    std::vector<double> m_numFluxZ; // Numerical flux value at the gauss points along Z.
    std::vector<double> m_parameters; // Contains the parameters for the computation of the different fluxes.
    std::string m_fluxName; // Name of the physical flux.
    std::string m_numFluxName; // flux name
    double m_time; // time at which the unknown is evaluated.
    std::vector<std::string> m_boundaryType; // Contains the boundary type of each node.

    public:

    // Initialize the unknown vector with the values at time 0 everywhere on the plane.
    Unknown(double totalGaussPointsNumber, 
            double totalNumNodes, 
            std::string fluxName, 
            std::string numFluxName, 
            std::vector<double> parameters);

    void getFluxes(const std::vector<double> normals,
                   const std::vector<std::pair<int, int>> neighbours,
                   int gaussPointsNumbers,
                   std::vector<std::vector<double>> allNodeValues,
                   std::vector<std::vector<std::pair<double,double>>> allGaussValues);

    void getGaussPointValues(std::vector<std::pair<int, int>> nodeCorrespondance, 
                             int elementNumNodes, 
                             std::vector<double> basisfunctions, 
                             int numGpPerFrontier);

    void getBoundaryConditions(std::vector<std::size_t> frontierNodes, 
                               int numGpPerFrontier,
                               int numNodesPerFrontier,
                               int elementType);

    void computeNextStep(std::vector<double> S,
                         std::vector<double> F, std::vector<double> M,
                         int numNodes,
                         double timeStep);

    void visualize(int viewTag);

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

    void incrementTime(double increment)
    {
        m_time += increment;
    }


};

void parametersLoading(const std::string paramPath,
                       std::string & mshPath,
                       int & numUnknown,
                       std::vector<std::string> & fluxName,
                       std::vector<std::string> & numFluxName,
                       std::vector<std::vector<double>> & parameters,
                       double & timeStep,
                       double & timeMax,
                       std::string & gaussType);

#endif