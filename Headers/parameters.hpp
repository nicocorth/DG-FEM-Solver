#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <iostream>
#include <gmsh.h>
#include <cmath>
#include <vector>
#include <fstream>

#include "gmsh.h"
#include "mesh.hpp"
#include "solver.hpp"

// Main of the program.

class Parameters{

    private :

    std::string m_mshPath; // path to the msh file
    int m_numUnknown; // Number of unknown of the problem.
    std::vector<std::string> m_fluxName; // All name of the fluxes in the unknown order.
    std::vector<std::string> m_numFluxName; // All name of the numerical fluxes.
    std::vector<std::vector<double>> m_parameters; // Various parameters.
    double m_timeStep; // timeStep of the simulation
    double m_timeMax; // maximum time of the simulation.
    std::string m_gaussType; // Type of Gauss integration.
    std::vector<std::string> m_viewNames; // get the name of the final views.
    std::string m_timeMethod; // Contains the name of the temporal method.
    int m_timeRegistration; // Contains the times at which registration occur.

    public:

    std::string mshPath()
    {
        return m_mshPath;
    }

    int numUnknown()
    {
        return m_numUnknown;
    }

    std::string fluxName(std::size_t i)
    {
        return m_fluxName[i];
    }

    std::string numFluxName(std::size_t i)
    {
        return m_numFluxName[i];
    }

    std::vector<double> parameters(std::size_t i)
    {
        return m_parameters[i];
    }

    double timeStep()
    {
        return m_timeStep;
    }

    double timeMax()
    {
        return m_timeMax;
    }

    std::string gaussType()
    {
        return m_gaussType;
    }

    std::string viewNames(std::size_t i)
    {
        return m_viewNames[i];
    }

    std::string timeMethod()
    {
        return m_timeMethod;
    }

    int timeRegistration()
    {
        return m_timeRegistration;
    }

    Parameters(std::string paramPath)
    {

        std::size_t i, j;

        int numParams;

        std::ifstream parametersFile;

        parametersFile.open(paramPath);

        // Getting the name of the msh file.
        getline(parametersFile, m_mshPath);

        parametersFile >> m_numUnknown;

        parametersFile.get();

        m_fluxName.resize(m_numUnknown);

        m_numFluxName.resize(m_numUnknown);

        m_parameters.resize(m_numUnknown);

        m_viewNames.resize(m_numUnknown);

        for(i = 0; i < m_numUnknown; ++i)
        {

            getline(parametersFile, m_fluxName[i]);

        }

        for(i = 0; i < m_numUnknown; ++i)
        {
            getline(parametersFile, m_numFluxName[i]);
        }

        parametersFile >> numParams;

        parametersFile.get();

        for(i = 0; i < m_numUnknown; ++i)
        {

            m_parameters[i].resize(numParams);

            for(j = 0; j < numParams; ++j)
            {

                parametersFile >> m_parameters[i][j];

                parametersFile.get();

            }

        }

        parametersFile >> m_timeStep;

        parametersFile.get();

        parametersFile >> m_timeMax;

        parametersFile.get();

        getline(parametersFile, m_gaussType);

        for(i = 0; i < m_numUnknown; ++i)
        {

            getline(parametersFile, m_viewNames[i]);
            
        }

        getline(parametersFile, m_timeMethod);

        parametersFile >> m_timeRegistration;

        parametersFile.get();

        parametersFile.close();

    }

};

#endif