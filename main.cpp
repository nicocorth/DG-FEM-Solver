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
    std::vector<std::string> viewNames; // get the name of the final views.
    std::string timeMethod; // Contains the name of the temporal method.


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
                      gaussType,
                      viewNames,
                      timeMethod);


    gmsh::initialize(argc, argv);
    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::option::setNumber("Mesh.SaveAll", 1);
    gmsh::open(mshPath);


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

    Element mainElements(gaussType, elementType[0], dimTags[0].second);

    mainElements.frontierAndNeighbouring(gaussType);


    solver(numUnknown,
           mainElements,
           fluxName,
           numFluxName,
           parameters,
           timeStep,
           timeMax,
           viewNames,
           timeMethod);

    gmsh::finalize();

    return 0;
}