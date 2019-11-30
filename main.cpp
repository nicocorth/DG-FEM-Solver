#include <cstdio>
#include <iostream>
#include <gmsh.h>
#include <cmath>
#include <vector>
#include <fstream>

#include "gmsh.h"
#include "mesh.hpp"
#include "solver.hpp"
#include "parameters.hpp"

// Main of the program.

int main(int argc, char **argv)
{

    Parameters simulInfos(argv[1]);

    std::vector<int> elementType;
    std::vector<std::vector<std::size_t>> bin, bin1;
    gmsh::vectorpair dimTags;

    gmsh::initialize(argc, argv);
    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::option::setNumber("Mesh.SaveAll", 1);
    gmsh::open(simulInfos.mshPath());


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

    Element mainElements(simulInfos.gaussType(), elementType[0], dimTags[0].second);

    mainElements.frontierAndNeighbouring(simulInfos.gaussType());

    solver(simulInfos, mainElements);

    gmsh::finalize();

    return 0;
}