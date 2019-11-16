This is the binary Software Development Kit (SDK) for Gmsh 4.4.1:

  * Operating system: Linux64-sdk (Linux)
  * C++ compiler: /usr/local/bin/g++
  * C++ compiler ID: GNU
  * C++ compiler version: 5.5.0
  * C++ compiler flags:-fPIC -std=c++11 -O2 -g -DNDEBUG
  * Build options: 64Bit ALGLIB Ann Bamg Blas[custom] Blossom Cgns DIntegration Dlopen DomHex Fltk Gmm Hxt Hxt3D Jpeg[fltk] Kbipack Lapack[custom] LinuxJoystick MathEx Med Mesh Metis Mmg3d Mpeg Netgen ONELAB ONELABMetamodel OpenCASCADE OpenCASCADE-CAF OpenGL OptHom Parser Plugins Png[fltk] Post QuadTri Solver TetGen/BR Voro++ Zlib

Gmsh is distributed under the terms of the GNU General Public License: see
share/doc/gmsh/LICENSE.txt and share/doc/gmsh/CREDITS.txt. For additional Gmsh
resources, see http://gmsh.info.

SDK layout:

  * lib/*gmsh*.{so,dll,dylib}: shared Gmsh library
  * lib/gmsh.lib: import library (Windows only)
  * lib/gmsh.py: Python module
  * lib/gmsh.jl: Julia module
  * include/gmsh.h: C++ API header
  * include/gmshc.h: C API header
  * include/gmsh.h_cwrap: C++ wrapper of the C API (see the `Notes' below)
  * bin/gmsh: gmsh executable (linked with the shared Gmsh library)
  * share/doc/gmsh/demos/api : API examples in C++, C, Python and Julia

Notes:

  * The C API should work with most compilers.

  * The C++ API will only work if your compiler has an Application Binary
    Interface (ABI) that is compatible with the ABI of the compiler used to
    build this SDK (see above for the compiler ID and version).

  * If your C++ compiler does not have a compatible ABI and if there are no
    compatibility flags available, you can rename `gmsh.h_cwrap' as `gmsh.h':
    this implementation redefines the C++ API in terms of the C API. Using this
    header will lead to (slightly) reduced performance compared to using the
    native Gmsh C++ API from the original `gmsh.h' header, as it entails
    additional data copies between this C++ wrapper, the C API and the native
    C++ code.

    For example, the Windows SDK is currently compiled using the GNU Compiler
    Collection (GCC). To compile a C++ example with Microsoft Visual Studio 2017
    in the Visual Studio shell and run it, you would do:

    C:\gmsh-git-Windows64-sdk> ren include\gmsh.h_cwrap gmsh.h
    C:\gmsh-git-Windows64-sdk> cl /Iinclude share\doc\gmsh\demos\api\simple.cpp lib\gmsh.lib
    C:\gmsh-git-Windows64-sdk> cd lib
    C:\gmsh-git-Windows64-sdk\lib> ..\simple.exe

  * To make it as portable and as easy-to-use as possible the shared Gmsh
    library embeds most Gmsh dependencies statically linked directly inside the
    share library: libgfortran, OpenBLAS, FreeType, OpenCASCADE, OpenBLAS, HDF5,
    MED, CGNS, FLTK, ... Linking your app with different versions of some of
    these dependencies (e.g. linking your app with the MKL BLAS) will cause
    issues: in this case you should rebuild the Gmsh shared library from source
    and manage the dependencies consistently.

  * The shared Gmsh library also references shared system libraries. If some of
    these libraries are missing you will need to install them.
