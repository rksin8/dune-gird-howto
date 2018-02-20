// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// Dune includes
#include <config.h>           // file constructed by ./configure script /*@\label{gs:inc0}@*/
#include <array>
#include <memory>
#include <dune/common/parallel/mpihelper.hh> // include mpi helper class
#include <dune/grid/yaspgrid.hh> // load Yasp grid definition
#include <dune/grid/common/gridinfo.hh> // definition of gridinfo /*@\label{gs:inc1}@*/

int main(int argc, char **argv)
{
  // initialize MPI, finalize is done automatically on exit
  Dune::MPIHelper::instance(argc,argv);

  // start try/catch block to get error messages from dune
  try{
    // make a grid
    const int dim=3;                                     /*@\label{gs:dim}@*/
    typedef Dune::YaspGrid<dim> GridType;                /*@\label{gs:gridtype}@*/
    Dune::FieldVector<double,dim> length(1.0);           /*@\label{gs:par0}@*/
    std::array<int,dim> elements;
    std::fill(elements.begin(), elements.end(), 4);      /*@\label{gs:par1}@*/
    GridType grid(length,elements);                      /*@\label{gs:grid}@*/

    // print some information about the grid
    Dune::gridinfo(grid);
  }
  catch (std::exception & e) {
    std::cout << "STL ERROR: " << e.what() << std::endl;
    return 1;
  }
  catch (Dune::Exception & e) {
    std::cout << "DUNE ERROR: " << e.what() << std::endl;
    return 1;
  }
  catch (...) {
    std::cout << "Unknown ERROR" << std::endl;
    return 1;
  }

  // done
  return 0;
}
