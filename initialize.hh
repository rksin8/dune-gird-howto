/*
 * initialize.hh
 *
 *  Created on: 16-Feb-2018
 *      Author: root
 */

#ifndef INITIALIZE_HH_
#define INITIALIZE_HH_

#include<dune/common/version.hh>
#include<dune/geometry/referenceelements.hh>
#include "transportproblem2.hh"

// initialize the vector of unknowns with initial value
template<class G, class M, class V>
void initialize(const G& grid, const M& mapper, V& c)
{
	//first we extract the dimensions of the grid
	const int dim = G::dimension;
	const int dimworld = G::dimensionworld;


	typedef typename G::LeafGridView GridView;

	GridView gv = grid.leafGridView();

	//loop through leaf grid an evulate c0 at cell center

	for(const auto& cell:elements(gv)){
//		//get geometry type
//		auto gt = cell.type();
//
		auto geo = cell.geometry();

		//get the global coordinate of cell center
		auto global = geo.center();

//		//get cell center in reference element
//		auto local = cell.local(global);

#if DUNE_VERSION_NEWER(DUNE_GRID,2,4)
		//initialize the cell concentration
		c[mapper.index(cell)]=c0(global);
#else
		//initialize the cell concentration
		c[mapper.map(cell)]=c0(global);
#endif
	}



}





#endif /* INITIALIZE_HH_ */
