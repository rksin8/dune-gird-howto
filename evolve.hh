/* evolve.hh
 *  Created on: 16-Feb-2018
 *      Author: ranjeet
 */
#ifndef EVOLVE_HH_
#define EVOLVE_HH_

#include <dune/common/version.hh>
#include <dune/geometry/referenceelements.hh>

template<class G, class M, class V>
void evolve(const G& grid, const M& mapper, V& c, double t,double& dt)
{
	//first we extract the dimensions of the grid
	const int dim = G::dimension;
	const int dimworld = G::dimensionworld;

	//type used for coordinates in the grid
	typedef typename G::ctype ct;

	//type of grid view on leaf part
	typedef typename G::LeafGridView GridView;

	//get grid view on leaf part
	GridView gv = grid.leafGridView();

	// allocate a temporary vector for update
	V update(c.size());

	for (typename V::size_type i = 0; i < c.size(); ++i) {
		update[i]=0;
	}

	//initialize dt very large
	dt=1E100;

	//compute update vector and optimum in one grid traversal
	for(const auto& cell:elements(gv)){
//		//get geometry type
		auto gt = cell.type();
		//cell geometry
      	auto geo = cell.geometry();

		//get the global coordinate of cell center
		auto global = geo.center();

		//get cell center in reference element
		auto local = geo.local(global);

		//cell volume, assume linear map here
		double volume = geo.volume();


#if DUNE_VERSION_NEWER(DUNE_GRID,2,4)
		//cell index
		int indexi = mapper.index(cell);
#else
		//cell index
		int indexi = mapper.map(cell);
#endif

		//variable to compute sum of positive factors
		double sumfactor=0.0;

		//run through all intersections with neighbors and boundary
		for(const auto& is: intersections(gv,cell))
		{
			//get the geometry type of face
            		auto is_geo = is.geometry();
			// get geometry type of face
			auto gtf = is_geo.type();

			//center in face's reference element
			auto faceglobal = is_geo.center();
            // center in the local coordinate
			auto facelocal = is_geo.local(faceglobal);

			//get normal vector scaled with volume
			auto integrationOuterNormal = is.integrationOuterNormal(facelocal);

			integrationOuterNormal *= is_geo.volume();

			//evaluate velocity at face center
			auto velocity = u(faceglobal,t);

			//compute factor occuring in flux formula
			double factor = velocity*integrationOuterNormal/volume;

			//for time step calculation
			if(factor>=0) sumfactor+=sumfactor;

			//handle interior face
			if(is.neighbor()){
				//Code for Discontinuous Galerkin or Finite Volume
				//access neighbor
				const auto& outside = is.outside();

#if DUNE_VERSION_NEWER(DUNE_GRID,2,4)
				int indexj = mapper.index(outside);
#else
				int indexj = mapper.map(outside);
#endif

			// compute flux from one side only
				// this should become easier with the new
				// IntersectionIterator functionality

				if(cell.level()>outside.level()||
						(cell.level()==outside.level()
						&& indexi < indexj))
				{
					//compute factor in neighbor
					Dune::GeometryType nbgt = outside.type();
					const Dune::FieldVector<ct, dim>&
					nblocal = Dune::ReferenceElements<ct, dim>::general(nbgt).
					position(0,0);

					double nbvolume = outside.geometry().integrationElement(nblocal)
					 *Dune::ReferenceElements<ct, dim>::general(nbgt).volume();

					double nbfactor = velocity*integrationOuterNormal/nbvolume;

					if(factor<0)// inflow
					{
						update[indexi]-=c[indexj]*factor;
						update[indexj]+=c[indexj]*nbfactor;
					}else //outflow
					{
						update[indexi]-=c[indexi]*factor;
						update[indexj]+=c[indexi]*nbfactor;
					}
				}
			}

			//handle boundary face
			if(is.boundary())
			{//handle potential Neumann boundary
				if(factor<0)//inflow, apply boundary condition
					update[indexi]-= b(faceglobal,t)*factor;
				else //outflow
					update[indexi]-=c[indexi]*factor;
			}

		} // end all intersections

		//compute dt restriction
		dt = std::min(dt, 1.0/sumfactor);
	} // end grid traversal

	//scale dt with safety factor
	dt *=0.99;

	//update the concentration vector
	for (int i = 0; i < c.size(); ++i) {
		c[i]+=dt*update[i];
	}
	return ;
}
#endif /* EVOLVE_HH_ */
