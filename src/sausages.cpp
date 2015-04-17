#include <iostream>
#include <vector>
#include <numeric>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <limits>
#include "io.h"
#include "main.h"
#include "maths.h"
#include "params.h"
#include "point.h"
#include "sausages.h"

using namespace std;

Sausage::Sausage(int ID){
	sausageID=ID;
	have_rotation_matrix = false;
	have_pobf = false;
}


/** Ordinary Flood-fill algorithm
 * 1) Find the first sausage-worthy point not yet sorted into a sausage
 * 2) Add it to a numbered sausage, and add all its unsorted neighbours to a stack
 * 3) While there are points in the stack, take the last one, add it to
 *     the sausage and its unsorted neighbours to the stack
 * 4) Eventually you'll have a continuous sausage, start with the next one.
 */
void flood_fill_separate(vector<Point*> allPoints, vector<Sausage> &allSausages){
	debug() << "begin flood_fill_separate" << endl << flush;

	int newSausageID=0;
	vector<size_t> stack; // Empty FILO stack, filled with indices of points to be coloured
	for (size_t firstPoint=0; firstPoint<allPoints.size(); firstPoint++){
		// Pick the first sausage not yet coloured
		if (allPoints[firstPoint]->isInASausage && allPoints[firstPoint]->sausageID==-1){
			verbose() << "Starting flood-fill of sausage #" << newSausageID << endl << flush;
			Sausage newSausage(newSausageID);
			stack.push_back(firstPoint);
			while (stack.size()>0){
				Point *curr = allPoints[stack.back()];
				stack.pop_back();

				curr->sausageID=newSausageID;
				newSausage.points.push_back(curr);
				curr->sausagePointsIndex=newSausage.points.size()-1;
				for (int iNeigh=0;iNeigh<6;iNeigh++){
					Point *neigh=curr->neighbours[iNeigh];
					if (neigh  && neigh->isInASausage && neigh->sausageID==-1 && // neigh is sausage-worthy but not yet sorted
					    !elementOf(stack,neigh->allPointsIndex)){ // Don't add to the stack if it's already there
						stack.push_back(neigh->allPointsIndex) ;
					}
				}
			}
			allSausages.push_back(newSausage);
			debug()<<allSausages[newSausageID].points.size()<<" points in the sausages."<<endl<<flush;
			newSausageID++;
		}
	}
}

/** Classifies a given sausage
 * Four points on the sausage are found, two on either side of each colloid.
 * From each point we flood-fill towards the centre, and see which of the other points we reach
 *
 * Class 5: If we reach only the other colloid's point of the same side (i.e. above-1 -> above-2) then we
 *  have a simple ring with no twists, such as an omega configuration.
 *  Also classified as class-5 will be 2nd-loop systems where the 2nd loop is touching the 1st on only one side.
 *
 * Class 1/2: If we end up only on the other side of the opposite colloid (above-1 -> below-2)
 *  then we have a figure-of-eight system. 1 is left-handed, 2 is right-handed
 *
 * Class 0: If we reach more than one other point we have junctions on both sides of the sausage, and are
 *  insufficiently resolved.
 *
 */
int Sausage::flood_fill_classify(const std::vector<Vector3d> colloidPos){
	/* We need to find the vector parallel to the sausage's PoBF which is perpendicular to the line joining the two colloids.
	 * Then we follow this vector out from a colloid to find suitable points on the sausage.
	 * If we arbitrarity set the vector to be length 1 we can determine it from the colloid positions and the PoBF.
	 */

	debug() << "begin flood_fill_classify" << endl << flush;

	// AB is vector joining colloids
	Vector3d AB = colloidPos[0] - colloidPos[1];
	double norm_AB = sqrt(dot(AB,AB));
	debug()<<"Vector AB is "<<AB<<" of length "<<norm_AB<<endl;
	if (norm_AB < params::epsilon){
		cerr << "Vector between colloids is very small, check inputs! Aborting." << endl;
		exit(EXIT_FAILURE);
	}

	// Plane of best fit is of form alpha*x+beta*y=z
	double alpha=-plane_of_best_fit[0];
	double beta=-plane_of_best_fit[1];

	// From being parallel to PoBF, perpendicular to AB, we can determine:
	Vector3d v (0,0,0);
	if (abs(AB.y+beta*AB.z) > params::epsilon){
		v.x=1; // No need to normalise
		v.y=(-alpha-AB.x)/(AB.y+beta*AB.z);
	}else{
		v.x=(-AB.y-beta*AB.z)/(AB.x+alpha); // = 0, more or less
		v.y=1.0/sqrt(1+beta*beta); // This pretty much normalises for free
	}
	v.z=alpha*v.x+beta*v.y;

	double norm_v = sqrt(dot(v,v));
	debug()<<"Vector v (perp. AB & || PoBF & unit) is "<<v<<" of length "<<norm_v<<endl;

	// For each point in the sausage, is it between two parallel planes extended from beside a colloid, perpendicular to the line between the colloids?
	// If so, add it to a region dependant on whether it is v-wards or anti-v-wards
	vector<Point*> above[2], below[2];
	for (vector<Point*>::iterator iter=points.begin(); iter!=points.end(); iter++){
		Point* it = *iter;
		for (int iColl=0; iColl<2; iColl++){
			Vector3d p (it->x, it->y, it->z); // xyz of this point

			Vector3d u = p - colloidPos[iColl]; // Vector from point to colloid

			// Eq. of plane perpendicular to AB going through point p is (x,y,z).AB - p.AB
			// Distance from point q to this plane P is |Px*qx + Py*qy + Pz*qz + P0|/norm(Px,Py,Pz)
			Vector3d first_plane  = colloidPos[iColl] + 0.5*params::flood_fill_classify_slice_size*params::pixel_size*(AB/norm_AB);
			Vector3d second_plane = colloidPos[iColl] - 0.5*params::flood_fill_classify_slice_size*params::pixel_size*(AB/norm_AB);
			double first_distance  = abs( dot(AB, p) - dot(first_plane,  AB) ) / norm_AB;
			double second_distance = abs( dot(AB, p) - dot(second_plane, AB) ) / norm_AB;

			// Is it between the two planes? If so, it is <= flood_fill_classify_slice_size from each plane
			if (first_distance < params::flood_fill_classify_slice_size*params::pixel_size &&
			   second_distance < params::flood_fill_classify_slice_size*params::pixel_size ){
				double projection = dot(u,v);
				if (projection>0){
					above[iColl].push_back(it);
					debug() << p << " added to region above colloid " << iColl << endl;
				} else {
					below[iColl].push_back(it);
					debug() << p << " added to region below colloid " << iColl << endl;
				}
			}
		}
	}

	verbose()<<above[0].size()<<" pixels in region above first colloid"<<endl;
	verbose()<<above[1].size()<<" pixels in region above second colloid"<<endl;
	verbose()<<below[0].size()<<" pixels in region below first colloid"<<endl;
	verbose()<<below[1].size()<<" pixels in region below second colloid"<<endl;

	if (below[0].size() * below[1].size() * above[0].size() * above[1].size() == 0){
		cerr << "One of the regions is empty, this is bad, exiting." << endl << flush;
		exit(EXIT_FAILURE);
	}

	/** Adjacency matrix of region
	 * 4-by-4 adjacency matrix. Each region is adjacent to itself.
	 * Regions are numbered:  0) below-0  1) above-0  2) below-1  3) above-1
	 * i.e. region is number 2*iColl+aboveOrBelow
	 */
	int adjacency[4][4] = {}; // initalise to all zero

	/** Set of all points reached in the flood-fills from regions above/below colloid 0.
	 * We use these to check whether a twist is left-handed (above-0 -> below-1 has higher z than above-1 -> below-0) or right-handed
	 */
	vector<Point*> fromAbove0;
	vector<Point*> fromBelow0;

	// For each region
	for (int iColl=0;iColl<2;iColl++){ for (int aboveOrBelow=0;aboveOrBelow<2;aboveOrBelow++){ // aboveOrBelow is 0 for below, 1 for above
		debug()<<"in region: colloid "<<iColl<<" aboveOrBelow "<<aboveOrBelow<<endl<<flush;
		vector<Point*> region = aboveOrBelow ? above[iColl] : below[iColl];
		int otherColl=(iColl+1)%2;

		// find all neighbours of the region which are inside this sausage and closer to the /other/ colloid
		vector<Point*> neighbours;
		for (vector<Point*>::iterator iter=region.begin(); iter!=region.end(); iter++){
			Point* it=*iter;
			double distanceToOtherColl = pow(it->x - colloidPos[otherColl].x,2) +
							pow(it->y - colloidPos[otherColl].y,2) +
							pow(it->z - colloidPos[otherColl].z,2) ;
			for (int iNeigh=0;iNeigh<6;iNeigh++){
				Point* neigh = it->neighbours[iNeigh];
				if (neigh && neigh->isInASausage && !elementOf(region,neigh)){
					double neighDistanceToOtherColl = pow(neigh->x - colloidPos[otherColl].x,2) +
									pow(neigh->y - colloidPos[otherColl].y,2) +
									pow(neigh->z - colloidPos[otherColl].z,2) ;
					if (neighDistanceToOtherColl < distanceToOtherColl) neighbours.push_back(it->neighbours[iNeigh]);
				}
			}
		}
		debug()<<"end of region"<<endl<<flush;

		// Get rid of duplicate neighbours
		debug()<<region.size()<<" points in this region, "<<neighbours.size()<<" neighbours inside the sausage."<<endl;
		sort(neighbours.begin(),neighbours.end());
		debug()<<"sorted neighbours"<<endl;
		vector<Point*>::iterator it = unique(neighbours.begin(),neighbours.end());
		neighbours.resize(distance(neighbours.begin(),it));
		debug()<<"Removed duplicates, now "<<neighbours.size()<<" neighbours."<<endl;
		debug()<<"Neighbours at:"<<endl;
		for (vector<Point*>::iterator it=neighbours.begin(); it!=neighbours.end(); it++){
			debug()<<(*it)->x<<","<<(*it)->y<<","<<(*it)->z<<endl;
		}

		// flood-fill from neighbours. If a point is in a region, make a note but don't add it to the flood-filled area.
		vector<Point*> stack=neighbours; // FILO stack, to be filled with contiguous points
		vector<Point*> visited; // Points we've already visited, and so should ignore
		while (stack.size()>0){
			//debug()<<"Another pixel on the stack"<<endl;
			Point* curr = stack.back();
			stack.pop_back();

			for (int iNeigh=0;iNeigh<6;iNeigh++){
				if (!curr->neighbours[iNeigh]) continue;
				Point* thisNeigh=curr->neighbours[iNeigh];
				if (elementOf(visited,thisNeigh)) continue;
				visited.push_back(thisNeigh);
				//debug()<<"neigh: "<<iNeigh<<" address: "<<curr.neighbours[iNeigh]<<" sID: "<<curr.neighbours[iNeigh]->sausageID<<endl;
				if (thisNeigh->isInASausage) {
					bool in_a_region=false;
					// For all of regions...
					for (int jColl=0;jColl<2;jColl++){ for (int inner_aboveOrBelow=0;inner_aboveOrBelow<2;inner_aboveOrBelow++){
						vector<Point*> inner_region = inner_aboveOrBelow ? above[jColl] : below[jColl];
						// If the neighbour is in the region, make a note but don't add to stack
						if (elementOf(inner_region,thisNeigh)){
							in_a_region=true;
							adjacency[2*iColl+aboveOrBelow][2*jColl+inner_aboveOrBelow]=1;
							//debug()<<"I'm in a region"<<endl;
							//debug()<<"From "<<2*iColl+aboveOrBelow<<" to "<<2*jColl+inner_aboveOrBelow<<endl;
						}
					}}
					if (!in_a_region && !elementOf(stack,thisNeigh)) stack.push_back(thisNeigh);
				}
			}
		}
		verbose()<<"Finished flood-fill from region "<<(aboveOrBelow?"above":"below")<<" colloid "<<iColl<<endl;
		//write_points("debug_aboverOrBelow"+to_string(aboveOrBelow)+"_colloid"+to_string(iColl)+".dat", visited);

		// Save paths for identifying twist handedness
		if (aboveOrBelow==0 && iColl==0){
			fromBelow0=visited;
		} else if (aboveOrBelow==1 && iColl==0){
			fromAbove0=visited;
		}

	}} // end of 'for each region'

	for (int i=0;i<4;i++){ for (int j=0;j<4;j++){
		debug()<<"Adjacency["<<i<<"]["<<j<<"] = "<<adjacency[i][j]<<endl;
	}}

	for (int i=0;i<4;i++){ for (int j=i;j<4;j++){
		if (adjacency[i][j]!=adjacency[j][i]){
			error()<<"Adjacency matrix is not symmetric!"<<endl;
			exit(EXIT_FAILURE);
		}
	}}
	debug()<<"Adjacency matrix is correctly symmetric."<<endl;
	for (int i=0;i<4;i++){
		int sum=0;
		for (int j=0;j<4;j++){
			sum+=adjacency[i][j];
		}
		if (sum<2){ // Including the region itself (i.e. adjacency[i][i]=1)
			error()<<"Region "<<(i%2?"above":"below")<<" colloid "<<(i/2)<<" doesn't appear to be linked to any others."<<endl;
			exit(EXIT_FAILURE);
		} else if (sum>2){
			error()<<"Region "<<(i%2?"above":"below")<<" colloid "<<(i/2)<<" appears to be connected to more than one other region, inspect manually."<<endl;
			exit(EXIT_FAILURE);
		}
	}
	if (adjacency[0][2] && adjacency[1][3]){
		info()<<"System appears to be untwisted (Class 5)"<<endl;
		return 5;
	} else if (adjacency[0][3] && adjacency[1][2]){
		info()<<"System appears to be twisted (Class 1/2)"<<endl;
		// Work out whether it's a RHS twist or a LHS one
		int handedness = find_twist_handedness(fromBelow0, fromAbove0);
		if (handedness==1){
			info()<<"Looks left-handed (Class 1)"<<endl;
			return 1;
		} else if (handedness==2) {
			info()<<"Looks right-handed (Class 2)"<<endl;
			return 2;
		} else {
			error()<<"Something went wrong finding handedness of the twist"<<endl;
			exit(EXIT_FAILURE);
		}
	} else {
		brief({1}) << "Structure undefined." << endl;
		brief({1}) << "0" << endl;
		error()<<"Unexpected system linkage, this should never happen."<<endl;
		exit(EXIT_FAILURE);
	}
}

/** Takes two vectors of Points, corresponding to the two sets of all points reached in the flood-fills from regions above/below colloid 0.
 * We use these to check whether a twist is left-handed (return 1) (above-0 -> below-1 has higher z than above-1 -> below-0) or right-handed (return 2)
 */
int Sausage::find_twist_handedness(std::vector<Point*> fromBelowPoints, std::vector<Point*> fromAbovePoints){

	// Make vector<Vector3d> of the points in each arm.
	/*
	 * This would be nice (implicit casting) but I CBA to go back over flood_fill_classify and replace <Point*> with <Point&> or whatev2r
	std::vector<Vector3d> fromBelow(fromBelowPoints.begin(), fromBelowPoints.end());
	std::vector<Vector3d> fromAbove(fromAbovePoints.begin(), fromAbovePoints.end());
	*/
	std::vector<Vector3d> fromBelow, fromAbove;
	for (vector<Point*>::iterator iter=fromBelowPoints.begin(); iter!=fromBelowPoints.end(); iter++){
		fromBelow.push_back((Vector3d)(**iter));
	}
	for (vector<Point*>::iterator iter=fromAbovePoints.begin(); iter!=fromAbovePoints.end(); iter++){
		fromAbove.push_back((Vector3d)(**iter));
	}

	rotate_to_xy_plane(fromBelow);
	rotate_to_xy_plane(fromAbove);


	// Points are within a cross-over region if their xy projection is 'close' to one of a point from the other arm
	// 'Close' is defined as 2*pixel_size
	double above_minz = numeric_limits<double>::infinity();
	double above_maxz = -numeric_limits<double>::infinity();
	double below_minz = numeric_limits<double>::infinity();
	double below_maxz = -numeric_limits<double>::infinity();
	bool overlap=false;
	for (vector<Vector3d>::iterator above=fromAbove.begin(); above!=fromAbove.end(); above++){
		for (vector<Vector3d>::iterator below=fromBelow.begin(); below!=fromBelow.end(); below++){
			if ( pow(above->x-below->x, 2) + pow(above->y-below->y, 2) < pow(2*params::pixel_size,2) ){
				//debug() << *below << " and " << *above << " are in the DANGER ZONE! <cough> crossing zone." << endl;
				overlap=true;
				above_minz = (above->z < above_minz) ? above->z : above_minz;
				above_maxz = (above->z > above_maxz) ? above->z : above_maxz;
				below_minz = (below->z < below_minz) ? below->z : below_minz;
				below_maxz = (below->z > below_maxz) ? below->z : below_maxz;
			}
		}
	}
	if (above_minz>above_maxz || below_minz>below_maxz || !overlap){
		error() << "Something's wrong with the crossing-zone, aborting" << endl;
		exit(EXIT_FAILURE);
	}
	if (above_minz > below_maxz){
		info()<< "Looks like a left-handed figure-of-eight" << endl;
		return 1;
	} else if (above_maxz < below_minz){
		info()<< "Looks like a right-handed figure-of-eight" << endl;
		return 2;
	} else {
		error() << "Something's wrong with the crossing-zone" << endl;
		return 0;
	}
}

void Sausage::find_endpoints(void){

	debug() << "begin find_endpoints of sausage " << sausageID << endl << flush;

	// There is no 'infinity' for integers, like there is for doubles.
	// However, there is a max(), which works out to be something like 2,147,000,000
	// This is (hopefully) much larger than any actual distance across the network
	// However, max+max gives an overflow, so we need something where inf+inf>inf
	int myInfinity = numeric_limits<int>::max()/2-1;
	debug()<<"myInf: "<<myInfinity<<", myInf+myInf: "<<myInfinity+myInfinity<<endl;

	// dist is |V| x |V| array of min distances, initialised to infinity
	int** dist = new int* [points.size()];
	for (size_t i=0; i<points.size(); i++){
		dist[i] = new int [points.size()];
		for (size_t j=0; j<points.size(); j++){
			dist[i][j]=myInfinity;
		}
	}
	if (points.size()>15) debug()<<"Arbitrary initial distance: "<<dist[10][15]<<endl;

	// For each vertex v, dist[v][v]=0
	// For each edge (u,v), dist[u][v]=w(u,v)=1 in this case
	for (size_t v=0; v<points.size(); v++){
		dist[v][v]=0;
		for (int iNeigh=0;iNeigh<6;iNeigh++){
			if (NULL==points[v]->neighbours[iNeigh] || points[v]->neighbours[iNeigh]->sausageID != sausageID) continue;
			size_t u=points[v]->neighbours[iNeigh]->sausagePointsIndex;
			dist[u][v]=1;
		}
	}

	// shortestPath(i,j,k+1) = min ( shortestPath(i,j,k), shortestPath(i,k+1,k)+shortestPath(k+1,j,k) )
	for (size_t k=0; k<points.size(); k++){
		if (k%100==0) debug()<<"k="<<k<<"/"<<points.size()<<endl<<flush;
		for (size_t i=0; i<points.size(); i++){
			for (size_t j=0; j<points.size(); j++){
				if (dist[i][j] > dist[i][k] + dist[k][j]){
					dist[i][j] = dist[i][k] + dist[k][j];
				}
			}
		}
	}

	// Find maximum minimum-distance
	int max_min_dist=0;
	endpoints[0] = endpoints[1] = 0; // In case there's only a single point
	for (size_t i=0; i<points.size(); i++){
		for (size_t j=i; j<points.size(); j++){
			if (dist[i][j]>max_min_dist){
				max_min_dist=dist[i][j];
				endpoints[0]=i;
				endpoints[1]=j;
				//debug()<<"dist: "<<max_min_dist<<", i: "<<i<<", j:"<<j<<endl;
			}
		}
	}

	debug() << "Endpoint 0 is #" << endpoints[0] << " at: " << *points[endpoints[0]] << endl;
	debug() << "Endpoint 1 is #" << endpoints[1] << " at: " << *points[endpoints[1]] << endl;
	debug() << "Minimum distance along the sausage between them is: " << max_min_dist << endl;

	for (size_t i=0; i<points.size(); i++){
		delete [] dist[i];
	}
	delete [] dist;
}

void join_endpoints(Model &model){
	// Work out distance between endpoints
	double dist[2*model.allSausages.size()][2*model.allSausages.size()];
	for (size_t i=0; i<model.allSausages.size(); i++){
		Sausage si = model.allSausages[i];
		for (size_t iEndpoint=0; iEndpoint<2; iEndpoint++){
			for (size_t j=i; j<model.allSausages.size(); j++){ // Only need half, by symmetry
				Sausage sj = model.allSausages[j];
				for (size_t jEndpoint=0; jEndpoint<2; jEndpoint++){
					if (i==j && iEndpoint>=jEndpoint) continue; // More symmetry
					float pi[3],pj[3],p2p[3];
					pi[0]=si.points[si.endpoints[iEndpoint]]->x;
					pi[1]=si.points[si.endpoints[iEndpoint]]->y;
					pi[2]=si.points[si.endpoints[iEndpoint]]->z;
					pj[0]=sj.points[sj.endpoints[jEndpoint]]->x;
					pj[1]=sj.points[sj.endpoints[jEndpoint]]->y;
					pj[2]=sj.points[sj.endpoints[jEndpoint]]->z;
					p2p[0]=pi[0]-pj[0];
					p2p[1]=pi[1]-pj[1];
					p2p[2]=pi[2]-pj[2];

					dist[2*i+iEndpoint][2*j+jEndpoint]=sqrt(inner_product(p2p,p2p+3,p2p,0.0));
					dist[2*j+jEndpoint][2*i+iEndpoint]=dist[2*i+iEndpoint][2*j+jEndpoint];
					debug()<<"Distance "<<i<<":"<<iEndpoint<<","<<j<<":"<<jEndpoint<<","<<dist[2*i+iEndpoint][2*j+jEndpoint]<<endl;
				}
			}
		}
	}

	// Join endpoints that are 'close' to each other
	size_t other_endpoint;
	bool nearby_endpoint;
	for (size_t i=0; i<model.allSausages.size(); i++){
		for (size_t iEndpoint=0; iEndpoint<2; iEndpoint++){
			nearby_endpoint=false;
			int num_nearby_endpoints=0;
			for (size_t j=i; j<model.allSausages.size(); j++){ // Only need half, by symmetry
				for (size_t jEndpoint=0; jEndpoint<2; jEndpoint++){
					if (i==j && iEndpoint>=jEndpoint) continue; // More symmetry
					// If we've already closed this sausage into a loop, the endpoints will be the same point
					if (dist[2*i+iEndpoint][2*j+jEndpoint]<params::epsilon) continue;
					// Are there any points close enough to join to?
					if (dist[2*i+iEndpoint][2*j+jEndpoint]<params::max_sausage_gap_size*params::pixel_size){
						debug()<<"Picked out "<<j<<":"<<jEndpoint<<" and "<<i<<":"<<iEndpoint<<endl;
						nearby_endpoint=true;
						other_endpoint=2*j+jEndpoint;
					}
					// Are there other points which are also nearby, which make this hard to decide on?
					if (dist[2*i+iEndpoint][2*j+jEndpoint]<params::min_sausage_gap_next_neigh_distance*params::pixel_size){
						num_nearby_endpoints++;
					}

				}
			}
			if (num_nearby_endpoints>1){
				error()<<"Too many endpoints near endpoint "<<iEndpoint<<" of sausage "<<i<<endl;
				exit(EXIT_FAILURE);
			}
			if (nearby_endpoint){
				// Join endpoints i:iEndpoint and j:jEndpoint
				// To do this:
				//  - Add sausagej.points to sausagei, changing their sausageID and sausagePointsIndex (unless sausagej==sausagei)
				//  - Recalculate sausagei's endpoints
				//  - Remove sausagej from relevantSausages (unless sausagej==sausagei)

				size_t j=other_endpoint/2;
				size_t jEndpoint=other_endpoint%2;
				Sausage &sausagei=model.allSausages[i];
				Sausage &sausagej=model.allSausages[j];

				//  - Add sausagej.points to sausagei, changing their sausageID and sausagePointsIndex (unless sausagej==sausagei)
				int original_sausagei_size = sausagei.points.size();
				if (sausagej.sausageID!=sausagei.sausageID){
					for (vector<Point*>::iterator iter=sausagej.points.begin(); iter!=sausagej.points.end(); ++iter){
						(*iter)->sausageID=sausagei.sausageID;
						(*iter)->sausagePointsIndex+=sausagei.points.size();
					}
					debug()<<"sausagei size: "<<sausagei.points.size()<<endl;
					debug()<<"sausagej size: "<<sausagej.points.size()<<endl;
					sausagei.points.insert(sausagei.points.end(),sausagej.points.begin(),sausagej.points.end());
					debug()<<"new sausagei size: "<<sausagei.points.size()<<endl;
				}

				// Artificially add points between nearby endpoints
				Point *curr = sausagei.points[sausagei.endpoints[iEndpoint]];
				Point *target = sausagej.points[sausagej.endpoints[jEndpoint]];

				debug() << "Joining, heading from point " << *curr << " to point " << *target << endl;

				// We can't just follow links anymore, since we don't save non-sausage points when using diot.
				// If a Point is missing, its neighbours will have Null links to it. We have to sort out the links here.
				// We don't want to have to check all of allPoints at each step, so we'll pre-compute the Points that might
				// be neighbours of a new Point. In this case it's the small area around the two endpoints.
				vector<Point*> possNeighs;
				double border = 2*params::pixel_size;
				debug() << "Border size: " << border << endl;
				for (Point* p : model.allPoints){
					if (p->x > min(curr->x, target->x)-border &&
					    p->x < max(curr->x, target->x)+border &&
					    p->y > min(curr->y, target->y)-border &&
					    p->y < max(curr->y, target->y)+border &&
					    p->z > min(curr->z, target->z)-border &&
					    p->z < max(curr->z, target->z)+border )
						possNeighs.push_back(p); // WARNING This is probably safe, allPoints won't get resized here
				}
				debug() << "possNeighs has size: " << possNeighs.size() << endl;
				for (Point *p : possNeighs){
					debug() << *p << endl;
				}



				// At each step we move along the dimension which takes us fastest to target.
				while (curr->x!=target->x || curr->y!=target->y || curr->z!=target->z){
					debug() << "Joining, now at point " << *curr << endl;
					if (curr->sausageID != sausagei.sausageID && curr->sausageID != sausagej.sausageID && curr->isInASausage){
						error()<<"Looks like the gap spans a third sausage, this is not good, aborting."<<endl;
						exit(EXIT_FAILURE);
					}
					// Add it to the sausage if needed
					if (curr->sausageID!=sausagei.sausageID){
						debug() << "Adding point to sausage" << endl;
						curr->isInASausage=true;
						curr->sausageID=sausagei.sausageID;
						sausagei.points.push_back(curr);
						curr->sausagePointsIndex=sausagei.points.size()-1;
					}
					double dx = abs(curr->x - target->x);
					double dy = abs(curr->y - target->y);
					double dz = abs(curr->z - target->z);

					enum direction {L,R,U,D,F,B};

					Point* next;
					direction dir;
					if (dx>dy && dx>dz){
						next = (curr->x > target->x) ? curr->left : curr->right;
						dir  = (curr->x > target->x) ? L : R;
					} else if (dy>dz) {
						next = (curr->y > target->y) ? curr->down : curr->up;
						dir  = (curr->y > target->y) ? D : U;
					} else {
						next = (curr->z > target->z) ? curr->back : curr->forward;
						dir  = (curr->z > target->z) ? B : F;
					}
					debug() << "next is: " << next << endl;

					// Do we need to create a new point?
					if (next==NULL){
						Point &p = *(new Point ());
						p.x=curr->x; p.y=curr->y; p.z=curr->z;
						// WARNING cl, cp, cs are unknown, we threw it away when reading
						p.cl = p.cp = p.cs = -1;

						switch (dir){
							case L: p.x -= params::pixel_size; break;
							case R: p.x += params::pixel_size; break;
							case D: p.y -= params::pixel_size; break;
							case U: p.y += params::pixel_size; break;
							case B: p.z -= params::pixel_size; break;
							case F: p.z += params::pixel_size; break;
						}

						// Look for all 6 neighbours by looking in possNeighs for Points in the right place
						for (direction d : {L,R,U,D,F,B} ){
							// Where should we look?
							Vector3d neighLoc = p;
							switch (d){
								case L: neighLoc.x -= params::pixel_size; break;
								case R: neighLoc.x += params::pixel_size; break;
								case D: neighLoc.y -= params::pixel_size; break;
								case U: neighLoc.y += params::pixel_size; break;
								case B: neighLoc.z -= params::pixel_size; break;
								case F: neighLoc.z += params::pixel_size; break;
							}
							// Which point is there?
							bool foundNeigh=false;
							for (Point* candidate : possNeighs){
								if (mag((Vector3d) *candidate - neighLoc) < (params::epsilon * params::epsilon)){
									foundNeigh=true;
									switch (d){
										case L:
											p.left = candidate;
											p.left->right = &p;
											p.left->neighbours[1] = &p;
											break;
										case R:
											p.right = candidate;
											p.right->left = &p;
											p.right->neighbours[0] = &p;
											break;
										case D:
											p.down = candidate;
											p.down->up = &p;
											p.down->neighbours[2] = &p;
											break;
										case U:
											p.up = candidate;
											p.up->down = &p;
											p.up->neighbours[3] = &p;
											break;
										case B:
											p.back = candidate;
											p.back->forward = &p;
											p.back->neighbours[4] = &p;
											break;
										case F:
											p.forward = candidate;
											p.forward->back = &p;
											p.forward->neighbours[5] = &p;
											break;
									}
									break; // We've found the neighbour
								}
							}
							if (!foundNeigh){
								debug() << "Couldn't find suitable neighbour at " << neighLoc << endl;
							}
						}

						// Set neighbours array
						p.neighbours[0]=p.left;
						p.neighbours[1]=p.right;
						p.neighbours[2]=p.up;
						p.neighbours[3]=p.down;
						p.neighbours[4]=p.forward;
						p.neighbours[5]=p.back;

						// Add it to allPoints
						model.allPoints.push_back(&p);
						p.allPointsIndex=model.allPoints.size()-1;

						// Also add it to possNeighs
						possNeighs.push_back(&p);


						debug() << "Made a new point at " << p << endl;

						// And use it!
						next=&p;
					};

					// Move to next point and continue
					curr=next;
				}
				debug() << "Sausagei size: " << sausagei.points.size() << endl;

				debug() << "Before endpoint-swap..." << endl;
				debug() << "i[0]: " << *(sausagei.points[sausagei.endpoints[0]]) << endl;
				debug() << "i[1]: " << *(sausagei.points[sausagei.endpoints[1]]) << endl;
				debug() << "j[0]: " << *(sausagej.points[sausagej.endpoints[0]]) << endl;
				debug() << "j[1]: " << *(sausagej.points[sausagej.endpoints[1]]) << endl;
				debug() << "join-point of i: " << sausagei.endpoints[iEndpoint] << endl;
				debug() << "join-point of i: " << *(sausagei.points[sausagei.endpoints[iEndpoint]]) << endl;
				debug() << "non-join-point of j: " << sausagej.endpoints[(jEndpoint+1)%2] << endl;
				debug() << "non-join-point of j: " << *(sausagej.points[sausagej.endpoints[(jEndpoint+1)%2]]) << endl;;

				//  Sort out new endpoints & remove sausagej if it's been merged into sausagei,
				if (sausagei.sausageID == sausagej.sausageID){
					// We've closed the sausage into a loop, make distance between endpoints = 0 to stop more joining
					sausagei.endpoints[1] = sausagei.endpoints[0];
				} else {
					// If a->b is longest path in saus i, and c->d is longest in saus j, then a->d will be longest if we link b&c
					// Need to make joining-endpoint of sausagei now equal to non-joining endpoint of sausagej
					sausagei.endpoints[iEndpoint] = sausagej.endpoints[(jEndpoint+1)%2] + original_sausagei_size;
					debug() << "Removing sausage " << j << " from allSausages" << endl;
					model.allSausages.erase(model.allSausages.begin()+j); // Awkward, but erase() requires an iterator
				}
				debug() << "After endpoint-swap..." << endl;
				debug() << "i[0]: " << *(sausagei.points[sausagei.endpoints[0]]) << endl;
				debug() << "i[1]: " << *(sausagei.points[sausagei.endpoints[1]]) << endl;
				debug() << "j[0]: " << *(sausagej.points[sausagej.endpoints[0]]) << endl;
				debug() << "j[1]: " << *(sausagej.points[sausagej.endpoints[1]]) << endl;
				debug() << "join-point of i: " << sausagei.endpoints[iEndpoint] << endl;
				debug() << "join-point of i: " << *(sausagei.points[sausagei.endpoints[iEndpoint]]) << endl;
				debug() << "non-join-point of j: " << sausagej.endpoints[(jEndpoint+1)%2] << endl;
				debug() << "non-join-point of j: " << *(sausagej.points[sausagej.endpoints[(jEndpoint+1)%2]]) << endl;;
				debug() << "Sausagei's endpoints are now at: " << *(sausagei.points[sausagei.endpoints[0]])
					<< " and " << *(sausagei.points[sausagei.endpoints[1]]) << endl;

				// Recurse, in case more sausages need joining, then exit
				join_endpoints(model);
				return;
			}
		}
	}
}

void Sausage::find_com(void){
	double x=0,y=0,z=0;
	for (vector<Point*>::iterator iter=points.begin(); iter!=points.end(); ++iter){
		Point* it=*iter;
		x+=it->x;
		y+=it->y;
		z+=it->z;
	}
	centre_of_mass[0] = x/points.size();
	centre_of_mass[1] = y/points.size();
	centre_of_mass[2] = z/points.size();
	info() <<  "Centre of mass " << centre_of_mass[0] << " " << centre_of_mass[1] << " " << centre_of_mass[2] << endl;
}

void Sausage::rotate_to_xy_plane(double** pointsArray){

	calculate_rotation_matrix();
	debug() << "Rotating to xy-plane" << endl;

	// loop over all coordinate points and rotate
	for(size_t i = 0; i != points.size(); i++) {
		double x = rotation_matrix(0)*pointsArray[i][0] + rotation_matrix(1)*pointsArray[i][1] + rotation_matrix(2)*pointsArray[i][2];
		double y = rotation_matrix(3)*pointsArray[i][0] + rotation_matrix(4)*pointsArray[i][1] + rotation_matrix(5)*pointsArray[i][2];
		double z = rotation_matrix(6)*pointsArray[i][0] + rotation_matrix(7)*pointsArray[i][1] + rotation_matrix(8)*pointsArray[i][2];
		pointsArray[i][0] = x;
		pointsArray[i][1] = y;
		pointsArray[i][2] = z;
	}
	return;
}

void Sausage::rotate_to_xy_plane(std::vector<Vector3d> pointsArray){

	calculate_rotation_matrix();
	debug() << "Rotating to xy-plane" << endl;

	// loop over all coordinate points and rotate
	/*
	for(size_t i = 0; i != points.size(); i++) {
		pointsArray[i][0] = rotation_matrix(0)*pointsArray[i][0] + rotation_matrix(1)*pointsArray[i][1] + rotation_matrix(2)*pointsArray[i][2];
		pointsArray[i][1] = rotation_matrix(3)*pointsArray[i][0] + rotation_matrix(4)*pointsArray[i][1] + rotation_matrix(5)*pointsArray[i][2];
		pointsArray[i][2] = rotation_matrix(6)*pointsArray[i][0] + rotation_matrix(7)*pointsArray[i][1] + rotation_matrix(8)*pointsArray[i][2];
	}
	*/
	for (vector<Vector3d>::iterator iter=pointsArray.begin(); iter!=pointsArray.end(); ++iter){
		double x = rotation_matrix(0)*iter->x + rotation_matrix(1)*iter->y + rotation_matrix(2)*iter->z;
		double y = rotation_matrix(3)*iter->x + rotation_matrix(4)*iter->y + rotation_matrix(5)*iter->z;
		double z = rotation_matrix(6)*iter->x + rotation_matrix(7)*iter->y + rotation_matrix(8)*iter->z;
		//iter = Vector3d (x,y,z);
		iter->x = x;
		iter->y = y;
		iter->z = z;
	}
	return;
}

// Calculate COM of all points contained within a sphere of radius 'radius' around point center
// radius is calculate if the input value is below R_min_sphere, radius is chosen so that 10 points lie within 
// if radius is above the threshold the value itself is used
void Sausage::calculate_com_sphere(Vector3d center, double radius, Vector3d &com){

    double R;
    int counter;
    Vector3d delta (0,0,0);

    //Calculate radius, so that sphere includes 10 points, unless radius is passed into the fn in which case it should be > R_min
    if (radius < params::R_min_sphere){
        R = params::R_min_sphere;
		counter =0;
		while (counter < 10) { //exit loop once R is sufficiently large to include at least 10 points
			counter = 0;
			// loop over all points
			for(size_t i = 0; i != points.size(); i++){
				delta.x = points[i]->x - center.x;
				delta.y = points[i]->y - center.y;
				delta.z = points[i]->z - center.z;
				if ( mag(delta) < R) counter++;
			}
			R += params::dR_sphere;    //gradually increase radius
		}
    }
    //set R to radius that was passed into the fn
    else{ R = radius;}

    // calcualte COM for all points within distance R of center point
    counter = 0;
	com.x= 0.0;
	com.y= 0.0;
	com.z= 0.0;
	for(size_t i = 0; i != points.size(); i++){
		delta.x = points[i]->x - center.x;
		delta.y = points[i]->y - center.y;
		delta.z = points[i]->z - center.z;
		// add all points within radius to work out com 
		if (mag(delta) < R) {
			com.x += points[i]->x;
			com.y += points[i]->y;
			com.z += points[i]->z;
			counter ++;
		}
	}
	if (counter != 0){
		com.x /= counter; com.y /= counter; com.z /= counter;}
	else {
		cerr << " Halfsphere algorithm failed. Sphere to calculate COM is empty! " << endl;
		exit(EXIT_FAILURE);}
}

void Sausage::sphere_tracking(){

    info() << "--------> Estimating sausage length using spheres... " << endl;

	Vector3d current_pos (0,0,0);                   //current grid position in sausage
    Vector3d com (0,0,0);                           //current com calculated      
    Vector3d dir (0,0,0),dir_hat (0,0,0);                   //current direction + unit vector direction
    Vector3d delta (0,0,0);

    double dist_max = 10000;
	int max_index;

	// first point is arbritary pick
	current_pos.x = points[0]->x;
	current_pos.y = points[0]->y;
	current_pos.z = points[0]->z;

	// calculate COM of sphere with radius R_gap around current_pos
	calculate_com_sphere(current_pos,params::R_gap,com);
	verbose() << "First COM for sphere tracking: " << com << std::endl;
    sphere_COMs.push_back(com);

	// pick 2nd point, pick point which is closest to be 2*R_gap away from first point
	for(size_t i = 0; i != points.size(); i++){
		delta.x = points[i]->x - current_pos.x;
		delta.y = points[i]->y - current_pos.y;
		delta.z = points[i]->z - current_pos.z;
		if ( fabs(mag(delta) - 2.0*params::R_gap) < dist_max){
			dist_max = fabs(mag(delta) - 2.0*params::R_gap);
			max_index = i;
		}
	}
    // move to 2nd point in sausage
    current_pos.x = points[max_index]->x;
	current_pos.y = points[max_index]->y;
	current_pos.z = points[max_index]->z;
	calculate_com_sphere(current_pos,params::R_gap,com);
	verbose() << "Second COM for sphere tracking: pos " << com << std::endl;
    sphere_COMs.push_back(com);

	double distsq_previous, distsq_preprevious;
	bool reached_start = false;
    int index = 2; //because first two points are already found

    // move along the disclination line step by step until we reach the start point
	while (reached_start == false){

		// work out moving direction from two previous coms
		dir = sphere_COMs[index-1] - sphere_COMs[index-2];
		dir_hat = norm (dir); 
		debug() << "Halfsphere tracking direction " << dir_hat <<std::endl;

		// find new current position in sausage by moving alon the direction dir 
		current_pos = sphere_COMs[index-1] + params::R_gap*dir_hat;
	
        //find sphere with radius R for current position that is large enough to include 10 points
        //then find COM of all included points
		calculate_com_sphere(current_pos,0.0,com);
		
		// check that we aren't reversing on ourselves, if so exit the program
		distsq_previous = distance(com,sphere_COMs[index-1]);
		distsq_preprevious = distance(com,sphere_COMs[index-2]);
        if (distsq_preprevious < distsq_previous){ 
            cerr << "The next point in sphere sausage tracking was closer to the preprevious point than the previous one. We are reversing back on ourselves." << endl;
        }
       
        //add COM found to final array of COMs if index > 5, because first few points will be slighly noisy
        sphere_COMs.push_back(com);
		
		// check whether we have reached our starting point once we have found the first few points
		if (distance(sphere_COMs[index],sphere_COMs[0]) < params::R_gap)
		{
			info() << "Halfsphere algorithm has successfully reached its starting point. \n" << endl;
			reached_start = true;
		}
		index++;

        // Exit program if sphere tracking takes too long.
		if ( index > 1000) {
			cerr << "Halfsphere algorithm hasn't reached its starting point after 1000 iterations. \n" << endl;
			exit(EXIT_FAILURE);
		}

	}//end of while loop over sausage

    //print all COMs of sphere
	verbose()<<"COMS sphere tracking:"<<std::endl;
	for (int k=0; k<sphere_COMs.size(); k++){
    	brief({2})<<"COM "<<sphere_COMs[k]<<std::endl;
    	verbose()<<"COM "<<sphere_COMs[k]<<std::endl;}

    return;
}

void Sausage::calculate_sausage_length_spheres(){

	// loop over all COMs found by sphere tracking to calculate length
	length = 0;
	for(int i = 0; i < sphere_COMs.size()-1; i++) { 
		length += distance(sphere_COMs[i],sphere_COMs[i+1]); 
        debug() << "Length "<<i<<" "<< length << endl;
    }
	length += distance(sphere_COMs[sphere_COMs.size()-1],sphere_COMs[0]);
    debug() << "Length last step"<< length << endl;

	info() << "The estimated length of the sausage using spheres is " << length << endl;
	brief({2}) << "Length "<<sausageID<<" "<< length << endl;
	return;
}

/** Check that all relevant sausages are closed loops
 * Add description
 */
void Sausage::flood_fill_closed_loops(const std::vector<Vector3d> colloidPos){

	debug() << "begin flood_fill_closed_loops" << endl << flush;
    
    // find centre of mass of this sausage
    find_com();
    Vector3d com_vec (0,0,0);
    com_vec.x = centre_of_mass[0];
    com_vec.y = centre_of_mass[1];
    com_vec.z = centre_of_mass[2];

    // check that COM isn't very close to the sausage
    Vector3d delta (0,0,0);
    for(size_t i = 0; i != points.size(); i++){
		delta.x = points[i]->x - com_vec.x;
		delta.y = points[i]->y - com_vec.y;
		delta.z = points[i]->z - com_vec.z;
		if ( fabs(mag(delta) < 3.0)){
            cerr << "FF closed loops, shows that COM is very close to sausage point." << endl;
			exit(EXIT_FAILURE);
		}
	}
    debug() << "Check that COM is at least 3.0 away from any sausge point."<< endl;

	// AB is vector joining colloids
	Vector3d AB = colloidPos[0] - colloidPos[1];
	double norm_AB = sqrt(dot(AB,AB));

	// Plane of best fit is of form alpha*x+beta*y=z
	double alpha=-plane_of_best_fit[0];
	double beta=-plane_of_best_fit[1];

	// From being parallel to PoBF, perpendicular to AB, we can determine: (This is an arbitrary pick.)
	Vector3d v (0,0,0);
	if (abs(AB.y+beta*AB.z) > params::epsilon){
		v.x=1; // No need to normalise
		v.y=(-alpha-AB.x)/(AB.y+beta*AB.z);
	}else{
		v.x=(-AB.y-beta*AB.z)/(AB.x+alpha); // = 0, more or less
		v.y=1.0/sqrt(1+beta*beta); // This pretty much normalises for free
	}
	v.z=alpha*v.x+beta*v.y;

	double norm_v = sqrt(dot(v,v));
	verbose()<<"Vector v (perp. AB & || PoBF & unit) is "<<v<<" of length "<<norm_v<<endl;

	// Loop over all points in sausage and add the ones within tubes around v starting in COM
	vector<Point*> region1,region2,region;
	for (vector<Point*>::iterator iter=points.begin(); iter!=points.end(); iter++){
		Point* it = *iter;
		Vector3d p (0,0,0); // xyz of this point
		p.x=it->x; p.y=it->y; p.z=it->z;

        Vector3d u = p - com_vec; // Vector from point to COM

		// Eq. of plane perpendicular to AB going through point p is (x,y,z).AB - p.AB
		// Distance from point q to this plane P is |Px*qx + Py*qy + Pz*qz + P0|/norm(Px,Py,Pz)
		Vector3d first_plane  = com_vec + 0.5*params::flood_fill_classify_slice_size*params::pixel_size*(AB/norm_AB);
		Vector3d second_plane = com_vec - 0.5*params::flood_fill_classify_slice_size*params::pixel_size*(AB/norm_AB);
		double first_distance  = abs( dot(AB, p) - dot(first_plane,  AB) ) / norm_AB;
		double second_distance = abs( dot(AB, p) - dot(second_plane, AB) ) / norm_AB;

		// Is it between the two planes? If so add to region 
		if (first_distance < params::flood_fill_classify_slice_size*params::pixel_size &&
		    second_distance < params::flood_fill_classify_slice_size*params::pixel_size ){
		    double projection = dot(u,v);
			if (projection>0){
				region1.push_back(it);
				verbose()<<it->x<<" "<<it->y<<" "<<it->z<<"TESTIN1"<<endl;
            }else{
      			region2.push_back(it);
				verbose()<<it->x<<" "<<it->y<<" "<<it->z<<"TESTIN2"<<endl;
            }
         }else{verbose()<<it->x<<" "<<it->y<<" "<<it->z<<"TESTOUT"<<endl;}
    }

	verbose()<<region1.size()<<" pixels in region 1"<<endl;
	verbose()<<region2.size()<<" pixels in region 2"<<endl;
	if (region1.size() == 0 && region2.size() == 0){
		cerr << "Both regions in FF closed loops check is empty, this is bad, exiting." << endl << flush;
		exit(EXIT_FAILURE);
	}

    //Pick larger region
    if (region1.size() >= region2.size() ) region = region1;
    verbose()<<region.size()<<" pixels in region chosen (larger region)"<<endl;

    //Find vector u which is parallel to PoBF, perpendicular to v
	Vector3d u (0,0,0);
	if (abs(v.y+beta*v.z) > params::epsilon){
		u.x=1; // No need to normalise
		u.y=(-alpha-v.x)/(v.y+beta*v.z);
	}else{
		u.x=(-v.y-beta*v.z)/(v.x+alpha); // = 0, more or less
		u.y=1.0/sqrt(1+beta*beta); // This pretty much normalises for free
	}
	u.z=alpha*u.x+beta*u.y;

	verbose()<<"Vector u (perp. v & || PoBF ) is "<<u<<endl;


	//Find com of region selected
	double x=0,y=0,z=0;
    Vector3d com_vec_reg (0,0,0); //com of region selected
   	for (vector<Point*>::iterator it=region.begin(); it!=region.end(); it++){
		x+=(*it)->x;
		y+=(*it)->y;
		z+=(*it)->z;
	}
	com_vec_reg.x = x/region.size();
	com_vec_reg.y = y/region.size();
	com_vec_reg.z = z/region.size();
	info() <<  "Centre of mass of region " << com_vec_reg << endl;
   	// Find neighbours of region
    vector<Point*> neighbours;
    vector<Point*> neighbours_otherside;
    Vector3d p (0,0,0); //Vector from COM of region to point p
	for (vector<Point*>::iterator iter=region.begin(); iter!=region.end(); iter++){
		Point* it=*iter;
		for (int iNeigh=0;iNeigh<6;iNeigh++){
			Point* neigh = it->neighbours[iNeigh];
    		if (neigh->isInASausage && !elementOf(region,neigh)){
                p.x = neigh->x; 
                p.y = neigh->y; 
                p.z = neigh->z; 
                //Only pick neighbour points along u away from com of region
                debug() << " p  "<<p<< endl;
                debug() << " p - com_vec_reg "<<p-com_vec_reg<< endl;
                debug() << "dot(u,p-com_vec_reg) "<<dot(u,p-com_vec_reg)<<endl;
                if (dot(u,p-com_vec_reg) > 0){
			        neighbours.push_back(it->neighbours[iNeigh]);
                }else{
                    neighbours_otherside.push_back(it->neighbours[iNeigh]);
                }
			}
		}
	}
	debug()<<"end of region"<<endl<<flush;

	// Get rid of duplicate neighbours
	debug()<<region.size()<<" points in this region, "<<neighbours.size()<<" neighbours inside the sausage."<<endl;
	sort(neighbours.begin(),neighbours.end());
	debug()<<"sorted neighbours"<<endl;
	vector<Point*>::iterator it = unique(neighbours.begin(),neighbours.end());
	neighbours.resize(distance(neighbours.begin(),it));
	debug()<<"Removed duplicates, now "<<neighbours.size()<<" neighbours."<<endl;
	debug()<<"Neighbours at:"<<endl;
	for (vector<Point*>::iterator it=neighbours.begin(); it!=neighbours.end(); it++){
		verbose()<<(*it)->x<<" "<<(*it)->y<<" "<<(*it)->z<<"NEIGH"<<endl;
	}
    verbose()<<neighbours.size()<<" pixels in neighbour region"<<endl;

    // Get rid of duplicate neighbours on the other side of region
	debug()<<region.size()<<" points in this region, "<<neighbours_otherside.size()<<" neighbours otherside inside the sausage."<<endl;
	sort(neighbours_otherside.begin(),neighbours_otherside.end());
	debug()<<"sorted neighbours otherside"<<endl;
	vector<Point*>::iterator it2 = unique(neighbours_otherside.begin(),neighbours_otherside.end());
	neighbours_otherside.resize(distance(neighbours_otherside.begin(),it2));
	debug()<<"Removed duplicates, now "<<neighbours_otherside.size()<<" neighbours otherside."<<endl;
	debug()<<"Neighbours at:"<<endl;
	for (vector<Point*>::iterator it=neighbours_otherside.begin(); it!=neighbours_otherside.end(); it++){
		verbose()<<(*it)->x<<" "<<(*it)->y<<" "<<(*it)->z<<"NEIGH2"<<endl;
	}
    verbose()<<neighbours_otherside.size()<<" pixels in neighbour otherside region"<<endl;


	// flood-fill from neighbours. If the point is in the region, don't add it to the flood-filled area.
	vector<Point*> stack=neighbours; // FILO stack, to be filled with contiguous points
	vector<Point*> visited; // Points we've already visited, and so should ignore
	while (stack.size()>0){
		//debug()<<"Another pixel on the stack"<<endl;
		Point* curr = stack.back();
        debug()<<"STACK "<< curr->x << " "<<curr->y<<" "<<curr->z<< endl;
		stack.pop_back();
		for (int iNeigh=0;iNeigh<6;iNeigh++){
			if (!curr->neighbours[iNeigh]) continue;
			Point* thisNeigh=curr->neighbours[iNeigh];
			if (elementOf(visited,thisNeigh)) continue;
				visited.push_back(thisNeigh);
				if (thisNeigh->isInASausage) {
				    //debug()<<"neigh: "<<iNeigh<<" address: "<<curr.neighbours[iNeigh]<<" sID: "<<curr.neighbours[iNeigh]->sausageID<<endl;
                    //check whether current point is inside the region
				    bool in_region=false;
				    if (elementOf(neighbours_otherside,thisNeigh)){
					    in_region=true;
					    verbose()<<"I'm a neighbour on the other side."<<endl;
                        //remove from neighbours_otherside list
                        neighbours_otherside.erase(std::remove(neighbours_otherside.begin(), neighbours_otherside.end(), thisNeigh), neighbours_otherside.end());
                        verbose()<<"REGION "<< thisNeigh->x << " "<<thisNeigh->y<<" "<<thisNeigh->z<< endl;
				    }
                    if (!in_region && !elementOf(stack,thisNeigh)) stack.push_back(thisNeigh);
                }
	    }
    }

    if (neighbours_otherside.size() != 0){
        verbose()<<neighbours_otherside.size()<<" pixels in neighbour otherside region that were undetected by FF"<<endl;
        cerr << "It seems there are open loops in the system. Exit program ...";
        exit(EXIT_FAILURE);
    }
    verbose()<<"Finished flood-fill closed loops."<<endl;
}

