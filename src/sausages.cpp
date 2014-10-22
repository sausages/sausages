#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include <limits>
#include "io.h"
#include "point.h"
#include "sausages.h"
#include "params.h"

using namespace std;

Sausage::Sausage(int ID){
	sausageID=ID;
}

/** Points are in a sausge if their cl value is < threshold
 */
int threshold(vector<Point> &allPoints){
	int num_below_threshold=0;
	vector<Point>::iterator it;
	for (it=allPoints.begin(); it!=allPoints.end(); ++it){
		if (it->cl < params::threshold){
			it->isInASausage=true;
			num_below_threshold++;
		}
	}
	return num_below_threshold;
}

/** Ordinary Flood-fill algorithm
 * 1) Find the first sausage-worthy point not yet sorted into a sausage
 * 2) Add it to a numbered sausage, and add all its unsorted neighbours to a stack
 * 3) While there are points in the stack, take the last one, add it to
 *     the sausage and its unsorted neighbours to the stack
 * 4) Eventually you'll have a continuous sausage, start with the next one.
 */
void flood_fill_separate(vector<Point> &allPoints, vector<Sausage> &allSausages){

	int newSausageID=0;
	vector<size_t> stack; // Empty FILO stack, filled with indices of points to be coloured
	for (size_t firstPoint=0; firstPoint<allPoints.size(); firstPoint++){
		// Pick the first sausage not yet coloured
		if (allPoints[firstPoint].isInASausage && allPoints[firstPoint].sausageID==-1){
			verbose() << "Starting flood-fill of sausage #" << newSausageID << endl;
			Sausage newSausage(newSausageID);
			stack.push_back(firstPoint);
			while (stack.size()>0){
				Point *curr = &(allPoints[stack.back()]);
				stack.pop_back();

				curr->sausageID=newSausageID;
				newSausage.points.push_back(curr);
				curr->sausagePointsIndex=newSausage.points.size()-1;
				for (int iNeigh=0;iNeigh<6;iNeigh++){
					Point *neigh=curr->neighbours[iNeigh];
					if (neigh  && neigh->isInASausage && neigh->sausageID==-1 && // neigh is sausage-worthy but not yet sorted
					    find(stack.begin(),stack.end(),neigh->allPointsIndex)==stack.end()){ // Don't add to the stack if it's already there
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
 * Class 1: If we reach only the other colloid's point of the same side (i.e. above-1 -> above-2) then we
 *  have a simple ring with no twists, such as a theta configuration.
 *  Also classified as class-1 will be 2nd-loop systems where the 2nd loop is touching the 1st on only one side.
 *
 * Class 2: If we end up only on the other side of the opposite colloid (above-1 -> below-2)
 *  then we have a figure-of-eight system.
 *
 * Class 3: If we reach more than one other point we have junctions on both sides of the sausage, and are
 *  insufficiently resolved.
 *
 */
void Sausage::flood_fill_classify(void){
	/* We need to find the vector parallel to the sausage's PoBF which is perpendicular to the line joining the two colloids.
	 * Then we follow this vector out from a colloid to find suitable points on the sausage.
	 * If we arbitrarity set the vector to be length 1 we can determine it from the colloid positions and the PoBF.
	 */
	double v[3];
	double AB[3];
	double alpha,beta;

	// AB is vector joining colloids
	AB[0]=params::colloids[0][0]-params::colloids[1][0];
	AB[1]=params::colloids[0][1]-params::colloids[1][1];
	AB[2]=params::colloids[0][2]-params::colloids[1][2];
	double norm_AB=sqrt(inner_product(AB,AB+3,AB,0.0));
	debug()<<"Vector AB is {"<<AB[0]<<","<<AB[1]<<","<<AB[2]<<"} of length "<<norm_AB<<endl;

	// Plane of best fit is of form alpha*x+beta*y=z
	alpha=-plane_of_best_fit[0];
	beta=-plane_of_best_fit[1];

	// From being parallel to PoBF, perpendicular to AB, we can determine:
	if (abs(AB[1]+beta*AB[2]) > params::epsilon){
		v[0]=1; // No need to normalise
		v[1]=(-alpha-AB[0])/(AB[1]+beta*AB[2]);
	}else{
		v[0]=(-AB[1]-beta*AB[2])/(AB[0]+alpha); // = 0, more or less
		v[1]=1.0/sqrt(1+beta*beta); // This pretty much normalises for free
	}
	v[2]=alpha*v[0]+beta*v[1];

	double norm_v=sqrt(inner_product(v,v+3,v,0.0));
	debug()<<"Vector v (perp. AB & || PoBF & unit) is {"<<v[0]<<","<<v[1]<<","<<v[2]<<"} of length "<<norm_v<<endl;

	// For each point in the sausage, is it between two parallel planes extended from beside a colloid, perpendicular to the line between the colloids?
	// If so, add it to a region dependant on whether it is v-wards or anti-v-wards
	vector<Point*> above[2], below[2];
	for (vector<Point*>::iterator iter=points.begin(); iter!=points.end(); iter++){
		Point* it = *iter;
		for (int iColl=0; iColl<2; iColl++){
			double p[3]; // xyz of this point
			p[0]=it->x; p[1]=it->y; p[2]=it->z;

			double u[3]; // Vector from point to colloid
			u[0]=it->x - params::colloids[iColl][0];
			u[1]=it->y - params::colloids[iColl][1];
			u[2]=it->z - params::colloids[iColl][2];
			//debug()<<"Vector u is {"<<u[0]<<","<<u[1]<<","<<u[2]<<"} of length "<<norm_u<<endl;

			// Eq. of plane perpendicular to AB going through point p is (x,y,z).AB - p.AB
			// Distance from point q to this plane P is |Px*qx + Py*qy + Pz*qz + P0|/norm(Px,Py,Pz)
			double first_plane[3], second_plane[3];
			for (int i=0; i<3; i++){
				first_plane[i]  = params::colloids[iColl][i] + 0.5*params::flood_fill_classify_slice_size*params::pixel_size*(AB[i]/norm_AB);
				second_plane[i] = params::colloids[iColl][i] - 0.5*params::flood_fill_classify_slice_size*params::pixel_size*(AB[i]/norm_AB);
			}
			double first_distance  = abs( inner_product(AB, AB+3, p, 0.0) - inner_product(first_plane,first_plane+3,AB,0.0) ) / norm_AB;
			double second_distance = abs( inner_product(AB, AB+3, p, 0.0) - inner_product(second_plane,second_plane+3,AB,0.0) ) / norm_AB;
			//debug()<<"distances: "<<first_distance<<","<<second_distance<<endl;

			// Is it between the two planes? If so, it is <= flood_fill_classify_slice_size from each plane
			if (first_distance < params::flood_fill_classify_slice_size*params::pixel_size &&
			   second_distance < params::flood_fill_classify_slice_size*params::pixel_size ){
				double projection = inner_product(u,u+3,v,0.0);
				if (projection>0){
					above[iColl].push_back(it);
					debug()<<it->x<<","<<it->y<<","<<it->z<<" added to region above colloid "<<iColl<<endl;
				} else {
					below[iColl].push_back(it);
					debug()<<it->x<<","<<it->y<<","<<it->z<<" added to region below colloid "<<iColl<<endl;
				}
			}
		}
	}

	verbose()<<above[0].size()<<" pixels in region above first colloid"<<endl;
	verbose()<<above[1].size()<<" pixels in region above second colloid"<<endl;
	verbose()<<below[0].size()<<" pixels in region below first colloid"<<endl;
	verbose()<<below[1].size()<<" pixels in region below second colloid"<<endl;

	/* Debug: print out all points in sausage, and points in regions */
	/*
	cout << "XXXX" << endl;
	for (vector<Point>::iterator it=points.begin(); it!=points.end(); it++){ cout << 0 <<","<<it->x<<","<<it->y<<","<<it->z<<endl; }
	for (vector<Point>::iterator it=above[0].begin(); it!=above[0].end(); it++){ cout << 1 <<","<<it->x<<","<<it->y<<","<<it->z<<endl; }
	for (vector<Point>::iterator it=below[0].begin(); it!=below[0].end(); it++){ cout << 2 <<","<<it->x<<","<<it->y<<","<<it->z<<endl; }
	for (vector<Point>::iterator it=above[1].begin(); it!=above[1].end(); it++){ cout << 3 <<","<<it->x<<","<<it->y<<","<<it->z<<endl; }
	for (vector<Point>::iterator it=below[1].begin(); it!=below[1].end(); it++){ cout << 4 <<","<<it->x<<","<<it->y<<","<<it->z<<endl; }
	cout << "XXXX" << endl;
	*/

	/** Adjacency matrix of region
	 * 4-by-4 adjacency matrix. Each region is adjacent to itself.
	 * Regions are numbered:  0) below-0  1) above-0  2) below-1  3) above-1
	 * i.e. region is number 2*iColl+aboveOrBelow
	 */
	int adjacency[4][4] = {}; // initalise to all zero
	// For each region
	for (int iColl=0;iColl<2;iColl++){ for (int aboveOrBelow=0;aboveOrBelow<2;aboveOrBelow++){
		debug()<<"in region, colloid "<<iColl<<" aboveOrBelow "<<aboveOrBelow<<endl<<flush;
		vector<Point*> region = aboveOrBelow ? above[iColl] : below[iColl];
		int otherColl=(iColl+1)%2;

		// find all neighbours of the region which are inside this sausage and closer to the /other/ colloid
		vector<Point*> neighbours;
		for (vector<Point*>::iterator iter=region.begin(); iter!=region.end(); iter++){
			Point* it=*iter;
			double distanceToOtherColl = pow(it->x - params::colloids[otherColl][0],2) +
							pow(it->y - params::colloids[otherColl][1],2) +
							pow(it->z - params::colloids[otherColl][2],2) ;
			for (int iNeigh=0;iNeigh<6;iNeigh++){
				Point* neigh = it->neighbours[iNeigh];
				if (neigh->isInASausage && find(region.begin(),region.end(),neigh)==region.end()){
					double neighDistanceToOtherColl = pow(neigh->x - params::colloids[otherColl][0],2) +
									pow(neigh->y - params::colloids[otherColl][1],2) +
									pow(neigh->z - params::colloids[otherColl][2],2) ;
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
				if (find(visited.begin(),visited.end(),thisNeigh)!=visited.end()) continue;
				visited.push_back(thisNeigh);
				//debug()<<"neigh: "<<iNeigh<<" address: "<<curr.neighbours[iNeigh]<<" sID: "<<curr.neighbours[iNeigh]->sausageID<<endl;
				if (thisNeigh->isInASausage) {
					bool in_a_region=false;
					// For all of regions...
					for (int jColl=0;jColl<2;jColl++){ for (int inner_aboveOrBelow=0;inner_aboveOrBelow<2;inner_aboveOrBelow++){
						vector<Point*> inner_region = inner_aboveOrBelow ? above[jColl] : below[jColl];
						// If the neighbour is in the region, make a note but don't add to stack
						if (find(inner_region.begin(),inner_region.end(),thisNeigh)!=inner_region.end()){
							in_a_region=true;
							adjacency[2*iColl+aboveOrBelow][2*jColl+inner_aboveOrBelow]=1;
							//debug()<<"I'm in a region"<<endl;
							//debug()<<"From "<<2*iColl+aboveOrBelow<<" to "<<2*jColl+inner_aboveOrBelow<<endl;
						}
					}}
					if (!in_a_region && find(stack.begin(),stack.end(),thisNeigh)==stack.end()) stack.push_back(thisNeigh);
				}
			}
		}
		verbose()<<"Finished flood-fill from region "<<(aboveOrBelow?"above":"below")<<" colloid "<<iColl<<endl;

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
		info()<<"System appears to be untwisted (Class 1)"<<endl;
	} else if (adjacency[0][3] && adjacency[1][2]){
		info()<<"System appears to be twisted (Class 2)"<<endl;
	} else {
		error()<<"Unexpected system linkage, this should never happen."<<endl;
		exit(EXIT_FAILURE);
	}
}

void Sausage::find_endpoints(void){
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
	debug()<<"Arbitrary initial distance: "<<dist[10][15]<<endl;

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

	debug()<<"Endpoint 0 is #"<<endpoints[0]<<" at: "<<points[endpoints[0]]->x<<","<<points[endpoints[0]]->y<<","<<points[endpoints[0]]->z<<","<<endl;
	debug()<<"Endpoint 1 is #"<<endpoints[1]<<" at: "<<points[endpoints[1]]->x<<","<<points[endpoints[1]]->y<<","<<points[endpoints[1]]->z<<","<<endl;
	debug()<<"Minimum distance along the sausage between them is: "<<max_min_dist<<endl;

	for (size_t i=0; i<points.size(); i++){
		delete [] dist[i];
	}
	delete [] dist;
}

void join_endpoints(vector<Sausage> &allSausages, vector<int> &relevantSausages){
	// Work out distance between endpoints
	double dist[2*relevantSausages.size()][2*relevantSausages.size()];
	for (size_t i=0; i<relevantSausages.size(); i++){
		Sausage si = allSausages[relevantSausages[i]];
		for (size_t iEndpoint=0; iEndpoint<2; iEndpoint++){
			for (size_t j=i; j<relevantSausages.size(); j++){ // Only need half, by symmetry
				Sausage sj = allSausages[relevantSausages[j]];
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
	bool another_endpoint;
	for (size_t i=0; i<relevantSausages.size(); i++){
		for (size_t iEndpoint=0; iEndpoint<2; iEndpoint++){
			another_endpoint=false;
			int num_nearby_endpoints=0;
			for (size_t j=i; j<relevantSausages.size(); j++){ // Only need half, by symmetry
				for (size_t jEndpoint=0; jEndpoint<2; jEndpoint++){
					if (i==j && iEndpoint>=jEndpoint) continue; // More symmetry
					if (dist[2*i+iEndpoint][2*j+jEndpoint]<params::max_sausage_gap_size*params::pixel_size){
						debug()<<"Picked out "<<j<<":"<<jEndpoint<<" and "<<i<<":"<<iEndpoint<<endl;
						another_endpoint=true;
						other_endpoint=2*j+jEndpoint;
					}
					if (dist[2*i+iEndpoint][2*j+jEndpoint]<params::min_sausage_gap_next_neigh_distance*params::pixel_size){
						num_nearby_endpoints++;
					}

				}
			}
			if (num_nearby_endpoints>1){
				error()<<"Too many endpoints near endpoint "<<iEndpoint<<" of sausage "<<relevantSausages[i]<<endl;
				exit(EXIT_FAILURE);
			}
			if (another_endpoint){
				// Join endpoints i:iEndpoint and j:jEndpoint
				// To do this:
				//  - Add sausagej.points to sausagei, changing their sausageID and sausagePointsIndex (unless sausagej==sausagei)
				//  - Recalculate sausagei's endpoints
				//  - Remove sausagej from relevantSausages (unless sausagej==sausagei)

				size_t j=other_endpoint/2;
				size_t jEndpoint=other_endpoint%2;
				Sausage *sausagei=&allSausages[relevantSausages[i]];
				Sausage *sausagej=&allSausages[relevantSausages[j]];

				//  - Add sausagej.points to sausagei, changing their sausageID and sausagePointsIndex (unless sausagej==sausagei)
				if (sausagej!=sausagei){
					for (vector<Point*>::iterator iter=sausagej->points.begin(); iter!=sausagej->points.end(); ++iter){
						(*iter)->sausageID=relevantSausages[i];
						(*iter)->sausagePointsIndex+=sausagei->points.size();
					}
					debug()<<"sausagei size: "<<sausagei->points.size()<<endl;
					debug()<<"sausagej size: "<<sausagej->points.size()<<endl;
					sausagei->points.insert(sausagei->points.end(),sausagej->points.begin(),sausagej->points.end());
					debug()<<"new sausagei size: "<<sausagei->points.size()<<endl;
				}

				// Artificially add points between nearby endpoints
				Point *curr = sausagei->points[sausagei->endpoints[iEndpoint]];
				Point *target = sausagej->points[sausagej->endpoints[jEndpoint]];
				// Draw line along x
				while (curr->x!=target->x){
					if (curr->sausageID != sausagei->sausageID && curr->sausageID != sausagej->sausageID && curr->isInASausage){
						error()<<"Looks like the gap spans a third sausage, this is not good, aborting."<<endl;
						exit(EXIT_FAILURE);
					}
					if (curr->sausageID!=sausagei->sausageID){
						curr->isInASausage=true;
						curr->sausageID=sausagei->sausageID;
						sausagei->points.push_back(curr);
						curr->sausagePointsIndex=sausagei->points.size()-1;
					}
					if (curr->x > target->x){
						curr=curr->left;
					} else {
						curr=curr->right;
					}
				}
				// Draw line along y
				while (curr->y!=target->y){
					if (curr->sausageID != sausagei->sausageID && curr->sausageID != sausagej->sausageID && curr->isInASausage){
						error()<<"Looks like the gap spans a third sausage, this is not good, aborting."<<endl;
						exit(EXIT_FAILURE);
					}
					if (curr->sausageID!=sausagei->sausageID){
						curr->isInASausage=true;
						curr->sausageID=sausagei->sausageID;
						sausagei->points.push_back(curr);
						curr->sausagePointsIndex=sausagei->points.size()-1;
					}
					if (curr->y > target->y){
						curr=curr->down;
					} else {
						curr=curr->up;
					}
				}
				// Draw line along z
				while (curr->z!=target->z){
					if (curr->sausageID != sausagei->sausageID && curr->sausageID != sausagej->sausageID && curr->isInASausage){
						error()<<"Looks like the gap spans a third sausage, this is not good, aborting."<<endl;
						exit(EXIT_FAILURE);
					}
					if (curr->sausageID!=sausagei->sausageID){
						curr->isInASausage=true;
						curr->sausageID=sausagei->sausageID;
						sausagei->points.push_back(curr);
						curr->sausagePointsIndex=sausagei->points.size()-1;
					}
					if (curr->z > target->z){
						curr=curr->back;
					} else {
						curr=curr->forward;
					}
				}

				// Sausage has changed, need to recalculate endpoints
				sausagei->find_endpoints();

				//  - Remove sausagej from relevantSausages if it's been merged into sausagei,
				//    but check we haven't just closed a gap in sausagei (i.e. joined i:0 and i:0) 
				if (sausagei!=sausagej) relevantSausages.erase(relevantSausages.begin()+j); // Awkward, but erase() requires an iterator

				// Recurse, in case more sausages need joining, then exit
				join_endpoints(allSausages,relevantSausages);
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

	debug() << "Rotating to xy-plane" << endl;

	// loop over all coordinate points and rotate
	for(size_t i = 0; i != points.size(); i++) {
		pointsArray[i][0] = rotation_matrix(0)*pointsArray[i][0] + rotation_matrix(1)*pointsArray[i][1] + rotation_matrix(2)*pointsArray[i][2];
		pointsArray[i][1] = rotation_matrix(3)*pointsArray[i][0] + rotation_matrix(4)*pointsArray[i][1] + rotation_matrix(5)*pointsArray[i][2];
		pointsArray[i][2] = rotation_matrix(6)*pointsArray[i][0] + rotation_matrix(7)*pointsArray[i][1] + rotation_matrix(8)*pointsArray[i][2];
	}
	return;
}

void Sausage::calculate_sausage_length_pobf_coms(double **pos_coms_pobf_slice){

	// loop over all COM of slices and calculate length
	length = 0;
	for(int i = 0; i != nSlices-1; i++) {
		length += sqrt( pow(pos_coms_pobf_slice[i][0]-pos_coms_pobf_slice[i+1][0],2)+
				pow(pos_coms_pobf_slice[i][1]-pos_coms_pobf_slice[i+1][1],2)+
				pow(pos_coms_pobf_slice[i][2]-pos_coms_pobf_slice[i+1][2],2));
	}
	length += sqrt( pow(pos_coms_pobf_slice[nSlices-1][0]-pos_coms_pobf_slice[0][0],2)+
			pow(pos_coms_pobf_slice[nSlices-1][1]-pos_coms_pobf_slice[0][1],2)+
			pow(pos_coms_pobf_slice[nSlices-1][2]-pos_coms_pobf_slice[0][2],2));

	info() << "The estimated length of the sausage using pobf is " << length << endl;
	return;
}

void Sausage::estimate_sausage_length_using_pobf(){

	info() << "--------> Estimating sausage length using pobf slice... " << endl;

	// Make array of points, shifted so that COM is in origin
	double **rotatedPoints = new double*[points.size()];
	for (size_t i=0; i<points.size(); i++){
		rotatedPoints[i] = new double[3];
		rotatedPoints[i][0]=points[i]->x - centre_of_mass[0];
		rotatedPoints[i][1]=points[i]->y - centre_of_mass[1];
		rotatedPoints[i][2]=points[i]->z - centre_of_mass[2];
	}

	// calculate rotation matrix and its inverse, so that plane of best fit projects onto xy-plane
	calculate_rotation_matrix();
	rotate_to_xy_plane(rotatedPoints);

	// Find number of slices
	nSlices = points.size()/params::points_per_slice;
	info() << "Sausage was divided in " << nSlices << " slices." << endl;

	double **pos_coms_pobf_slice = new double*[nSlices];
	for (int i=0; i<nSlices; i++){
		pos_coms_pobf_slice[i] = new double[3];
	}
	double *slice_counter = new double[nSlices];


	// find what slice we are in (loop over all points)
	double theta;
	int sliceId;
	for(size_t i = 0; i != points.size(); i++) {

		theta=atan2(rotatedPoints[i][1],rotatedPoints[i][0])+M_PI;
		sliceId = (int) nSlices*theta/(2.0*M_PI);
		if (sliceId < 0 || sliceId >= nSlices ) {
			cerr << "SliceId is out of bound" << endl;
			exit(EXIT_FAILURE);
		}

		//add x,y,z coordinates to that slice
		pos_coms_pobf_slice[sliceId][0]+=rotatedPoints[i][0];
		pos_coms_pobf_slice[sliceId][1]+=rotatedPoints[i][1];
		pos_coms_pobf_slice[sliceId][2]+=rotatedPoints[i][2];
		slice_counter[sliceId]++;
	}

	// work out centre of mass for each slice
	for(int k = 0; k != nSlices; k++) {
		if (slice_counter[k] == 0) {
			cerr << "Slice is empty!!" << endl;
			exit(EXIT_FAILURE);}
		else{
			pos_coms_pobf_slice[k][0] = pos_coms_pobf_slice[k][0]/slice_counter[k];
			pos_coms_pobf_slice[k][1] = pos_coms_pobf_slice[k][1]/slice_counter[k];
			pos_coms_pobf_slice[k][2] = pos_coms_pobf_slice[k][2]/slice_counter[k];

			debug() << slice_counter[k] << " COM " << pos_coms_pobf_slice[k][0] << " " << pos_coms_pobf_slice[k][1] << " " << pos_coms_pobf_slice[k][2] << endl;}
	}

	// convert COM's of slice back to initial frame
	//rotate_from_xy_plane(pos_coms_pobf_slice);

	// calculate length
	calculate_sausage_length_pobf_coms(pos_coms_pobf_slice);

	// clean up
	debug() << "deleting rotated_points" << endl;
	for (size_t i=0; i<points.size(); i++){
		delete rotatedPoints[i];
	}
	delete rotatedPoints;

	debug() << "deleting pos_coms_pobf_slice" << endl;
	for (int i=0; i<nSlices; i++){
		delete pos_coms_pobf_slice[i];
	}
	delete pos_coms_pobf_slice;

	return;
}

void Sausage::calculate_sausage_length_halfsphere_coms(double **pos_coms_halfsphere){

	// loop over all COM of slices and calculate length
/*	length = 0;
	for(int i = 0; i != nSlices-1; i++) {
		length += sqrt( pow(pos_coms_halfsphere[i][0]-pos_coms_halfsphere[i+1][0],2)+
				pow(pos_coms_halfsphere[i][1]-pos_coms_halfsphere[i+1][1],2)+
				pow(pos_coms_halfsphere[i][2]-pos_coms_halfsphere[i+1][2],2));
	}
	length += sqrt( pow(pos_coms_halfsphere[nSlices-1][0]-pos_coms_halfsphere[0][0],2)+
			pow(pos_coms_halfsphere[nSlices-1][1]-pos_coms_halfsphere[0][1],2)+
			pow(pos_coms_halfsphere[nSlices-1][2]-pos_coms_halfsphere[0][2],2));

	info() << "The estimated length of the sausage using halfspheres is " << length << endl;*/
	return;
}

void Sausage::estimate_sausage_length_using_halfsphere(){

/*	info() << "--------> Estimating sausage length using halfspheres... " << endl;

	// Make array of points, shifted so that COM is in origin
	double **rotatedPoints = new double*[points.size()];
	for (size_t i=0; i<points.size(); i++){
		rotatedPoints[i] = new double[3];
		rotatedPoints[i][0]=points[i]->x - centre_of_mass[0];
		rotatedPoints[i][1]=points[i]->y - centre_of_mass[1];
		rotatedPoints[i][2]=points[i]->z - centre_of_mass[2];
	}

	// calculate rotation matrix and its inverse, so that plane of best fit projects onto xy-plane
	calculate_rotation_matrix();
	rotate_to_xy_plane(rotatedPoints);

	// Find number of slices
	nSlices = points.size()/params::points_per_slice;
	info() << "Sausage was divided in " << nSlices << " slices." << endl;

	double **pos_coms_halfsphere = new double*[nSlices];
	for (int i=0; i<nSlices; i++){
		pos_coms_halfsphere[i] = new double[3];
	}
	double *slice_counter = new double[nSlices];


	// find what slice we are in (loop over all points)
	double theta;
	int sliceId;
	for(size_t i = 0; i != points.size(); i++) {

		theta=atan2(rotatedPoints[i][1],rotatedPoints[i][0])+M_PI;
		sliceId = (int) nSlices*theta/(2.0*M_PI);
		if (sliceId < 0 || sliceId >= nSlices ) {
			cerr << "SliceId is out of bound" << endl;
			exit(EXIT_FAILURE);
		}

		//add x,y,z coordinates to that slice
		pos_coms_halfsphere[sliceId][0]+=rotatedPoints[i][0];
		pos_coms_halfsphere[sliceId][1]+=rotatedPoints[i][1];
		pos_coms_halfsphere[sliceId][2]+=rotatedPoints[i][2];
		slice_counter[sliceId]++;
	}

	// work out centre of mass for each slice
	for(int k = 0; k != nSlices; k++) {
		if (slice_counter[k] == 0) {
			cerr << "Slice is empty!!" << endl;
			exit(EXIT_FAILURE);}
		else{
			pos_coms_halfsphere[k][0] = pos_coms_halfsphere[k][0]/slice_counter[k];
			pos_coms_halfsphere[k][1] = pos_coms_halfsphere[k][1]/slice_counter[k];
			pos_coms_halfsphere[k][2] = pos_coms_halfsphere[k][2]/slice_counter[k];

			debug() << slice_counter[k] << " COM " << pos_coms_halfsphere[k][0] << " " << pos_coms_halfsphere[k][1] << " " << pos_coms_halfsphere[k][2] << endl;}
	}

	// convert COM's of slice back to initial frame
	//rotate_from_xy_plane(pos_coms_halfsphere);

	// calculate length
	calculate_sausage_length_halfsphere(pos_coms_halfsphere);

	// clean up
	debug() << "deleting rotated_points" << endl;
	for (size_t i=0; i<points.size(); i++){
		delete rotatedPoints[i];
	}
	delete rotatedPoints;

	debug() << "deleting pos_coms_halfsphere" << endl;
	for (int i=0; i<nSlices; i++){
		delete pos_coms_halfsphere[i];
	}
	delete pos_coms_halfsphere;*/

	return;
}



