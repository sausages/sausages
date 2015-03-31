#include <iostream>
#include <vector>
#include <numeric>
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
void flood_fill_separate(vector<Point> &allPoints, vector<Sausage> &allSausages){
	debug() << "begin flood_fill_separate" << endl << flush;

	int newSausageID=0;
	vector<size_t> stack; // Empty FILO stack, filled with indices of points to be coloured
	for (size_t firstPoint=0; firstPoint<allPoints.size(); firstPoint++){
		// Pick the first sausage not yet coloured
		if (allPoints[firstPoint].isInASausage && allPoints[firstPoint].sausageID==-1){
			verbose() << "Starting flood-fill of sausage #" << newSausageID << endl << flush;
			Sausage newSausage(newSausageID);
			bool writePoints = false;
			std::ofstream sausageFile;
			if (!params::sausage_filename.empty()){
				writePoints = true;
				sausageFile.open(params::sausage_filename+std::to_string(newSausageID));
			}
			stack.push_back(firstPoint);
			while (stack.size()>0){
				Point *curr = &(allPoints[stack.back()]);
				stack.pop_back();

				curr->sausageID=newSausageID;
				newSausage.points.push_back(curr);
				curr->sausagePointsIndex=newSausage.points.size()-1;
				if (writePoints) sausageFile << curr->x << ' ' << curr->y << ' ' << curr->z << endl;
				for (int iNeigh=0;iNeigh<6;iNeigh++){
					Point *neigh=curr->neighbours[iNeigh];
					if (neigh  && neigh->isInASausage && neigh->sausageID==-1 && // neigh is sausage-worthy but not yet sorted
					    //find(stack.begin(),stack.end(),neigh->allPointsIndex)==stack.end()){ // Don't add to the stack if it's already there
					    !elementOf(stack,neigh->allPointsIndex)){ // Don't add to the stack if it's already there
						stack.push_back(neigh->allPointsIndex) ;
					}
				}
			}
			allSausages.push_back(newSausage);
			debug()<<allSausages[newSausageID].points.size()<<" points in the sausages."<<endl<<flush;
			if (writePoints) sausageFile.close();
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
void Sausage::flood_fill_classify(const std::vector<vector3d> colloidPos){
	/* We need to find the vector parallel to the sausage's PoBF which is perpendicular to the line joining the two colloids.
	 * Then we follow this vector out from a colloid to find suitable points on the sausage.
	 * If we arbitrarity set the vector to be length 1 we can determine it from the colloid positions and the PoBF.
	 */

	debug() << "begin flood_fill_classify" << endl << flush;

	// AB is vector joining colloids
	vector3d AB = colloidPos[0] - colloidPos[1];
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
	vector3d v;
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
			vector3d p; // xyz of this point
			p.x=it->x; p.y=it->y; p.z=it->z;

			vector3d u = p - colloidPos[iColl]; // Vector from point to colloid

			// Eq. of plane perpendicular to AB going through point p is (x,y,z).AB - p.AB
			// Distance from point q to this plane P is |Px*qx + Py*qy + Pz*qz + P0|/norm(Px,Py,Pz)
			vector3d first_plane  = colloidPos[iColl] + 0.5*params::flood_fill_classify_slice_size*params::pixel_size*(AB/norm_AB);
			vector3d second_plane = colloidPos[iColl] - 0.5*params::flood_fill_classify_slice_size*params::pixel_size*(AB/norm_AB);
			double first_distance  = abs( dot(AB, p) - dot(first_plane,  AB) ) / norm_AB;
			double second_distance = abs( dot(AB, p) - dot(second_plane, AB) ) / norm_AB;

			// Is it between the two planes? If so, it is <= flood_fill_classify_slice_size from each plane
			if (first_distance < params::flood_fill_classify_slice_size*params::pixel_size &&
			   second_distance < params::flood_fill_classify_slice_size*params::pixel_size ){
				double projection = dot(u,v);
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

	if (below[0].size() * below[1].size() * above[0].size() * above[1].size() == 0){
		cerr << "One of the regions is empty, this is bad, exiting." << endl << flush;
		exit(EXIT_FAILURE);
	}

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
				//if (neigh->isInASausage && find(region.begin(),region.end(),neigh)==region.end()){
				if (neigh->isInASausage && !elementOf(region,neigh)){
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
				//if (find(visited.begin(),visited.end(),thisNeigh)!=visited.end()) continue;
				if (elementOf(visited,thisNeigh)) continue;
				visited.push_back(thisNeigh);
				//debug()<<"neigh: "<<iNeigh<<" address: "<<curr.neighbours[iNeigh]<<" sID: "<<curr.neighbours[iNeigh]->sausageID<<endl;
				if (thisNeigh->isInASausage) {
					bool in_a_region=false;
					// For all of regions...
					for (int jColl=0;jColl<2;jColl++){ for (int inner_aboveOrBelow=0;inner_aboveOrBelow<2;inner_aboveOrBelow++){
						vector<Point*> inner_region = inner_aboveOrBelow ? above[jColl] : below[jColl];
						// If the neighbour is in the region, make a note but don't add to stack
						//if (find(inner_region.begin(),inner_region.end(),thisNeigh)!=inner_region.end()){
						if (elementOf(inner_region,thisNeigh)){
							in_a_region=true;
							adjacency[2*iColl+aboveOrBelow][2*jColl+inner_aboveOrBelow]=1;
							//debug()<<"I'm in a region"<<endl;
							//debug()<<"From "<<2*iColl+aboveOrBelow<<" to "<<2*jColl+inner_aboveOrBelow<<endl;
						}
					}}
					//if (!in_a_region && find(stack.begin(),stack.end(),thisNeigh)==stack.end()) stack.push_back(thisNeigh);
					if (!in_a_region && !elementOf(stack,thisNeigh)) stack.push_back(thisNeigh);
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
		brief({1}) << "Figure of omega structure." << endl;
		brief({1}) << "3" << endl;
	} else if (adjacency[0][3] && adjacency[1][2]){
		info()<<"System appears to be twisted (Class 2)"<<endl;
		brief({1}) << "Figure of eight structure." << endl;
		brief({1}) << "2" << endl;
	
	} else {
		brief({1}) << "Structure undefined." << endl;
		brief({1}) << "0" << endl;
		error()<<"Unexpected system linkage, this should never happen."<<endl;
		exit(EXIT_FAILURE);
	}
}

void Sausage::find_endpoints(void){

	debug() << "begin find_endpoints" << endl << flush;

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

				// At each step we move along the dimension which takes us fastest to target.
				while (curr->x!=target->x || curr->y!=target->y || curr->z!=target->z){
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
					double dx = abs(curr->x - target->x);
					double dy = abs(curr->y - target->y);
					double dz = abs(curr->z - target->z);

					if (dx>dy && dx>dz){
						curr = (curr->x > target->x) ? curr->left : curr->right;
					} else if (dy>dz) {
						curr = (curr->y > target->y) ? curr->down : curr->up;
					} else {
						curr = (curr->z > target->z) ? curr->back : curr->forward;
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

// Calculate COM of all points contained within a sphere of radius 'radius' around point center
// radius is calculate if the input value is below R_min_sphere, radius is chosen so that 10 points lie within 
// if radius is above the threshold the value itself is used
void Sausage::calculate_com_sphere(vector3d center, double radius, vector3d &com){

    double R;
    int counter;
    vector3d delta;

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
				if ( magnitude(delta) < R) counter++;
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
		if (magnitude(delta) < R) {
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

	vector3d current_pos;                   //current grid position in sausage
    vector3d com;                           //current com calculated      
    vector3d dir,dir_hat;                   //current direction + unit vector direction
    vector3d delta;

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
		if ( fabs(magnitude(delta) - 2.0*params::R_gap) < dist_max){
			dist_max = fabs(magnitude(delta) - 2.0*params::R_gap);
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
       
        //add COM found to final array of COMs
        sphere_COMs.push_back(com);
		
		// check whether we have reached our starting point
		if ( distance(sphere_COMs[index],sphere_COMs[0]) < params::R_gap)
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
    	verbose()<<sphere_COMs[k]<<std::endl;}

    return;
}

void Sausage::calculate_sausage_length_spheres(){

	// loop over all COMs found by sphere tracking to calculate length
	length = 0;
	for(int i = 0; i < sphere_COMs.size(); i++) { 
		length += distance(sphere_COMs[i],sphere_COMs[i+1]); 
    }
	length += distance(sphere_COMs[sphere_COMs.size()-1],sphere_COMs[0]);

	info() << "The estimated length of the sausage using spheres is " << length << endl;
	return;
}


