#include <iostream>
#include <vector>
#include <stdexcept>
#include "point.h"
#include "find_sausages.h"
#include "main.h"

using namespace std;

/** Sets the sausageID of each point.
 * Set to 1 if the cl value is below the threshold (in a sausage),
 * 0 otherwise (not in a sausage)
 */
void threshold(vector<Point> &allPoints, double threshold_level){
	vector<Point>::iterator it;
	for (it=allPoints.begin(); it!=allPoints.end(); ++it){
		if (it->cl < threshold_level){
			it->sausageID=1;
		} else {
			it->sausageID=0;
		}
	}
}

/** Returns a vector where the nth element contains
 * the number of points where sausageID==n
 */
vector<int> count_sausages(const vector<Point> &allPoints){
	vector<Point>::const_iterator it;
	vector<int> sausage_count;
	for (it=allPoints.begin(); it!=allPoints.end(); ++it){
		try{
			sausage_count.at(it->sausageID)++;
		}catch(out_of_range o){
			sausage_count.resize(it->sausageID+1);
			sausage_count.at(it->sausageID)++;
		}
	}
	return sausage_count;
}

/** Ordinary Flood-fill algorithm
 * 1) Find the first point in an unnumbered sausage
 * 2) Add it to a numbered sausage, and add all its unsorted neighbours to a stack
 * 3) While there are points in the stack, take the last one, add it to
 *     the sausage and its unsorted neighbours to the stack
 * 4) Eventually you'll have a continuous sausage, start with the next one.
 */
void flood_fill(vector<Point> &allPoints){
	int newSausage=2; // Start at 2 (0 is no-sausage and 1 is unsorted)
	vector<Point> stack; // Empty FILO stack, filled with points to be coloured

	for (size_t firstPoint=0; firstPoint<allPoints.size(); firstPoint++){
		// Pick the first sausage not yet coloured
		if (allPoints[firstPoint].sausageID==1){ 
			if (params::verbosity>=VERBOSE) cout << "Starting sort of sausage #" << newSausage << endl;
			stack.push_back(allPoints[firstPoint]);
			while (stack.size()>0){
				Point curr = stack.back();
				stack.pop_back();
				// There's got to be a nice way of doing this
				curr.self->sausageID=newSausage;
				if (curr.left    && curr.left->sausageID==1)    stack.push_back(*(curr.left)) ;
				if (curr.right   && curr.right->sausageID==1)   stack.push_back(*(curr.right)) ;
				if (curr.up      && curr.up->sausageID==1)      stack.push_back(*(curr.up)) ;
				if (curr.down    && curr.down->sausageID==1)    stack.push_back(*(curr.down)) ;
				if (curr.forward && curr.forward->sausageID==1) stack.push_back(*(curr.forward)) ;
				if (curr.back    && curr.back->sausageID==1)    stack.push_back(*(curr.back)) ;
			}
			newSausage++;
		}
	}
}
