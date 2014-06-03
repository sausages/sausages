#include <iostream>
#include <vector>
#include <stdexcept>
#include "point.h"
#include "find_sausages.h"

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
	cout << "begin flood-fill" << endl;
	int newSausage=2; // Start at 2 (0 is no-sausage and 1 is unsorted)
	cout << "newSausage is " << newSausage << endl;
	vector<Point> stack; // Empty FILO stack, filled with points to be coloured

	for (size_t firstPoint=0; firstPoint<allPoints.size(); firstPoint++){
		// Pick the first sausage not yet coloured
		if (allPoints[firstPoint].sausageID==1){ 
			cout << "Starting sort of sausage #" << newSausage << endl;
			stack.push_back(allPoints[firstPoint]);
			while (stack.size()>0){
				cout << "stack size is " << stack.size() << endl;
				Point curr = stack.back();
				stack.pop_back();
				cout << "post-pop stack size is " << stack.size() << endl;
				//curr.sausageID=newSausage;
				curr.self->sausageID=newSausage;
				// There's got to be a nice way of doing this
				if (curr.left->sausageID==1)    {stack.push_back(*(curr.left)) ; cout << "left" << endl;}
				if (curr.right->sausageID==1)   {stack.push_back(*(curr.right)) ; cout << "right" << endl;}
				if (curr.up->sausageID==1)      {stack.push_back(*(curr.up)) ; cout << "up" << endl;}
				if (curr.down->sausageID==1)    {stack.push_back(*(curr.down)) ; cout << "down" << endl;}
				if (curr.forward->sausageID==1) {stack.push_back(*(curr.forward)) ; cout << "forward" << endl;}
				if (curr.back->sausageID==1)    {stack.push_back(*(curr.back)) ; cout << "back" << endl;}
				cout << "post-pushes stack size is " << stack.size() << endl;
			}
			newSausage++;
		}
	}
}
