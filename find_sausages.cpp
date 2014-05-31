#include <vector>
#include "main.h"
#include "find_sausages.h"

using namespace std;

void threshold(vector<Point> &allPoints, double threshold_level){
	vector<Point>::iterator it;
	for (it=allPoints.begin(); it!=allPoints.end(); ++it){
		if (it->cl > threshold_level){
			it->sausageID=1; // In a sausage
		} else {
			it->sausageID=0; // Not in a sausage
		}
	}


}
