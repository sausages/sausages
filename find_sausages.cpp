#include <vector>
#include <stdexcept>
#include "main.h"
#include "find_sausages.h"

using namespace std;

void threshold(vector<Point> &allPoints, double threshold_level){
	vector<Point>::iterator it;
	for (it=allPoints.begin(); it!=allPoints.end(); ++it){
		if (it->cl < threshold_level){
			it->sausageID=1; // In a sausage
		} else {
			it->sausageID=0; // Not in a sausage
		}
	}
}

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

