/*
 * DDIRCTruthPoint.h
 *
 *  Created on: Oct 7, 2013
 *      Author: yqiang
 */

#ifndef DDIRCTRUTHPOINT_H_
#define DDIRCTRUTHPOINT_H_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
using namespace jana;

class DDIRCTruthPoint: public JObject {

public:
	JOBJECT_PUBLIC (DDIRCTruthPoint);

	float x, y, z;	// true point of intersection
	float px, py, pz; // true 3 momentum of particle
	float t;	// time
	float E;	// energy
	int track;	///< Track number
	int itrack;	///< MCThrown track index
	int primary;	///< primary track=1    not primary track=0
	int ptype;    /// particle type
	int barnum; // number of the radiator hit
	int bboxnum; // number of the bar box hit

	void toStrings(vector<pair<string, string> >&items) const {
		AddString(items, "x", "%1.3f", x);
		AddString(items, "y", "%1.3f", y);
		AddString(items, "z", "%1.3f", z);
		AddString(items, "px", "%1.3f", px);
		AddString(items, "py", "%1.3f", py);
		AddString(items, "pz", "%1.3f", pz);
		AddString(items, "t", "%1.3f", t);
		AddString(items, "E", "%1.3f", E);
		AddString(items, "track", "%d", track);
		AddString(items, "itrack", "%d", itrack);
		AddString(items, "primary", "%d", primary);
		AddString(items, "ptype", "%d", ptype);
		//AddString(items, "barnum", "%d", barnum);
		//AddString(items, "bboxnum", "%d", bboxnum);
	}
};

#endif /* DDIRCTRUTHPOINT_H_ */
