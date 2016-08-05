/*
 * DDIRCMirHit.h
 *
 *  Created on: Oct 7, 2013
 *      Author: mpatsyuk
 */

#ifndef DDIRCMIRHIT_H_
#define DDIRCMIRHIT_H_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
using namespace jana;

class DDIRCMirHit: public JObject {

public:
	JOBJECT_PUBLIC (DDIRCMirHit);

	float x, y, z;	// true point of intersection
	//float px, py, pz; // true 3 momentum of particle
	float t;	// time
	//float E;	// energy
	int track;	///< Track number
	//int parentID; ///< Charged track which produced CKOV (plan to add with Geant4)
	int reflectionID;
	int volumeID;

	void toStrings(vector<pair<string, string> >&items) const {
		AddString(items, "x", "%1.3f", x);
		AddString(items, "y", "%1.3f", y);
		AddString(items, "z", "%1.3f", z);
		//AddString(items, "px", "%1.3f", px);
		//AddString(items, "py", "%1.3f", py);
		//AddString(items, "pz", "%1.3f", pz);
		AddString(items, "t", "%1.3f", t);
		//AddString(items, "E", "%1.3f", E);
		AddString(items, "track", "%d", track);
		//AddString(items, "parentID", "%d", parentID);
		AddString(items, "reflectionID", "%d", reflectionID);
		AddString(items, "volumeID", "%d" , volumeID);
	}
};

#endif /* DDIRCMIRHIT_H_ */
