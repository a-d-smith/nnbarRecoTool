#ifndef ERTOOL_ERALGOFINDVERTEX_CXX
#define ERTOOL_ERALGOFINDVERTEX_CXX

#include "ERAlgofindVertex.h"

namespace ertool {

	ERAlgofindVertex::ERAlgofindVertex(const std::string& name) : AlgoBase(name)
	{}

	void ERAlgofindVertex::Reset()
	{}

	void ERAlgofindVertex::AcceptPSet(const ::fcllite::PSet& cfg)
	{}

	void ERAlgofindVertex::ProcessBegin()
	{}

	bool ERAlgofindVertex::Reconstruct(const EventData &data, ParticleGraph& graph){
		// Here is how this code works:
		//

		// Extract the data and push it onto a vector of part structures
		std::vector<part> partList;
		for(auto const& graphPart : graph.GetParticleArray()){
			part thisPart;
			thisPart.nodeID = graphPart.ID();
			thisPart.recoType = graphPart.RecoType();	
			if (thisPart.recoType == ertool::kTrack){
				auto const& track = data.Track(graphPart.RecoID());
				thisPart.start = TVector3(track[0][0], track[0][1], track[0][2]);
				thisPart.end   = TVector3(track[track.size()-1][0], track[track.size()-1][1], track[track.size()-1][2]);
			}
			if (thisPart.recoType == ertool::kShower){
				auto const& shower = data.Shower(graphPart.RecoID());

				thisPart.start = TVector3(shower.Start()[0], shower.Start()[1], shower.Start()[2]);
      	TVector3 thisp(graphPart.Momentum()[0], graphPart.Momentum()[1], graphPart.Momentum()[2]);
	      TVector3 thisDir = (1/thisp.Mag())*thisp;
	      double thisL     = shower.Length();
				thisPart.end = thisPart.start + thisL*thisDir;
			}
			partList.push_back(thisPart);
		}

		// Define the distance tolerance which we use as the radius of a sphere which holds
		// all the points in a vertex
		double distTol = 1;

		
		
		return true;
	}

	void ERAlgofindVertex::ProcessEnd(TFile* fout)
	{}

}

#endif
