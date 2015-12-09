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

	void ERAlgofindVertex::ProcessBegin(){
		nt = new TNtuple("vertex","vertex","evIndex:vertexID:X:Y:Z:nodeID:partX:partY:partZ:trackshower");
	}

	bool ERAlgofindVertex::Reconstruct(const EventData &data, ParticleGraph& graph){
		// Here is how this code works:
		// - First we extract all the relevent data from the EventData and the ParticleGraph and add it to a vector
		//   of particles called partList.
		// - We choose the first particle in the event and add a "vertex" at it's start.
		//   This is done by pushing a vertex object on to a a vector of vertex objects called vertexList
		//   A vertex object holds the following:
		//   	pos      : The position vector of the vertex
		//   	outList  : A list of the particles which exit the vertex
		//		R        : Radius of the smallest sphere which will hold all of the particles in the vertex
		// - Then we then move to the next particle and add it to the outList of an existing vertex in
		//   the vertexList if its start position is within a distance tolerance (distTol) of the vertex, 
		//   else we add a new vertex to the vertexList.
		// - Iterate...
		// - Once we have looked at all of the particles in the partList we calcualte the R values for each
		//   vertex (this allows us to see the suitability of our distTol choice)
		//   

		evIndex++;
		std::cout << "Event : " << evIndex << " ---------------------------" << std::endl;

		// Extract the data and push it onto a vector of part structures
		std::vector<part> partList;
		for(auto const& graphPart : graph.GetParticleArray()){
			part thisPart;
			thisPart.nodeID = graphPart.ID();
			thisPart.recoType = graphPart.RecoType();	
			if (thisPart.recoType == ertool::kTrack){
				auto const& track = data.Track(graphPart.RecoID());
				thisPart.start = TVector3(track[0][0], track[0][1], track[0][2]);
			}
			if (thisPart.recoType == ertool::kShower){
				auto const& shower = data.Shower(graphPart.RecoID());

				thisPart.start = TVector3(shower.Start()[0], shower.Start()[1], shower.Start()[2]);
			}
			partList.push_back(thisPart);
		}

		// Define the distance tolerance which we use as the radius of a sphere which holds
		// all the points in a vertex
		double distTol = 4;

		std::vector<vertex> vertexList;	
		for (part const &p : partList){
			bool startUsed = false;
			for (vertex &v : vertexList){
				// See if the start position matches with any existing vertex
				if ( (v.pos - p.start).Mag() < distTol ){
					if (startUsed){
						std::cout << "Warning: Particle " << p.nodeID << "'s start point has already been assigned to a vertex." << std::endl;
					}
					else{
						v.outList.push_back(p);
						startUsed = true;
							// Re-calculate the vertex position as the average of all of the points in the position
						TVector3 newPos(0,0,0);
						for (part const &outPart : v.outList){
							newPos = newPos + outPart.start;
						}
						v.pos = (1/v.outList.size())*newPos;
					}
				}
			}
			// If the particle does not match with any existing vertex then make a new vertex
			if (!startUsed){
				vertex newVertex;
				newVertex.pos = p.start;
				newVertex.outList.push_back(p);
				vertexList.push_back(newVertex);
			}
		}
		// Calculate the smallest radius, R which holds all the points in each vertex
		for (vertex &v : vertexList){
			v.R = 2*distTol;
			for (part const &p : v.outList){
				if ((v.pos - p.start).Mag() < v.R){
					v.R = (v.pos - p.start).Mag();
				}
			}	
		}

		// vertexID:X:Y:Z:nodeID:partX:partY:partZ:inout:trackshower
		// Let's see some stuff
		int i = 0;
		for (vertex &v : vertexList){
			std::cout << "Vertex : " << i << std::endl;
			std::cout << "     Out list : " << std::endl;
			for (part const &p : v.outList){
				std::cout << "          Particle : " << p.nodeID << std::endl;
				if (p.recoType == ertool::kTrack){
					nt->Fill(evIndex,i,v.pos.X(),v.pos.Y(),v.pos.Z(),p.nodeID,p.start.X(),p.start.Y(),p.start.Z(),0);
				}
				if (p.recoType == ertool::kShower){
					nt->Fill(evIndex,i,v.pos.X(),v.pos.Y(),v.pos.Z(),p.nodeID,p.start.X(),p.start.Y(),p.start.Z(),1);
				}
			}
			std::cout << "     R        : " << v.R << std::endl; 
			
			i++; 
		}
		if (evIndex == 50){std::cin.get();}
		return true;
	}

	void ERAlgofindVertex::ProcessEnd(TFile* fout){
		nt->Write();
	}

}

#endif
