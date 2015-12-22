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
			else{
				if (thisPart.recoType == ertool::kShower){
					auto const& shower = data.Shower(graphPart.RecoID());
					thisPart.start = TVector3(shower.Start()[0], shower.Start()[1], shower.Start()[2]);
				}
				else{
					// This is something that has been added by the pi0 reco
					// It is an extra type of particle to analyse
					thisPart.start = TVector3(graphPart.Vertex()[0], graphPart.Vertex()[1], graphPart.Vertex()[2]);
				}
			}
			partList.push_back(thisPart);
		}

		// Define the distance tolerance which we use as the radius of a sphere which holds
		// all the points in a vertex
		double distTol = 0.05;
			
		std::vector<vertex> vertexList;
		for (part &p : partList){
			bool startUsed = false;
			int vi = 0;
			for (vertex &v : vertexList){
				// See if the start position matches with any existing vertex
				if ( (v.pos - p.start).Mag() < distTol ){
					if (startUsed){
						// See if this new vertex is better than the current on
						if ( (v.pos - p.start).Mag() < (vertexList[p.vtx].pos - p.start).Mag() ){
							// Add the particle to the new vertex
							v.outList.push_back(p);

							// Re-calculate the new vertex position as the average of all of the points in the vertex
							TVector3 newPos(0,0,0);
							for (part const &outPart : v.outList){
								newPos = newPos + outPart.start;
							}
							v.pos = (1/double(v.outList.size()))*newPos;

							// Remove the particle from the old vertex
							int i = 0;
							for (part const &outPart : vertexList[p.vtx].outList){
								if (outPart.nodeID == p.nodeID){
									vertexList[p.vtx].outList.erase(vertexList[p.vtx].outList.begin()+i);
									break;
								}
								i++;
							}

							// Re-calculate the old vertex position as the average of all of the points in the vertex
							newPos = TVector3(0,0,0);
							for (part const &outPart : vertexList[p.vtx].outList){
								newPos = newPos + outPart.start;
							}
							vertexList[p.vtx].pos = (1/double(vertexList[p.vtx].outList.size()))*newPos;

							// Set the new vertex as THE vertex.
							p.vtx = vi;
						}
					}
					else{
						v.outList.push_back(p);
						p.vtx = vi;
						startUsed = true;
						// Re-calculate the vertex position as the average of all of the points in the vertex
						TVector3 newPos(0,0,0);
						for (part const &outPart : v.outList){
							newPos = newPos + outPart.start;
						}
						v.pos = (1/double(v.outList.size()))*newPos;
					}
				}
				vi++;
			}
			// If the particle does not match with any existing vertex then make a new vertex
			if (!startUsed){
				vertex newVertex;
				newVertex.pos = p.start;
				newVertex.outList.push_back(p);
				vertexList.push_back(newVertex);
				p.vtx = vertexList.size()-1;
			}
		}

		// Set the vertices in the particle graph
		for (auto const &v : vertexList){
			for (auto const &p : v.outList){
				auto &graphPart  = graph.GetParticle(p.nodeID);
				auto const &PDG  = graphPart.PdgCode();

				// Calculate the momentum and mass of tracks
				if (p.recoType == ertool::kTrack){
					auto const& track = data.Track(graphPart.RecoID());
					if (track._pid != 1 && track._pid != 2 && track._pid != 3 && track._pid != 4){
						// We don't know what this is!
						// Use the defaults
						auto const &mass = graphPart.Mass();
						auto const &mom  = graphPart.Momentum();
						graphPart.SetParticleInfo(PDG, mass, v.pos, mom);
					}
					else{
						double mass;
						switch (track._pid){
							case 1:
								// Proton
								mass = 938.272;
								break;
							case 2:
								// Kaon
								mass = 493.677;
								break;
							case 3:
								// Pion
								mass = 139.570;
								break;
							case 4:
								// Muon
								mass = 105.658;
								break;
						}
						TVector3 end(track[track.size()-1][0], track[track.size()-1][1], track[track.size()-1][2]);
						TVector3 dir = (1/(end-p.start).Mag())*(end-p.start);
						TVector3 mom = std::pow( std::pow(track._energy, 2) + 2*mass*track._energy, 0.5)*dir;
						graphPart.SetParticleInfo(PDG, mass, v.pos, mom);
					}
				}
				else{
					auto const &mass = graphPart.Mass();
					auto const &mom  = graphPart.Momentum();
					graphPart.SetParticleInfo(PDG, mass, v.pos, mom);
				}
			}
		}
		

		return true;
	}

	void ERAlgofindVertex::ProcessEnd(TFile* fout){
	}

}

#endif
