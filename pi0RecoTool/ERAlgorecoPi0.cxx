#ifndef ERTOOL_ERALGORECOPI0_CXX
#define ERTOOL_ERALGORECOPI0_CXX

#include "ERAlgorecoPi0.h"

namespace ertool {

	ERAlgorecoPi0::ERAlgorecoPi0(const std::string& name) : AlgoBase(name)
	{}

	void ERAlgorecoPi0::Reset()
	{}

	void ERAlgorecoPi0::AcceptPSet(const ::fcllite::PSet& cfg)
	{}

	void ERAlgorecoPi0::ProcessBegin(){
		// Setup an NTuple that holds: cut, Nc, Nm, Ni, Nn
		// M and D cuts chosen from the efficiency plot
		Mcut = 0.001;
		Dcut = 0.01;

		// PDF data
		muM    = 128.134;
		sigmaM = 25.5467;
		AM     = 1/64.0361;
		muD    = -0.174207;
		sigmaD = 1.50479;
		AD     = 1/(sigmaD*std::pow(2*3.14159265358979323,0.5));

	}

	bool ERAlgorecoPi0::Reconstruct(const EventData &data, ParticleGraph &graph){
		std::vector<showerData> shower;

		// Loop over all showers in the event and add them to the shower vector
		for (auto const& sh : graph.GetParticleNodes(RecoType_t::kShower)){
			auto const& sh1 = data.Shower(graph.GetParticle(sh).RecoID());
			auto const& part = graph.GetParticle(sh);

			showerData thisShower;

			TVector3 thisStart(sh1.Start()[0], sh1.Start()[1], sh1.Start()[2]);
			TVector3 thisp(part.Momentum()[0], part.Momentum()[1], part.Momentum()[2]);
			TVector3 thisDir = (1/thisp.Mag())*thisp;

			double           thisE      = part.Energy();
			int              thisID     = part.ID();
			ertool::NodeID_t thisNodeID = sh;
			int    thisPdg   = part.PdgCode();

			thisShower.start  = thisStart;
			thisShower.dir    = thisDir;
			thisShower.p      = thisp;
			thisShower.E      = thisE;
			thisShower.ID     = thisID;
			thisShower.nodeID = thisNodeID;
			thisShower.PDG    = thisPdg;

			shower.push_back(thisShower);
		}
		// Loop over a range of cut values {
		// Loop over all unique photon shower pairs 
		for (showerData const& s1 : shower){
			for (showerData const& s2 : shower){
				if (s1.ID > s2.ID){
					if (s1.PDG == 22 && s2.PDG == 22){
						if (!isnan(s1.E) && !isnan(s2. E)){ 
							// Calculate some useful numbers for the geometry of the problem
							double costheta  = s1.dir.Dot(s2.dir);
							double sin2theta = 1-std::pow(costheta,2);
							TVector3 deltar = s1.start - s2.start;
							double n1 = deltar.Dot(s1.dir);
							double n2 = deltar.Dot(s2.dir);
							double lambda1 = (costheta*n2-n1)/sin2theta;
							double lambda2 = (n2-costheta*n1)/sin2theta;

							// Calcualte their distance of closest approach, D
							double D = (s1.start - s2.start + lambda1*s1.dir - lambda2*s2.dir).Mag();
							
							// Calculate their invariant mass, M
							double M = std::pow(2 * s1.E * s2.E * (1-costheta), 0.5);
		
							// Calculate the PDF for M
							double PM     = AM*std::exp(-0.5*std::pow((M-muM)/sigmaM,2));
							// Calculate the PDF for D
							double PD     = AD*std::exp(-0.5*std::pow((std::log(D)-muD)/sigmaD,2));
	
							// We choose the pair as having a pi0 parent if:
							// 		PM >= Mcut
							// 		PD >= Dcut
							if (PM >= Mcut && PD >= Dcut){
								// Find the point of closest approach
								// Treat this as the pi0 vertex
								TVector3 P = 0.5*(s1.start + s2.start + lambda1*s1.dir + lambda2*s2.dir);

								// Find the momentum 
								auto &pi0 = graph.CreateParticle();
								pi0.SetParticleInfo(111, 134.977, P, s1.p + s2.p);
							}
						}
						else {
							std::cout << "Warning: Shower energy is not a number" << std::endl;
						}
					}
				}
			}
		}
		return true;
	}

	void ERAlgorecoPi0::ProcessEnd(TFile* fout){
	}

}

#endif
