/**
 * \file ERAlgorecoPi0.h
 *
 * \ingroup ERToolsPi0Reco
 * 
 * \brief Class def header for a class ERAlgorecoPi0
 *
 * @author rsjones
 */

/** \addtogroup ERToolsPi0Reco

    @{*/

#ifndef ERTOOL_ERALGORECOPI0_H
#define ERTOOL_ERALGORECOPI0_H

#include "ERTool/Base/AlgoBase.h"
#include "TNtuple.h"
namespace ertool {

  /**
     \class ERAlgorecoPi0
     User custom Analysis class made by kazuhiro
   */
  class ERAlgorecoPi0 : public AlgoBase {
  
  public:

    /// Default constructor
    ERAlgorecoPi0(const std::string& name="ERAlgorecoPi0");

    /// Default destructor
    virtual ~ERAlgorecoPi0(){}

    /// Reset function
    virtual void Reset();

    /// Function to accept fclite::PSet
    void AcceptPSet(const ::fcllite::PSet& cfg);

    /// Called @ before processing the first event sample
    void ProcessBegin();

    /// Function to evaluate input showers and determine a score
    bool Reconstruct(const EventData &data, ParticleGraph &graph);

    /// Called after processing the last event sample
    void ProcessEnd(TFile* fout=nullptr);
	
		// Structure to hold data on a shower
		struct showerData{
      TVector3 start, dir, p;
      double E;
      int ID, PDG;
			ertool::NodeID_t nodeID;
    };
		double Mcut, Dcut;
		double muM, sigmaM, AM, muD, sigmaD, AD;

 };
}
#endif

/** @} */ // end of doxygen group 
