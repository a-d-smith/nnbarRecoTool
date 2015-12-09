/**
 * \file ERAlgofindVertex.h
 *
 * \ingroup findVertex
 * 
 * \brief Class def header for a class ERAlgofindVertex
 *
 * @author rsjones
 */

/** \addtogroup findVertex

    @{*/

#ifndef ERTOOL_ERALGOFINDVERTEX_H
#define ERTOOL_ERALGOFINDVERTEX_H

#include "ERTool/Base/AlgoBase.h"
#include "TNtuple.h"

namespace ertool {

  /**
     \class ERAlgofindVertex
     User custom Algorithm class made by kazuhiro
   */
  class ERAlgofindVertex : public AlgoBase {
  
  public:

    /// Default constructor
    ERAlgofindVertex(const std::string& name="ERAlgofindVertex");

    /// Default destructor
    virtual ~ERAlgofindVertex(){};

    /// Reset function
    void Reset();

    /// Function to accept fclite::PSet
    void AcceptPSet(const ::fcllite::PSet& cfg);

    /// Called @ before processing the first event sample
    void ProcessBegin();

    /// Function to evaluate input showers and determine a score
    bool Reconstruct(const EventData &data, ParticleGraph& graph);

    /// Called after processing the last event sample
    void ProcessEnd(TFile* fout=nullptr);

		struct part{
			ertool::NodeID_t nodeID;
			ertool::RecoType_t recoType;
			TVector3 start;
		};

		struct vertex{
			TVector3 pos;
			std::vector<part> outList;
			double R;
		};
		
		int evIndex = 0;
		TNtuple *nt;
  };
}
#endif

/** @} */ // end of doxygen group 
