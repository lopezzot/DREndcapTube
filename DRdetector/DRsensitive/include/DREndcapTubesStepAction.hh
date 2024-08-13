//**************************************************************************
// \file DREndcapTubesStepAction.hh
// \brief:  definition of DREndcapTubesStepAction class
// \author: Lorenzo Pezzotti (CERN) @lopezzot
// \start date: 12 August 2024
//**************************************************************************

#ifndef DREndcapTubesStepAction_h
#define DREndcapTubesStepAction_h 1

// Includers from Geant4
#include "G4UserSteppingAction.hh"
#include "G4Step.hh"

// Includers from project files
#include "DREndcapTubesEvtAction.hh"

namespace dd4hep {
  namespace sim {
    
    class DREndcapTubesStepAction final : public G4UserSteppingAction {

	// Constructor and destructor
	//
	public:
	  DREndcapTubesStepAction(DREndcapTubesEvtAction* evtAction);
	  virtual ~DREndcapTubesStepAction() = default;
	
	// Virtual methods from base class
	//
	public:
	  virtual void UserSteppingAction(const G4Step* step) override;
	
	// Fields
	//
	private:
	  DREndcapTubesEvtAction* fEvtAction;
    };

  } // namespace sim
} // namespace dd4hep

#endif // DREndcapTubesStepAction_h

//**************************************************************************
