//**************************************************************************
// \file DREndcapTubesStepAction.cpp
// \brief:  implementation of DREndcapTubesStepAction class
// \author: Lorenzo Pezzotti (CERN) @lopezzot
// \start date: 12 August 2024
//**************************************************************************

// Includers from project files
#include "DREndcapTubesStepAction.hh"

namespace dd4hep {
  namespace sim {

      // Constructor
      //
      DREndcapTubesStepAction::DREndcapTubesStepAction(DREndcapTubesEvtAction* evtAction)
	: G4UserSteppingAction(),
	  fEvtAction(evtAction) {}

      // UserSteppingAction method
      //
      void DREndcapTubesStepAction::UserSteppingAction(const G4Step* step) {
        fEvtAction->AddEdepScin(step->GetTotalEnergyDeposit()); 
      }
 
  } // namespace dd4hep
} // namespace sim

//**************************************************************************
