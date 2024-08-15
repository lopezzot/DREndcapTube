//**************************************************************************
// \file DREndcapTubesEvtAction.cpp
// \brief:  implementation of DREndcapTubesEvtAction class
// \author: Lorenzo Pezzotti (CERN) @lopezzot
// \start date: 12 August 2024
//**************************************************************************

// Includers from project files
#include "DREndcapTubesEvtAction.hh"

// Includers from Geant4
#include "G4AnalysisManager.hh"
#include "G4ParticleDefinition.hh"

namespace dd4hep {
  namespace sim {
       
    // Define constructor 
    //
    DREndcapTubesEvtAction::DREndcapTubesEvtAction()
      : G4UserEventAction(),
        EnergyScin(0.),
        SglScin(0) {}
    
    // Define BeginOfEventAction virtual method 
    //
    void DREndcapTubesEvtAction::BeginOfEventAction(const G4Event*) {
        
      EnergyScin = 0.;
      SglScin = 0.;
    }

    // Define EndOfEventAction virtual method
    // 
    void DREndcapTubesEvtAction::EndOfEventAction(const G4Event*) {

      // Fill ntuple of output file
      G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
      analysisManager->FillNtupleDColumn(0, EnergyScin);
      analysisManager->FillNtupleIColumn(1, SglScin);
      analysisManager->AddNtupleRow();
    }

  } // namespace sim
} // namespace dd4hep

//**************************************************************************
