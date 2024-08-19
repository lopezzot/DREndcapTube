//**************************************************************************
// \file DREndcapTubesEvtAction.hh
// \brief:  definition of DREndcapTubesEvtAction class
// \author: Lorenzo Pezzotti (CERN) @lopezzot
// \start date: 12 August 2024
//**************************************************************************

#ifndef DREndcapTubesEvtAction_h
#define DREndcapTubesEvtAction_h 1

// Includers from Geant4
#include "G4UserEventAction.hh"
#include "globals.hh"

namespace dd4hep {
  namespace sim {
    class DREndcapTubesEvtAction final : public G4UserEventAction {
      
      public:
	// Constructor and destructor
	//
        DREndcapTubesEvtAction();
        virtual ~DREndcapTubesEvtAction() = default;

	// Virtual methods from base class
	//
        virtual void BeginOfEventAction(const G4Event*) override;
        virtual void EndOfEventAction(const G4Event*) override;
 
        // Inline methods
	inline void AddEdepScin(G4double edep) { EnergyScin += edep; }
	inline void AddSglScin(G4int scinsgl) { SglScin += scinsgl; }
	inline void AddEdepCher(G4double edep) { EnergyCher += edep; }
	inline void AddSglCher(G4int chersgl) { SglCher += chersgl; }

      private:
	// Fields
	//
	G4double EnergyScin;
	G4int SglScin;
	G4double EnergyCher;
	G4int SglCher;
    };
  } // namespace sim
} // namespace dd4hep

#endif // DREndcapTubesEvtAction_h

//**************************************************************************
