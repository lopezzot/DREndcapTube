//**************************************************************************
// \file DREndcapTubesRunAction.hh
// \brief:  definition of DREndcapTubesRunAction class
// \author: Lorenzo Pezzotti (CERN) @lopezzot
// \start date: 12 August 2024
//**************************************************************************

#ifndef DREndcapTubesRunAction_h
#define DREndcapTubesRunAction_h 1


//Includers from Geant4
//
#include "G4UserRunAction.hh"
#include "globals.hh"
#include "G4Timer.hh"

// Includers from project files
//#include "DREndcapTubesEventAction.h"

// Forward declarations from Geant4
class G4Run;

namespace dd4hep {
  namespace sim {
    class DREndcapTubesRunAction final : public G4UserRunAction {
            
      public:
        // Constructor and destructor
        //
        DREndcapTubesRunAction(/*DREndcapTubesEventAction* eventAction*/);
        virtual ~DREndcapTubesRunAction() = default;

        // Virtual methods from base class
        //
        virtual void BeginOfRunAction(const G4Run*) override;
        virtual void EndOfRunAction(const G4Run*) override;

	// Fields
        private:
          //DREndcapTubesEventAction* fEventAction;
	  G4Timer fTimer;

    };
  }
}

#endif // DREndcapTubesRunAction_h

//**************************************************************************
