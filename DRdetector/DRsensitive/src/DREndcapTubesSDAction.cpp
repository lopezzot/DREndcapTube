//**************************************************************************
// \file DREndcapTubesSDAction.cpp
// \brief: Implementation of Geant4SensitiveAction template class for
//         DREndcapTubes
// \author: Lorenzo Pezzotti (CERN) @lopezzot
// \start date: 11 August 2024
//**************************************************************************

// Includers from DD4HEP
//#include <DDG4/Geant4SensDetAction.inl>
//#include <DDG4/Geant4ParticleInformation.h>
//#include <DDG4/Factories.h>

#include "DD4hep/Version.h"
#include "DDG4/Geant4SensDetAction.inl"
#include "DDG4/Factories.h"
#include "DDG4/Geant4EventAction.h"
#include "DDG4/Geant4RunAction.h"
#include "DDG4/Geant4GeneratorAction.h"
#include "DDG4/Geant4Mapping.h"
#include "G4OpticalPhoton.hh"
#include "G4VProcess.hh"
#include "G4UserSteppingAction.hh"
#include "G4Types.hh"
#include "G4Poisson.hh"
#include "globals.hh"
#include "G4ProcessManager.hh"
#include "G4OpBoundaryProcess.hh"

// Includers from project files
#include "DREndcapTubesRunAction.hh"
#include "DREndcapTubesEvtAction.hh"
#include "DREndcapTubesStepAction.hh"

//#define DREndcapTubesSDDebug

namespace dd4hep {
  namespace sim {
    class DREndcapTubesSD {

    // Constructor and destructor
    //
    public:
      DREndcapTubesSD() = default;
      ~DREndcapTubesSD() = default;  

    // Methods
    //
    public:
      void clar() {
        // nothing to clear
      }

      // Method for generating hit(s)
      G4bool process(G4Step* /*const step*/, G4TouchableHistory*) {
        return true;
      }

      void beginRun(const G4Run* run){
        fRunAction = new DREndcapTubesRunAction;
	fEvtAction = new DREndcapTubesEvtAction;
	fStepAction = new DREndcapTubesStepAction(fEvtAction);
        fRunAction->BeginOfRunAction(run);
      }

      void endRun(const G4Run* run){
        fRunAction->EndOfRunAction(run);
      }

      void beginEvent(const G4Event* event){
        fEvtAction->BeginOfEventAction(event);
      }

      void endEvent(const G4Event* event){
        fEvtAction->EndOfEventAction(event);
      }

      void process(G4Step const * step, G4TouchableHistory* ) {
        fStepAction->UserSteppingAction(step);
      }

    // Fields
    //
    public:
      double eDep = 0.;
      DREndcapTubesRunAction* fRunAction;
      DREndcapTubesEvtAction* fEvtAction;
      DREndcapTubesStepAction* fStepAction;
      Geant4Sensitive*  sensitive{};
    };
  } // namespace sim
} // namespace dd4hep

namespace dd4hep {
  namespace sim   {

    // Function template specialization of Geant4SensitiveAction class.
    // Define actions
    template <> void Geant4SensitiveAction<DREndcapTubesSD>::initialize() {

      eventAction().callAtBegin(&m_userData, &DREndcapTubesSD::beginEvent);
      eventAction().callAtEnd(&m_userData, &DREndcapTubesSD::endEvent);

      runAction().callAtBegin(&m_userData,&DREndcapTubesSD::beginRun);
      runAction().callAtEnd(&m_userData,&DREndcapTubesSD::endRun);

      m_userData.sensitive = this;
    }

    // Function template specialization of Geant4SensitiveAction class.
    // Define collections created by this sensitivie action object
    template <> void Geant4SensitiveAction<DREndcapTubesSD>::defineCollections()    {
      //m_collectionID = declareReadoutFilteredCollection<MyTrackerSD::Hit>();
    }

    // Function template specialization of Geant4SensitiveAction class.
    // Method that accesses the G4Step object at each track step.
    template <> bool Geant4SensitiveAction<DREndcapTubesSD>::process(const G4Step* aStep, G4TouchableHistory* history ) {
    
      #ifdef DREndcapTubesSDDebug
      //Print out some info step-by-step in sensitive volumes
      //
      std::cout<<"Track #: "<< aStep->GetTrack()->GetTrackID()<< " " <<
                 "Step #: " << aStep->GetTrack()->GetCurrentStepNumber()<< " "<<
                 "Volume: " << aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName()<< " " << std::endl;
      std::cout<<"x: "<< aStep->GetPreStepPoint()->GetPosition().x() <<
                 "y: "<< aStep->GetPreStepPoint()->GetPosition().y() <<
                 "z: "<< aStep->GetPreStepPoint()->GetPosition().z() << std::endl;
      std::cout<<"Particle "<< aStep->GetTrack()->GetParticleDefinition()->GetParticleName()<< " " <<
                 "Dep(MeV) "<< aStep->GetTotalEnergyDeposit() << " " <<
                 "Mat "     << aStep->GetPreStepPoint()->GetMaterial()->GetName() << " " << 
                 "Vol "     << aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName() << std::endl; 
      #endif
      m_userData.process(aStep, history);
      /*Geant4StepHandler h(step);
      Position  prePos    = h.prePos();
      Position  postPos   = h.postPos();
      Position  direction = postPos - prePos;
      Position  pos       = mean_direction(prePos,postPos);
      Direction mom       = 0.5 * (h.preMom() + h.postMom());
      double    hit_len   = direction.R();
      double    depo      = h.deposit();
      double    tim       = h.track->GetGlobalTime();
      // Somehow extract here the physics you want
      MyTrackerSD::Hit* hit = 
	new MyTrackerSD::Hit(h.trkID(), h.trkPdgID(), depo, tim, hit_len, pos, mom);
      Geant4HitData::MonteCarloContrib contrib = Geant4HitData::extractContribution(step);
      hit->cellID        = cellID(step);
      hit->step_length   = hit_len;
      hit->prePos        = prePos;
      hit->postPos       = postPos;
      collection(m_collectionID)->add(hit);
      mark(h.track);
      if ( 0 == hit->cellID )  {
        hit->cellID = volumeID(step);
        except("+++ Invalid CELL ID for hit!");
      }
      printP1("Hit with deposit:%f  Pos:%f %f %f ID=%016X",
	      depo, pos.X(), pos.Y(), pos.Z(), (void*)hit->cellID);
      Geant4TouchableHandler handler(step);
      print("    Geant4 path:%s", handler.path().c_str());

      // Do something with my personal data (can be also something more clever ;-):
      m_userData.integratedDeposit += contrib.deposit;
      ++m_userData.mumDeposits;

      /// Let's play with the Geant4TrackInformation
      /// See issue https://github.com/AIDASoft/DD4hep/issues/1073
      if ( nullptr == h.track->GetUserInformation() )   {
	auto data = std::make_unique<ParticleUserData>();
	data->absolute_momentum = h.track->GetMomentum().mag();
	h.track->SetUserInformation(new Geant4ParticleInformation(std::move(data)));
      }*/
      return true;
    } // end of Geant4SensitiveAction::process() method specialization

  } // namespace sim
} // namespace dd4hep

//--- Factory declaration
namespace dd4hep { namespace sim {
    typedef Geant4SensitiveAction<DREndcapTubesSD> DREndcapTubesSDAction;
  }}
DECLARE_GEANT4SENSITIVE(DREndcapTubesSDAction)

//**************************************************************************
