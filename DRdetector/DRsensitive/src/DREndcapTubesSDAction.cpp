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
#include "DD4hep/Segmentations.h"

// Includers from Geant4
#include "G4OpBoundaryProcess.hh"

// Includers from project files
#include "DREndcapTubesRunAction.hh"
#include "DREndcapTubesEvtAction.hh"
#include "DREndcapTubesSglHpr.hh"

//#define DREndcapTubesSDDebug

namespace dd4hep {
  namespace sim {
    class DREndcapTubesSDData {

    // Constructor and destructor
    //
    public:
      DREndcapTubesSDData() = default;
      ~DREndcapTubesSDData() = default;  

    // Methods
    //
    public:

      void beginRun(const G4Run* run){
        fRunAction = new DREndcapTubesRunAction;
	fEvtAction = new DREndcapTubesEvtAction;
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

    // Fields
    //
    public:
      DREndcapTubesRunAction* fRunAction;
      DREndcapTubesEvtAction* fEvtAction;
      Geant4Sensitive*  sensitive{};
      int collection_cher;
    };
  } // namespace sim
} // namespace dd4hep

namespace dd4hep {
  namespace sim   {

    // Function template specialization of Geant4SensitiveAction class.
    // Define actions
    template <> void Geant4SensitiveAction<DREndcapTubesSDData>::initialize() {

      eventAction().callAtBegin(&m_userData, &DREndcapTubesSDData::beginEvent);
      eventAction().callAtEnd(&m_userData, &DREndcapTubesSDData::endEvent);

      runAction().callAtBegin(&m_userData,&DREndcapTubesSDData::beginRun);
      runAction().callAtEnd(&m_userData,&DREndcapTubesSDData::endRun);

      m_userData.sensitive = this;
    }

    // Function template specialization of Geant4SensitiveAction class.
    // Define collections created by this sensitivie action object
    template <> void Geant4SensitiveAction<DREndcapTubesSDData>::defineCollections()    {
      std::string ROname = m_sensitive.readout().name();
      m_collectionID = defineCollection<Geant4Calorimeter::Hit>(ROname+"Scin");
      m_userData.collection_cher = defineCollection<Geant4Calorimeter::Hit>(ROname+"Cher");
    }

    // Function template specialization of Geant4SensitiveAction class.
    // Method that accesses the G4Step object at each track step.
    template <> bool Geant4SensitiveAction<DREndcapTubesSDData>::process(const G4Step* aStep, G4TouchableHistory* /*history*/ ) {
    
      G4double Edep = aStep->GetTotalEnergyDeposit();

      dd4hep::BitFieldCoder decoder("tank:1,endcap:1,stave:10,tower:8,air:1,col:10,row:7,clad:1,core:1,cherenkov:1");
      auto VolID = volumeID(aStep);
      auto CherenkovID = decoder.get(VolID,"cherenkov");
      // Skip this step if edep is 0 and it is a scintillating fiber
      if(CherenkovID==0 && Edep==0.) return true;
      [[maybe_unused]] auto EndcapID = decoder.get(VolID,"endcap");
      [[maybe_unused]] auto StaveID = decoder.get(VolID,"stave");
      [[maybe_unused]] auto TowerID = decoder.get(VolID,"tower");
      [[maybe_unused]] auto AirID = decoder.get(VolID,"air");
      [[maybe_unused]] auto ColID = decoder.get(VolID,"col");
      [[maybe_unused]] auto RowID = decoder.get(VolID,"row");
      [[maybe_unused]] auto CladID = decoder.get(VolID,"clad");
      [[maybe_unused]] auto CoreID = decoder.get(VolID,"core");

      G4bool IsCherenkov = CherenkovID; // 1 for cherenkov 0 for scintillating fibers

      // If it is an optical photon inside the CLADDING of Cherenkov fibers kill it
      if (CladID==1 && CoreID==0 && CherenkovID==1
	  && aStep->GetTrack()->GetParticleDefinition() == G4OpticalPhoton::Definition() ) {
        aStep->GetTrack()->SetTrackStatus( fStopAndKill );
      }

      #ifdef DREndcapTubesSDDebug
      //Print out some info step-by-step in sensitive volumes
      //
      std::cout<<"-------------------------------"<<std::endl;
      std::cout<<"--> DREndcapTubes: track info: "<<std::endl;
      std::cout<<"----> Track #: "<< aStep->GetTrack()->GetTrackID()<< " " <<
                 "Step #: " << aStep->GetTrack()->GetCurrentStepNumber()<< " "<<
                 "Volume: " << aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName()<< " " << std::endl;
      std::cout<<"--> DREndcapTubes:: position info: "<<std::endl;
      std::cout<<"----> x: "<< aStep->GetPreStepPoint()->GetPosition().x() <<
                 "y: "<< aStep->GetPreStepPoint()->GetPosition().y() <<
                 "z: "<< aStep->GetPreStepPoint()->GetPosition().z() << std::endl;
      std::cout<<"--> DREndcapTubes: particle info: "<<std::endl;
      std::cout<<"----> Particle "<< aStep->GetTrack()->GetParticleDefinition()->GetParticleName()<< " " <<
                 "Dep(MeV) "<< aStep->GetTotalEnergyDeposit() << " " <<
                 "Mat "     << aStep->GetPreStepPoint()->GetMaterial()->GetName() << " " << 
                 "Vol "     << aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName() << std::endl; 
      std::cout<<"--> DREndcapTubes: physVolID info: "<<std::endl;
      std::cout<<"----> Volume ID "<<volumeID(aStep)<<" Endcap ID "<<EndcapID<<" Stave ID "<<StaveID<<std::endl;
      std::cout<<"----> Tower ID "<<TowerID<<" Air ID "<<AirID<<std::endl;
      std::cout<<"----> Col ID "<<ColID<<" Row ID "<<RowID<<std::endl;
      std::cout<<"----> Clad ID "<<CladID<<" Core ID "<<CoreID<<" Cherenkov ID "<<CherenkovID<<std::endl;
      #endif 

      // We now calculate the signal in S and C fiber according to the step contribution
      //
      G4double steplength = aStep->GetStepLength();
      G4int signalhit = 0;
      Geant4HitCollection* coll = (IsCherenkov ? collection(m_userData.collection_cher) : collection(m_collectionID));

      if(!IsCherenkov){ // it is a scintillating fiber
 
	m_userData.fEvtAction->AddEdepScin(Edep);
        if ( aStep->GetTrack()->GetParticleDefinition() == G4OpticalPhoton::Definition() ) {
          aStep->GetTrack()->SetTrackStatus( fStopAndKill );
        }

        if ( aStep->GetTrack()->GetDefinition()->GetPDGCharge() == 0 || steplength == 0. ) {
	    return true; // not ionizing particle
	}
        G4double distance_to_sipm = DREndcapTubesSglHpr::GetDistanceToSiPM(aStep);
        signalhit = DREndcapTubesSglHpr::SmearSSignal( DREndcapTubesSglHpr::ApplyBirks( Edep, steplength ) );
        signalhit = DREndcapTubesSglHpr::AttenuateSSignal(signalhit, distance_to_sipm);
	if(signalhit == 0) return true;
	m_userData.fEvtAction->AddSglScin(signalhit);
      } // end of scintillating fibre sigal calculation

      else{ // it is a Cherenkov fiber
	// save mc truth info in analysismanager auxiliary outputfile
	m_userData.fEvtAction->AddEdepCher(Edep);
        // calculate the signal in terms of Cherenkov photo-electrons
        if ( aStep->GetTrack()->GetParticleDefinition() == G4OpticalPhoton::Definition() ) {
          
	  G4OpBoundaryProcessStatus theStatus = Undefined;
          G4ProcessManager* OpManager = G4OpticalPhoton::OpticalPhoton()->GetProcessManager();

          if (OpManager) {
            G4int MAXofPostStepLoops = OpManager->GetPostStepProcessVector()->entries();
            G4ProcessVector* fPostStepDoItVector = OpManager->GetPostStepProcessVector(typeDoIt);

            for ( G4int i=0; i<MAXofPostStepLoops; i++) {
              G4VProcess* fCurrentProcess = (*fPostStepDoItVector)[i];
              G4OpBoundaryProcess* fOpProcess = dynamic_cast<G4OpBoundaryProcess*>(fCurrentProcess);
              if (fOpProcess) { theStatus = fOpProcess->GetStatus(); break; }
            }
          }

          switch ( theStatus ){
                                 
            case TotalInternalReflection: {
              
              G4double distance_to_sipm = DREndcapTubesSglHpr::GetDistanceToSiPM(aStep);
              G4int c_signal = DREndcapTubesSglHpr::SmearCSignal( );
              signalhit = DREndcapTubesSglHpr::AttenuateCSignal(c_signal, distance_to_sipm);
	      if(signalhit == 0) return true;
	      // save mc truth info in analysismanager auxiliary outputfile
              m_userData.fEvtAction->AddSglCher(signalhit);
              aStep->GetTrack()->SetTrackStatus( fStopAndKill );
              break;
            }
            default:
              aStep->GetTrack()->SetTrackStatus( fStopAndKill );
	      return true;
          } //end of swich cases
        } //end of optical photon
	else {return true;}
      } //end of Cherenkov fiber

      // We are going to create an hit per each fiber with a signal above 0
      // Each fiber is identified with a unique volID 
      //
      Geant4Calorimeter::Hit* hit  = coll->findByKey<Geant4Calorimeter::Hit>(VolID); // the hit
      if(!hit){ // if the hit does not exist yet, create it
        hit = new Geant4Calorimeter::Hit();
        hit->cellID = VolID; // this should be assigned only once
	hit->flag = CherenkovID; // this should be assigned only once
	G4ThreeVector SiPMVec = DREndcapTubesSglHpr::CalculateSiPMPosition(aStep);
	Position SiPMPos(SiPMVec.x(),SiPMVec.y(),SiPMVec.z());
	hit->position = SiPMPos; // this should be assigned only once
        // Note, when the hit is saved in edm4hep format the energyDeposit is
        // divided by 1000, i.e. it translates from MeV (Geant4 unit) to GeV (DD4hep unit).
        // Here I am using this field to save photo-electrons, so I multiply it by 1000
        hit->energyDeposit = signalhit*1000;
        coll->add(VolID, hit); // add the hit to the hit collection
      }
      else { // if the hit exists already, increment its fields
        hit->energyDeposit += signalhit*1000;
      }
 
      /*Position  prePos    = h.prePos();
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
    typedef Geant4SensitiveAction<DREndcapTubesSDData> DREndcapTubesSDAction;
  }}
DECLARE_GEANT4SENSITIVE(DREndcapTubesSDAction)

//**************************************************************************
