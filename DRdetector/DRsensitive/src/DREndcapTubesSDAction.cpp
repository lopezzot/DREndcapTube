//**************************************************************************
// \file DREndcapTubesSDAction.cpp
// \brief: Implementation of Geant4SensitiveAction template class for
//         DREndcapTubes
// \author: Lorenzo Pezzotti (CERN) @lopezzot
// \start date: 11 August 2024
//**************************************************************************

// Includers from DD4HEP
#include <DDG4/Geant4SensDetAction.inl>
#include <DDG4/Geant4ParticleInformation.h>
#include <DDG4/Factories.h>

#define DREndcapTubesSDDebug

class DREndcapTubesSD {
    public:
      double eDep = 0.;
};

namespace dd4hep {

  namespace sim   {

    // Function template specialization of Geant4SensitiveAction class.
    // Define collections created by this sensitivie action object
    template <> void Geant4SensitiveAction<DREndcapTubesSD>::defineCollections()    {
      //m_collectionID = declareReadoutFilteredCollection<MyTrackerSD::Hit>();
    }

    // Function template specialization of Geant4SensitiveAction class.
    // Method that accesses the G4Step object at each track step.
    template <> bool Geant4SensitiveAction<DREndcapTubesSD>::process(const G4Step* aStep, G4TouchableHistory* /*hist*/ ) {
    
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
