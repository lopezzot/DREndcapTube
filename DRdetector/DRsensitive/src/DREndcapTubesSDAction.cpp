//**************************************************************************
// \file DREndcapTubesSDAction.cpp
// \brief: Implementation of Geant4SensitiveAction template class for
//         DREndcapTubes
// \author: Lorenzo Pezzotti (CERN) @lopezzot
// \start date: 11 August 2024
//**************************************************************************

// Includers from DD4HEP
#include "DD4hep/Segmentations.h"
#include "DD4hep/Version.h"
#include "DDG4/Factories.h"
#include "DDG4/Geant4EventAction.h"
#include "DDG4/Geant4GeneratorAction.h"
#include "DDG4/Geant4Mapping.h"
#include "DDG4/Geant4RunAction.h"
#include "DDG4/Geant4SensDetAction.inl"

// Includers from Geant4
#include "G4OpBoundaryProcess.hh"
#include "G4OpticalPhoton.hh"
#include "G4Poisson.hh"
#include "G4ProcessManager.hh"
#include "G4Types.hh"
#include "G4UserSteppingAction.hh"
#include "G4VProcess.hh"
#include "globals.hh"

// Includers from project files
#include "DREndcapTubesEvtAction.hh"
#include "DREndcapTubesRunAction.hh"
#include "DREndcapTubesSglHpr.hh"

// #define DREndcapTubesSDDebug

namespace dd4hep
{
namespace sim
{
class DREndcapTubesSDData
{
    // Constructor and destructor
    //
  public:
    DREndcapTubesSDData() = default;
    ~DREndcapTubesSDData() = default;

    // Methods
    //
  public:
    void beginRun(const G4Run* run)
    {
      fRunAction = new DREndcapTubesRunAction;
      fEvtAction = new DREndcapTubesEvtAction;
      fRunAction->BeginOfRunAction(run);
    }

    void endRun(const G4Run* run) { fRunAction->EndOfRunAction(run); }

    void beginEvent(const G4Event* event) { fEvtAction->BeginOfEventAction(event); }

    void endEvent(const G4Event* event) { fEvtAction->EndOfEventAction(event); }

    // Fields
    //
  public:
    DREndcapTubesRunAction* fRunAction;
    DREndcapTubesEvtAction* fEvtAction;
    Geant4Sensitive* sensitive{};
    int collection_cher_right;
    int collection_cher_left;
    int collection_scin_left;
};
}  // namespace sim
}  // namespace dd4hep

namespace dd4hep
{
namespace sim
{

// Function template specialization of Geant4SensitiveAction class.
// Define actions
template<>
void Geant4SensitiveAction<DREndcapTubesSDData>::initialize()
{
  eventAction().callAtBegin(&m_userData, &DREndcapTubesSDData::beginEvent);
  eventAction().callAtEnd(&m_userData, &DREndcapTubesSDData::endEvent);

  runAction().callAtBegin(&m_userData, &DREndcapTubesSDData::beginRun);
  runAction().callAtEnd(&m_userData, &DREndcapTubesSDData::endRun);

  m_userData.sensitive = this;
}

// Function template specialization of Geant4SensitiveAction class.
// Define collections created by this sensitivie action object
template<>
void Geant4SensitiveAction<DREndcapTubesSDData>::defineCollections()
{
  std::string ROname = m_sensitive.readout().name();
  m_collectionID = defineCollection<Geant4Calorimeter::Hit>(ROname + "ScinRight");
  m_userData.collection_cher_right = defineCollection<Geant4Calorimeter::Hit>(ROname + "CherRight");
  m_userData.collection_scin_left = defineCollection<Geant4Calorimeter::Hit>(ROname + "ScinLeft");
  m_userData.collection_cher_left = defineCollection<Geant4Calorimeter::Hit>(ROname + "CherLeft");
}

// Function template specialization of Geant4SensitiveAction class.
// Method that accesses the G4Step object at each track step.
template<>
bool Geant4SensitiveAction<DREndcapTubesSDData>::process(const G4Step* aStep,
                                                         G4TouchableHistory* /*history*/)
{
  // NOTE: Here we do manipulation of the signal in each fiber (Scintillating and Cherenkov)
  // to compute the calorimeter signal and populate the corresponding hit.
  // Sensitive volumes are associated in the steering file using the DD4hep regexSensitiveDetector
  // and matching the substring DRETS. Usually DD4hep volIDs are retrieved with something like this
  // ---
  // dd4hep::BitFieldCoder
  // decoder("tank:1,endcap:1,stave:10,tower:8,air:1,col:10,row:7,clad:1,core:1,cherenkov:1"); auto
  // VolID = volumeID(aStep); auto CherenkovID = decoder.get(VolID,"cherenkov");
  // ---
  // However using regexSensitiveDetector does not populate the cache that allows using the
  // volumeID() method. One can still access volIDs in this sensitive detector action with something
  // like this
  // ---
  // G4VPhysicalVolume* PhysVol = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  // Geant4Mapping& Mapping = Geant4Mapping::instance();
  // PlacedVolume PlacedVol = Mapping.placement(PhysVol);
  // const PlacedVolumeExtension::VolIDs& TestIDs = PlacedVol.volIDs();
  // auto it = TestIDs.find("name");
  // std::cout<<it->first<<" "<<it->second<<std::endl;
  // ---
  // but this brute force method makes the simulaiton slower by more than two order of magnitudes
  // (see https://github.com/AIDASoft/DD4hep/issues/1319).
  // Therefore we use copynumbers instead of volIDs. The copynumbers used are
  // Scintillating core (0), Cherenkov core (2), Cherenkov cladding (3), Tube (1000*row+column)
  // and (phi)Stave (NbOfZRot).

  auto Edep = aStep->GetTotalEnergyDeposit();
  auto cpNo = aStep->GetPreStepPoint()->GetTouchable()->GetCopyNumber();
  bool IsScin = (cpNo == 0);
  bool IsCher = (cpNo == 2);
  bool IsCherClad = (cpNo == 3);
  bool IsRight = (aStep->GetPreStepPoint()->GetPosition().z() > 0.);

  // Skip this step if edep is 0 and it is a scintillating fiber
  if (IsScin && Edep == 0.) return true;
  // If it is a track inside the cherenkov CLADDING skip this step,
  // if it is an optical photon kill it first
  if (IsCherClad) {
    if (aStep->GetTrack()->GetParticleDefinition() == G4OpticalPhoton::Definition()) {
      aStep->GetTrack()->SetTrackStatus(fStopAndKill);
    }
    return true;
  }

  auto TubeID = aStep->GetPreStepPoint()->GetTouchable()->GetCopyNumber(2);
  auto TowerID = aStep->GetPreStepPoint()->GetTouchable()->GetCopyNumber(3);
  auto StaveID = aStep->GetPreStepPoint()->GetTouchable()->GetCopyNumber(4);
  int VolID = 15000000 * StaveID + 300000 * TowerID + TubeID;

#ifdef DREndcapTubesSDDebug
  // Print out some info step-by-step in sensitive volumes
  //
  DREndcapTubesSglHpr::PrintStepInfo(aStep);
#endif

  // We now calculate the signal in S and C fiber according to the step contribution
  //
  G4double steplength = aStep->GetStepLength();
  G4int signalhit = 0;
  Geant4HitCollection* coll = (IsRight && IsScin)    ? collection(m_collectionID)
                              : (IsRight && !IsScin) ? collection(m_userData.collection_cher_right)
                              : (!IsRight && IsScin) ? collection(m_userData.collection_scin_left)
                                                     : collection(m_userData.collection_cher_left);

  if (!IsCher) {  // it is a scintillating fiber

    m_userData.fEvtAction->AddEdepScin(Edep);

    if (aStep->GetTrack()->GetDefinition()->GetPDGCharge() == 0 || steplength == 0.) {
      return true;  // not ionizing particle
    }
    G4double distance_to_sipm = DREndcapTubesSglHpr::GetDistanceToSiPM(aStep);
    signalhit =
      DREndcapTubesSglHpr::SmearSSignal(DREndcapTubesSglHpr::ApplyBirks(Edep, steplength));
    signalhit = DREndcapTubesSglHpr::AttenuateSSignal(signalhit, distance_to_sipm);
    if (signalhit == 0) return true;
    m_userData.fEvtAction->AddSglScin(signalhit);
  }  // end of scintillating fibre sigal calculation

  else {  // it is a Cherenkov fiber
    // save mc truth info in analysismanager auxiliary outputfile
    m_userData.fEvtAction->AddEdepCher(Edep);
    // calculate the signal in terms of Cherenkov photo-electrons
    if (aStep->GetTrack()->GetParticleDefinition() == G4OpticalPhoton::Definition()) {
      G4OpBoundaryProcessStatus theStatus = Undefined;
      G4ProcessManager* OpManager = G4OpticalPhoton::OpticalPhoton()->GetProcessManager();

      if (OpManager) {
        G4int MAXofPostStepLoops = OpManager->GetPostStepProcessVector()->entries();
        G4ProcessVector* fPostStepDoItVector = OpManager->GetPostStepProcessVector(typeDoIt);

        for (G4int i = 0; i < MAXofPostStepLoops; i++) {
          G4VProcess* fCurrentProcess = (*fPostStepDoItVector)[i];
          G4OpBoundaryProcess* fOpProcess = dynamic_cast<G4OpBoundaryProcess*>(fCurrentProcess);
          if (fOpProcess) {
            theStatus = fOpProcess->GetStatus();
            break;
          }
        }
      }

      switch (theStatus) {
        case TotalInternalReflection: {
          // Kill Cherenkov photons inside fibers travelling towards the inner tip
          if (!DREndcapTubesSglHpr::IsReflectedForward(aStep)) return true;
          G4double distance_to_sipm = DREndcapTubesSglHpr::GetDistanceToSiPM(aStep);
          G4int c_signal = DREndcapTubesSglHpr::SmearCSignal();
          signalhit = DREndcapTubesSglHpr::AttenuateCSignal(c_signal, distance_to_sipm);
          if (signalhit == 0) return true;
          // save mc truth info in analysismanager auxiliary outputfile
          m_userData.fEvtAction->AddSglCher(signalhit);
          aStep->GetTrack()->SetTrackStatus(fStopAndKill);
          break;
        }
        default:
          aStep->GetTrack()->SetTrackStatus(fStopAndKill);
          return true;
      }  // end of swich cases
    }  // end of optical photon
    else {
      return true;
    }
  }  // end of Cherenkov fiber

  // We are going to create an hit per each fiber with a signal above 0
  // Each fiber is identified with a unique volID
  //
  Geant4Calorimeter::Hit* hit = coll->findByKey<Geant4Calorimeter::Hit>(VolID);  // the hit
  if (!hit) {  // if the hit does not exist yet, create it
    hit = new Geant4Calorimeter::Hit();
    hit->cellID = VolID;  // this should be assigned only once
    G4ThreeVector FiberVec = DREndcapTubesSglHpr::CalculateFiberPosition(aStep);
    Position FiberPos(FiberVec.x(), FiberVec.y(), FiberVec.z());
    hit->position = FiberPos;  // this should be assigned only once
    // Note, when the hit is saved in edm4hep format the energyDeposit is
    // divided by 1000, i.e. it translates from MeV (Geant4 unit) to GeV (DD4hep unit).
    // Here I am using this field to save photo-electrons, so I multiply it by 1000
    hit->energyDeposit = signalhit * 1000;
    coll->add(VolID, hit);  // add the hit to the hit collection
  }
  else {  // if the hit exists already, increment its fields
    hit->energyDeposit += signalhit * 1000;
  }

  return true;
}  // end of Geant4SensitiveAction::process() method specialization

}  // namespace sim
}  // namespace dd4hep

//--- Factory declaration
namespace dd4hep
{
namespace sim
{
typedef Geant4SensitiveAction<DREndcapTubesSDData> DREndcapTubesSDAction;
}
}  // namespace dd4hep
DECLARE_GEANT4SENSITIVE(DREndcapTubesSDAction)

//**************************************************************************
