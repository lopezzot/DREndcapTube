//**************************************************************************
// \file DREndcapTubesSglHpr.hh
// \brief:  definition of DREndcapTubesSglHpr class
// \author: Lorenzo Pezzotti (CERN) @lopezzot
// \start date: 13 August 2024
//**************************************************************************

// Header-only utility class to store methods needed when computing
// the DREndcapTubes signals in scintillating and Cherenkov fibers

#ifndef DREndcapTubesSglHpr_h
#define DREndcapTubesSglHpr_h 1

// Includers from Geant4
//
#include "globals.hh"
#include "G4Tubs.hh"
#include "CLHEP/Units/SystemOfUnits.h"

class DREndcapTubesSglHpr {

  public:
    DREndcapTubesSglHpr() = delete;

  private:
    // Fields
    //
    static constexpr G4double fk_B = 0.126; // Birks constant
    static constexpr G4double fSAttenuationLength = 1000.0*CLHEP::m;
    static constexpr G4double fCAttenuationLength = 1000.0*CLHEP::m;
	
  public:
    // Methods
    //
    // Apply Briks Law
    static constexpr G4double ApplyBirks(const G4double& de, const G4double& steplength){
      return (de/steplength) / ( 1+fk_B*(de/steplength) ) * steplength;
    }

    // Smear S signal according to Poissonian fluctuations and light yield
    static G4int SmearSSignal( const G4double& satde ){
      return G4Poisson(satde*9.5);
    }

    // Smear C signal according to Poissonian fluctuations and light yield
    static G4int SmearCSignal( ){
      return G4Poisson(0.153);
    }
    
    // Calculate distance from step in fiber to SiPM
    inline static G4double GetDistanceToSiPM(const G4Step* step);
    
    // Calculate how many photons survived after light attenuation
    inline static G4int AttenuateHelper(const G4int& signal, const G4double& distance, const G4double& attenuation_length); 

    // Attenuate light in fibers
    static G4int AttenuateSSignal(const G4int& signal, const G4double& distance) {
      return AttenuateHelper(signal, distance, fSAttenuationLength);
    }
    
    static G4int AttenuateCSignal(const G4int& signal, const G4double& distance) {
      return AttenuateHelper(signal, distance, fCAttenuationLength);
    }

    inline static G4ThreeVector CalculateSiPMPosition(const G4Step* step);

};

inline G4double DREndcapTubesSglHpr::GetDistanceToSiPM(const G4Step* step) {

  // Get the pre-step point
  const G4StepPoint* preStepPoint = step->GetPreStepPoint();
  // Get the global position of the pre-step point
  G4ThreeVector globalPos = preStepPoint->GetPosition();
  // Get the local position of the pre-step point in the current volume's coordinate system
  G4ThreeVector localPos = preStepPoint->GetTouchableHandle()->GetHistory()->GetTopTransform().TransformPoint(globalPos);

  // Get the logical volume of the current step
  G4LogicalVolume* currentVolume = preStepPoint->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  // Get the solid associated with the logical volume
  G4Tubs* solid = dynamic_cast<G4Tubs*>(currentVolume->GetSolid());
  // Get the dimensions of the solid (size of the volume)
  G4double size = solid->GetZHalfLength();

  G4double distance_to_sipm = size - localPos.z();
  return distance_to_sipm;
}

inline G4int DREndcapTubesSglHpr::AttenuateHelper(const G4int& signal, const G4double& distance, const G4double& attenuation_length) {
  
  double probability_of_survival = exp(-distance/attenuation_length);

  G4int survived_photons = 0;
  for (int i=0; i<signal; i++) {
    // Simulate drawing between 0 and 1 with probability x of getting 1
    if (G4UniformRand() <= probability_of_survival) survived_photons++;
  }
  return survived_photons;
}

inline G4ThreeVector DREndcapTubesSglHpr::CalculateSiPMPosition(const G4Step* step) {

  G4TouchableHandle theTouchable = step->GetPreStepPoint()->GetTouchableHandle();
  G4ThreeVector origin(0.,0.,0.);
  G4ThreeVector zdir(0.,0.,1.);
  G4ThreeVector vectPos = theTouchable->GetHistory()-> GetTopTransform().Inverse().TransformPoint(origin);
  G4ThreeVector direction = theTouchable->GetHistory()->GetTopTransform().Inverse().TransformAxis(zdir);

  // Get the logical volume of the current step
  G4LogicalVolume* currentVolume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  // Get the solid associated with the logical volume
  G4Tubs* solid = dynamic_cast<G4Tubs*>(currentVolume->GetSolid());
  // Get the dimensions of the solid (size of the volume)
  G4double size = solid->GetZHalfLength();
  G4double lengthfiber = size*2.;
  G4ThreeVector Halffibervect = direction*lengthfiber/2;
  // Fibre tip position
  // G4ThreeVector vectPostip = vectPos-Halffibervect;
  // SiPM position
  G4ThreeVector SiPMvecPos = vectPos+Halffibervect;

  return SiPMvecPos;
}

#endif // DREndcapTubesSglHpr_h

//**************************************************************************