//**************************************************************************
// \file DREndcapTubesRunAction.cpp
// \brief:  implementation of DREndcapTubesRunAction class
// \author: Lorenzo Pezzotti (CERN) @lopezzot
// \start date: 12 August 2024
//**************************************************************************

// Includers from project files
//
#include "DREndcapTubesRunAction.hh"

// Includers from Geant4
//
#include "G4AnalysisManager.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

// Includers from C++
//
#include <iostream>
#include <string>

namespace dd4hep
{
namespace sim
{

// Define constructor
//
DREndcapTubesRunAction::DREndcapTubesRunAction(/*DREndcapTubesEventAction* eventAction*/)
  : G4UserRunAction() /*, fEventAction(eventAction)*/
{
  // Instantiate analysis manager
  //
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetVerboseLevel(1);
  // analysisManager->SetNtupleMerging( 1 ); not available in sequential mode

  // Create ntuple
  //
  analysisManager->CreateNtuple("DREndcapTubesout", "DREndcapTubesoutput");
  analysisManager->CreateNtupleDColumn("EnergyScin");  // 0
  analysisManager->CreateNtupleIColumn("SglScin");  // 1
  analysisManager->CreateNtupleDColumn("EnergyCher");  // 2
  analysisManager->CreateNtupleIColumn("SglCher");  // 3
  analysisManager->FinishNtuple();
}  // end of constructor

// Define BeginOfRunAction() and EndOfRunAction() methods
//
void DREndcapTubesRunAction::BeginOfRunAction(const G4Run* Run)
{
  // Create output file, one per Run
  //
  auto analysisManager = G4AnalysisManager::Instance();
  std::string runnumber = std::to_string(Run->GetRunID());
  // Since G4-11.0 the output file format is deducted from extension (.root)
  G4String outputfile = "DREndcapTubesout_Run" + runnumber + ".root";
  analysisManager->OpenFile(outputfile);

  fTimer.Start();  // start timer for event processing
}

void DREndcapTubesRunAction::EndOfRunAction(const G4Run* aRun)
{
  // Stop timer and printout time
  //
  fTimer.Stop();
  G4int events = aRun->GetNumberOfEvent();
  G4cout << "====================================================================== " << G4endl;
  G4cout << " --> DREndcapTubes: run terminated, " << events << " events transported" << G4endl;
  G4cout << " --> DREndcapTubes: time: " << fTimer << G4endl;
  G4cout << " --> DREndcapTubes: time per event: "
         << fTimer.GetUserElapsed() / static_cast<double>(events) << G4endl;
  G4cout << "====================================================================== " << G4endl;

  // Write and close root output file
  //
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();
}

}  // namespace sim
}  // namespace dd4hep

//**************************************************************************
