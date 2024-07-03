// -------------------------------------------------------------------------
// File: MyDetector.cpp
//
// Purpose: C++ implementation of a very simple dd4hep detector MyDetector:
// a sphere of air inside a box of water
//
// Author: Lorenzo Pezzotti, wildly using dd4hep examples
//
// Created: 2/5/2024
// -------------------------------------------------------------------------

// Include files
#include <DD4hep/DetFactoryHelper.h>
#include <DD4hep/OpticalSurfaces.h>
#include <DD4hep/Printout.h>
#include <DD4hep/Detector.h>

// C/C++ include files
#include <iostream>
#include <map>

// Includers from project files
#include "DREndcapTube.hh"

using namespace dd4hep;
using namespace dd4hep::detail;

// A sphere a air inside a box of water
static Ref_t create_detector(Detector &description, xml_h e, SensitiveDetector /* sens */)  {

  std::cout<<"--> MyDetector::create_detector() start"<<std::endl;

  // Get MyDetector info coded in MyDetector.xml file
  xml_det_t   x_det    = e;
  std::string det_name = x_det.nameStr();
  std::cout<<"--> Getting detector name from xml file: "<<det_name<<", with ID: "<<x_det.id()<<std::endl;

  // Get xml info also for subdetectors, tank and bubble
  xml_det_t   x_tank   = x_det.child(_Unicode(tank));

  // Create the box of water that contains the bubble
  std::cout<<"--> Going to create a box of water sphere with dimensions: "
      <<x_tank.x()/mm<<" x(mm), "<<x_tank.y()/mm<<" y(mm), "<<x_tank.z()/mm<<" z(mm)"<<std::endl;
  // Equivalent to G4Box
  Box    tank_box(x_tank.x(), x_tank.y(), x_tank.z());
  // Equivalent to G4LogicalVolume
  Volume tank_vol("Tank",tank_box,description.material(x_tank.attr<std::string>(_U(material))));
  // Assign vis attributes to the water box (logical) volume
  tank_vol.setVisAttributes(description, x_tank.visStr());

  //Part of endcap container
  const double innerR = 2500.*mm; // inner radius
  const double tower_height = 2000.*mm; // tower height/length
  const double length = tower_height; // alias for tower height/length
  const int NbOfZRot = 36; // number or rations aroung Z axis
  const double phi_unit = 2*M_PI/(double)NbOfZRot; // slice delta phi
  // number of eta towers for the barrel
  const int NbOfBarrel = 40;
  // number of eta towers for the endcap
  int NbOfEndcap = NbOfBarrel-1;

  // Logical volume that contains a slice of the endcap r
  // (I use the G4Generictrap class, so I define the 8 points needed)
  //
  // The first two points of the inner face collide on the beam pipe into (0,0)
  // The second two points of the inner face are at (x=tan(0.5*phi_unit, y=innerR)
  // x with plus or minus sign
  double vertices[16];
  vertices[0] = static_cast<double>(0.);
  vertices[1] = static_cast<double>(0.);
  vertices[2] = static_cast<double>(0.);
  vertices[3] = static_cast<double>(0.);
  vertices[4] = static_cast<double>(-innerR*tan(0.5*phi_unit));
  vertices[5] = static_cast<double>(innerR);
  vertices[6] = static_cast<double>(innerR*tan(0.5*phi_unit));
  vertices[7] = static_cast<double>(innerR);
  // The first two points of the outer face collide on the beam pipe into (0,0)
  // The second two poits of the outer face are same as before with innerR+tower_height
  vertices[8] = static_cast<double>(0.);
  vertices[9] = static_cast<double>(0.);
  vertices[10] = static_cast<double>(0.);
  vertices[11] = static_cast<double>(0.);
  vertices[12] = static_cast<double>(-(innerR+tower_height)*tan(0.5*phi_unit));
  vertices[13] = static_cast<double>(innerR+tower_height);
  vertices[14] = static_cast<double>((innerR+tower_height)*tan(0.5*phi_unit));
  vertices[15] = static_cast<double>(innerR+tower_height);
  // Equivalent of Geant4 GenericTrap shape constructor
  EightPointSolid phiER( "phiER", tower_height/2., vertices);
  // Equivalent to Geant4 Logical Volume
  Volume phiERLog("phiEL", phiER,description.material(x_tank.attr<std::string>(_U(material))));
  phiERLog.setVisAttributes(description, x_tank.visStr());

  // Place logical volume containing the Endcap R slice multiple times
  //
  double deltatheta_barrel[40] = {0};
  // Assign delta theta of each barrel tower
  for(int i=0;i<NbOfBarrel;i++) deltatheta_barrel[i] = M_PI/4/(NbOfBarrel);
  double thetaB = 0; // add up max theta at the end of barrel
  for(int i=0;i<NbOfBarrel;i++) thetaB += deltatheta_barrel[i];
  // Rotate the endcap R slice around the Z axis
  //for(int j=0;j<NbOfZRot;j++){
  for(int j=0;j<NbOfZRot;j++){

    RotationZ rotz1(M_PI/2.); // here I discovered that rotations around Z are inverted
    RotationZ rotz2(j*phi_unit); // w.r.t. Geant4 (+->-)
    RotationX rotx(M_PI);
    RotationY roty(M_PI);
    Transform3D slice_trnsform(rotz1*rotz2*rotx*roty, Position(0,0,(innerR)*tan(thetaB)+length/2.)); 

    // Place each slice into the World
    // I move it along Z of innerR + hald tower length
    //new G4PVPlacement(rmER,G4ThreeVector(0,0,(innerR)*tan(thetaB)+length/2.),phiERLog,"phiERPhys",logicWorld,false,j,false);
    PlacedVolume phiERPlaced = tank_vol.placeVolume(phiERLog,j,slice_trnsform);
    phiERPlaced.addPhysVolID("phiERPlace",j);
  }

  // Create a (DetElement) corresponding to MyDetector.
  // From DD4hep docs https://dd4hep.web.cern.ch/dd4hep/usermanuals/DD4hepManual/DD4hepManualch2.html
  // "Construct the main detector element of this subdetector.This will be the unique entry point
  //  to access any information of the subdetector."
  DetElement  sdet(det_name,x_det.id());
  // Then "Place the subdetector envelope into its mother (typically the top level (world) volume)."
  Volume motherVolume = description.pickMotherVolume(sdet);
  // Place the water box (tank) inside the mother volume
  PlacedVolume tankPlace = motherVolume.placeVolume(tank_vol);
  tankPlace.addPhysVolID("tank",x_det.id());
  sdet.setPlacement(tankPlace);

  std::cout<<"--> MyDetector::create_detector() end"<<std::endl;

  return sdet;
}

DECLARE_DETELEMENT(DD4hep_MyDetector,create_detector)
