//**************************************************************************
// \file MyDetector.cpp
// \brief: Implementation of the endcap geometry of the IDEA dual-readout
//         calorimeter using the capillary tubes technology
// \author: Lorenzo Pezzotti (CERN) @lopezzot
// \start date: 12 July 2024
//**************************************************************************

// Includers from DD4hep
#include <DD4hep/DetFactoryHelper.h>
#include "DDRec/Vector3D.h"

// Includers from stl
#include <iostream>
#include <array>

using namespace dd4hep;
using namespace dd4hep::rec; // for dd4hep::rec::Vector3D

// Includers from project files
#include "DREndcapTube.hh"

// Create the endcap calorimeter
//
static Ref_t create_detector(Detector &description, xml_h e, SensitiveDetector /* sens */)  {

  std::cout<<"--> DREndcapTube::create_detector() start"<<std::endl;

  // Get info from the xml file
  // 
  // Info for the main detector element DREndcapTube
  xml_det_t   x_det    = e;
  std::string det_name = x_det.nameStr(); // DREndcapTube volume
  std::cout<<"--> Going to create "<<det_name<<", with ID: "<<x_det.id()<<std::endl;
  xml_dim_t x_dim = x_det.dimensions();
  const double innerR = x_dim.inner_radius();   // inner radius at theta = 0. rad
  const double tower_height = x_dim.z_length(); // tower height/length
  const double length = tower_height;           // alias for tower height/length
  const int NbOfZRot = static_cast<int>(x_dim.deltaphi()); // number or rations aroung Z axis
  const double phi_unit = 2*M_PI/(double)NbOfZRot; // slice delta phi
  std::cout<<"--> From XML description: innerR "<<innerR/m<<" m, tower length "<<tower_height/m<<
      " m, z-rotations "<<NbOfZRot;
  // Info for subdetectors
  xml_det_t   x_tank   = x_det.child(_Unicode(tank));
  xml_det_t   x_stave  = x_det.child(_Unicode(stave));
  xml_det_t   x_tower  = x_det.child(_Unicode(tower));
  xml_det_t   x_capillary_S = x_det.child(_Unicode(tube_S));
  xml_det_t   x_capillary_C = x_det.child(_Unicode(tube_C));
  const double tubeRadius = x_capillary_S.outer_radius();
  std::cout<<", tube radius "<<tubeRadius/mm<<" mm"<<std::endl;

  // Hardcoded parameters (not supposed to be changed)
  constexpr double thetaB{M_PI/4}; // theta at the end of barrel (45 deg)
  constexpr int NbOfEndcap{39};    // number of eta towers for the endcap.
                                   // De facto I will stop at 35 to leave
                                   // room for the beam pipe.
  // Length correction to be applied to capillary tubes
  // with length reduced by the intersection
  // with the tower's planes.
  // It is needed for the intersections with
  // the lateral sides for which the intersecting
  // point is not exactly the left or right edge of the tube
  constexpr double TubeLengthOffset{0.1*mm};
  const double tubeDiameter = tubeRadius*2.;
  const double y_pitch = tubeDiameter*sqrt(3)/2.; // y-distance from tube center to center
                                                  // fixed by the tubes positioning
  // The 8 vertices of a tower to be defined later
  Vector3D pt[8];
  // Create the geometry helper and set paramters
  DREndcapTubeHelper Helper;
  Helper.SetInnerR(innerR);
  Helper.SetTowerHeight(tower_height);
  Helper.SetNumZRot(NbOfZRot);
  Helper.SetTubeRadius(tubeRadius);

  // Start building the geometry
  //

  // Create a tank to place the calorimeter
  Box    tank_box(x_tank.x(), x_tank.y(), x_tank.z());
  Volume tank_vol("Tank",tank_box,description.material(x_tank.attr<std::string>(_U(material))));
  tank_vol.setVisAttributes(description, x_tank.visStr());

  // Volume that contains a slice of the right endcap
  // I use the EightPointSolid/Arb8/G4Generictrap so I define directly its 8 points
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
  EightPointSolid phiER("phiER", tower_height/2., vertices);
  Volume phiERLog("phiER", phiER,description.material(x_stave.attr<std::string>(_U(material))));
  phiERLog.setVisAttributes(description, x_stave.visStr());

  // Place logical volume containing the Endcap R slice multiple times
  //
  // Rotate the endcap R slice around the Z axis
  for(std::size_t j=0;j<static_cast<std::size_t>(NbOfZRot);j++){

    RotationZ rotz1(M_PI/2.);    // here I discovered that dd4hep rotations around Z are inverted
    RotationZ rotz2(j*phi_unit); // w.r.t. Geant4 (+->-)
    RotationX rotx(M_PI);
    RotationY roty(M_PI);
    Transform3D slice_trnsform(rotz1*rotz2*rotx*roty, Position(0,0,(innerR)*tan(thetaB)+length/2.)); 

    //PlacedVolume phiERPlaced = tank_vol.placeVolume(phiERLog,j,slice_trnsform);
    //phiERPlaced.addPhysVolID("phiERPlace",j);
  } // end of slice/stave placement

  // Create an S tube with full tower length
  Tube capillary_S(0.*mm, tubeRadius, tower_height/2., 2*M_PI);
  Volume capillary_SLog("capillary_SLog",capillary_S,description.material(x_capillary_S.attr<std::string>(_U(material))));
  capillary_SLog.setVisAttributes(description, x_capillary_S.visStr());

  // Create a C tube with full tower length
  Tube capillary_C(0.*mm, tubeRadius, tower_height/2., 2*M_PI);
  Volume capillary_CLog("capillary_CLog",capillary_C,description.material(x_capillary_C.attr<std::string>(_U(material))));
  capillary_CLog.setVisAttributes(description, x_capillary_C.visStr());

  // Build the towers inside and endcap R slice
  //
  Helper.SetRbool(1); // Build the right endcap (will reflect it for left part)
  double thetaofcenter = 0.; // theta angle to center of tower being constructed
  double thetaofcenter2 = 0.; // theta angle to center of next tower
  // I fix delta theta of every tower to be the same for every tower
  double deltatheta_endcap[40] = {0.};
  for(int i=0;i<NbOfEndcap+1;i++) deltatheta_endcap[i] = M_PI/4/(NbOfEndcap+1);
  double thetaE = 0.;
  for(int i=0;i<NbOfEndcap+1;i++) thetaE += deltatheta_endcap[i];
  double fulltheta = thetaE; // 45 deg by construction

  // Loop over tower number, stop at tower 35 to leave room for beam pipe
  //for(int i=0;i<NbOfEndcap-1;i++){
  for(int i=0;i<NbOfEndcap-1;i++){

    // Center of first (highest) tower is 45 deg - deltatheta_endcap[1]/2
    thetaofcenter = fulltheta-deltatheta_endcap[i]/2.;
    // Center of the next tower
    thetaofcenter2=thetaofcenter-deltatheta_endcap[i]/2.-deltatheta_endcap[i+1]/2.;
    // Update Helper class parameters accordingly
    Helper.SetDeltaTheta(deltatheta_endcap[i]);
    Helper.SetThetaOfCenter(thetaofcenter);
    Helper.SetDeltaTheta2(deltatheta_endcap[i+1]);
    Helper.SetThetaOfCenter2(thetaofcenter2);
    Helper.CalBasic(); // Perform internal calculations
    Helper.Getpt(pt);  // Update 8 Vectors defining the tower edges
    
    // Create now the tower as a Trap
    //
    // Problem: G4Trap has a constructors using the 8 vertices of the trapezoid
    // but DD4hep Trap does not have such constructor.
    // Therefore I calculate here the parameters to be passed to the DD4hep Trap. 
    auto pDz        = (pt[7]).z();
    auto pDy1       = ((pt[2]).y()-(pt[1]).y())*0.5;
    auto pDx1       = ((pt[1]).x()-(pt[0]).x())*0.5;
    auto pDx2       = ((pt[3]).x()-(pt[2]).x())*0.5;
    auto fTalpha1   = ((pt[2]).x()+(pt[3]).x()-(pt[1]).x()-(pt[0]).x())*0.25/pDy1;
    auto pAlp1      = std::atan(fTalpha1);
    auto pDy2       = ((pt[6]).y()-(pt[5]).y())*0.5;
    auto pDx3       = ((pt[5]).x()-(pt[4]).x())*0.5;
    auto pDx4       = ((pt[7]).x()-(pt[6]).x())*0.5;
    auto fTalpha2   = ((pt[6]).x()+(pt[7]).x()-(pt[5]).x()-(pt[4]).x())*0.25/pDy2;
    auto pAlp2      = std::atan(fTalpha2);
    auto fThetaCphi = ((pt[4]).x()+pDy2*fTalpha2+pDx3)/pDz;
    auto fThetaSphi = ((pt[4]).y()+pDy2)/pDz;
    double pPhi     = 0.;
    double pTheta   = 0.;
    if(fThetaSphi==0. && fThetaCphi==0.){}
    else{
         pPhi = std::atan(fThetaSphi/fThetaCphi);
         pTheta = std::atan(fThetaCphi/std::cos(pPhi));
    } // end of Trap parameters calculation
    
    std::cout<<"--> Tower "<<i<<" being constructed"<<std::endl;

    Trap tower("phiER", pDz, pTheta, pPhi, pDy1, pDx1, pDx2, pAlp1, pDy2, pDx3, pDx4, pAlp2);
    Volume towerLog("towerLog", tower, description.material(x_tower.attr<std::string>(_U(material))));
    towerLog.setVisAttributes(description, x_tower.visStr());

    RotationX rotX(-thetaofcenter);   // Rotation matrix for this tower
    RotationY rotY(0.);
    RotationZ rotZ(0.);
    Vector3D c = Helper.GetOrigin(0); // Get origin of tower
    Vector3D c_new(-c.y(),c.z(),c.x()-(innerR+0.5*length));
    //if(i<35) { // Place towers up to 35, this "if" is just a protection in case the loop range is changed
    //    Transform3D tower_trnsform(rotX*rotY*rotZ, Position(c_new.x(),c_new.y(),c_new.z())); 
    //    PlacedVolume towerPlaced = phiERLog.placeVolume(towerLog,i,tower_trnsform);
    //    towerPlaced.addPhysVolID("tower",i);
    //}
    // Or, to debug, place towers one next to each other in tank volume
    if(i<35) {
        double z = static_cast<int>(i/15)*(length+40*cm);
        double x = (i-static_cast<int>(i/15)*15)*100*cm - 5*m;
        tank_vol.placeVolume(towerLog,i,Position(-1.*x,0.,-1.*z));
    }

    // Capillary placement inside tower (both S and C)
    //
    // Fill a tower along y (vertical direction)
    // per each row calculate x and y of the starting tube (left most)
    const double y_backplane = pt[4].y(); // y-coordinate of bottom left corner of back face
    const double x_backplane = pt[4].x(); // x-coordinate of bottom left corner of back face
    for(std::size_t k=0; k<200; k++){
      double x_start = x_backplane;
      double y_tube = 0.;
      double delta_x = 0.; // correction to tune x-coordinate
      // if it is the first row fix the y coordinate
      if(k==0) y_tube = y_backplane + tubeRadius + 1.0*mm; //add 1 mm otherwise the tube immediately intersect with a plane
      else {
        y_tube = y_backplane + tubeRadius + 1.0*mm + k*2.*y_pitch;

        // Adjust starting x given the opening angle of the trapezoidal back face
        double hypotenuse = sqrt(pow((-1*y_backplane)*2,2) + pow(pt[6].x()-pt[4].x(),2));
        double angle = acos(((-1*y_backplane)*2) / hypotenuse);
        delta_x = ((-1.*y_backplane) - (-1.*y_tube)) * tan(angle);
      }
      if(delta_x > tubeDiameter) {
          // Calculate how many fibers should be inserted in x.
          // Caveat: casting to int does not round takes only the 
          // integer part (i.e. 6.9 -> 6)
          auto newTubesNo = static_cast<int>(delta_x/tubeDiameter);
          x_start = x_start - newTubesNo*tubeDiameter;
      }

      //Fill a tower row along x (lateral direction)
      for(std::size_t j=0;j<1000;j++){
         auto x_tube = x_start  + tubeRadius + j*(tubeDiameter);       // locate the center of the tube in x 
         Vector3D capillaryPos(x_tube, y_tube, length/2.);             // locate tube on tower back face
         auto capillaryLength = Helper.GetTubeLength(pt,capillaryPos); // calculate tube length
         if(capillaryLength == length){
           PlacedVolume capillaryPlaced = towerLog.placeVolume(capillary_SLog, 1000*k+j, Position(x_tube, y_tube, 0.));
           capillaryPlaced.addPhysVolID("capillary_S", 1000*k+j);
         }
         else if(capillaryLength > 5.0*cm){
           // Note: the root visualizer does not display tubes with capillaryLength < 4.5 cm.
           // Such tubes are usually the ones closest to the left or right side of a tower.
           // They seem to be placed correctly but not displayed. 
           // I am adding a cut over tube length of 5.0 cm: tubes below it will not be placed.
           Tube capillaryShort(0.*mm, tubeRadius, (capillaryLength-TubeLengthOffset)/2., 2*M_PI); // reduced capillary length
                                                                                                  // by a fixed offset
           Volume capillaryShortLog("capillaryShortLog",capillaryShort,description.material(x_capillary_S.attr<std::string>(_U(material))));
           capillaryShortLog.setVisAttributes(description, x_capillary_S.visStr());
           PlacedVolume capillaryShortPlaced = towerLog.placeVolume(capillaryShortLog, 1000*k+j, Position(x_tube, y_tube, length/2.-capillaryLength/2.+TubeLengthOffset/2.));
           capillaryShortPlaced.addPhysVolID("capillary_S", 1000*k+j);
         }
         else {;}

         // Now C tubes are placed.
         // First check there is enough room in y-direction, if not do not exit tube placement
         if((-1.*y_backplane)-y_tube < (y_pitch+tubeRadius)) break; // do not place C tube row
         // This checks if there is enough room in x-direction to place a C tube after its S one
         bool IsLastTube_C = -1.*x_backplane+delta_x-x_tube < tubeDiameter ? true : false;

         // After the S tube placement I place the closest C tube
         // according to the fixed structure of the tubes placement (gluing)
         //
         // If the S tube below was not placed (too short) do not place the C either
         // && if there is not enough room on x-direction do not place the C tube
         if (capillaryLength > 5.0*cm && !IsLastTube_C){
           double x_tube_C = x_tube + tubeRadius;
           double y_tube_C = y_tube + y_pitch; 
           Vector3D capillaryPos_C(x_tube_C, y_tube_C, length/2.);
           auto capillaryLength_C = Helper.GetTubeLength(pt,capillaryPos_C);
           if(capillaryLength_C == length){
             PlacedVolume capillaryPlaced_C = towerLog.placeVolume(capillary_CLog, 100000*k+j, Position(x_tube_C, y_tube_C, 0.));
             capillaryPlaced_C.addPhysVolID("capillary_C", 100000*k+j);
           } 
           else if(capillaryLength_C > 5.0*cm){
               Tube capillaryShort_C(0.*mm, tubeRadius, (capillaryLength_C-TubeLengthOffset)/2., 2*M_PI); // reduced capillary length
                                                                                                          // by a fixed offset
               Volume capillaryShortLog_C("capillaryShortLog_C",capillaryShort_C,description.material(x_capillary_C.attr<std::string>(_U(material))));
               capillaryShortLog_C.setVisAttributes(description, x_capillary_C.visStr());
               PlacedVolume capillaryShortPlaced_C = towerLog.placeVolume(capillaryShortLog_C, 100000*k+j, Position(x_tube_C, y_tube_C, length/2.-capillaryLength_C/2.+TubeLengthOffset/2.));
               capillaryShortPlaced_C.addPhysVolID("capillary_C", 100000*k+j);
           }
           else {;}
         }
           
         // condition for stopping S capillary placement along x
         if(-1.*x_backplane+delta_x-x_tube < tubeDiameter+tubeRadius) break;
      } // end x loop
 
      // Condition for stopping C capillary placement along y.
      // y_backplane is equal up and down so I can keep the same for exiting loop
      if((-1.*y_backplane)-y_tube < (2.*y_pitch+tubeRadius)) break;
    } // End y loop and tube placement

    // Update parameters
    Helper.Getpt(pt);
    fulltheta = fulltheta-deltatheta_endcap[i];
  } //End of towers creation and placement

  // Create a (DetElement) corresponding to MyDetector.
  // From DD4hep docs https://dd4hep.web.cern.ch/dd4hep/usermanuals/DD4hepManual/DD4hepManualch2.html
  // "Construct the main detector element of this subdetector.This will be the unique entry point
  // to access any information of the subdetector."
  DetElement  sdet(det_name,x_det.id());
  // Then "Place the subdetector envelope into its mother (typically the top level (world) volume)."
  Volume motherVolume = description.pickMotherVolume(sdet);
  // Place the tank container inside the mother volume
  PlacedVolume tankPlace = motherVolume.placeVolume(tank_vol);
  tankPlace.addPhysVolID("tank",x_det.id());
  sdet.setPlacement(tankPlace);

  std::cout<<"--> DREndcapTube::create_detector() end"<<std::endl;

  return sdet;
}

DECLARE_DETELEMENT(DD4hep_MyDetector,create_detector)

//**************************************************************************
