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
#include <cmath>
#include <array>

// Includers from project files
#include "DREndcapTube.hh"

using namespace dd4hep;
using namespace dd4hep::detail;
using namespace dd4hep::rec;

// This struct represents a plane given the 4 vertices
struct Plane{
  Vector3D P1,P2,P3,P4;
  Plane(Vector3D P1, Vector3D P2, Vector3D P3, Vector3D P4) : P1(P1),P2(P2),P3(P3),P4(P4){};
};

// This struct represents a Zline given a point and a fixed direction
struct ZLine{
  Vector3D origin;
  Vector3D fuZ = Vector3D(0,0,-1);
  ZLine(Vector3D P) : origin(P){};
};

// Function to calculate the intersection of a Zline with a plane
bool intersectLinePlane(const ZLine& line, const Plane& plane, Vector3D& intersection) {
    Vector3D p0 = plane.P1;
    Vector3D p1 = plane.P2;
    Vector3D p2 = plane.P3;

    // Compute the plane normal
    Vector3D u = Vector3D(p1.x() - p0.x(), p1.y() - p0.y(), p1.z() - p0.z());
    Vector3D v = Vector3D(p2.x() - p0.x(), p2.y() - p0.y(), p2.z() - p0.z());
    Vector3D normal = Vector3D(u.y() * v.z() - u.z() * v.y(), u.z() * v.x() - u.x() * v.z(), u.x() * v.y() - u.y() * v.x());

    double denominator = normal.x() * line.fuZ.x() + normal.y() * line.fuZ.y() + normal.z() * line.fuZ.z();

    if (std::fabs(denominator) < 1e-6) {
        return false; // The line is parallel to the plane
    }

    double d = normal.x() * p0.x() + normal.y() * p0.y() + normal.z() * p0.z();
    double t = (d - normal.x() * line.origin.x() - normal.y() * line.origin.y() - normal.z() * line.origin.z()) / denominator;

    if (t < 0) {
        return false; // The intersection is behind the line origin
    }

    intersection = Vector3D(line.origin.x() + t * line.fuZ.x(),
                            line.origin.y() + t * line.fuZ.y(),
                            line.origin.z() + t * line.fuZ.z());

    return true;
};

// Function to calculate fiber length
double GetFiberLength(const Vector3D (&pt)[8], const Vector3D& point){
    Plane FirstPlane(pt[0],pt[1],pt[2],pt[3]);
    Plane SecondPlane(pt[1],pt[5],pt[3],pt[7]);
    Plane ThirdPlane(pt[2],pt[3],pt[6],pt[7]);
    Plane FourthPlane(pt[0],pt[4],pt[2],pt[6]);
    Plane FifthPlane(pt[0],pt[1],pt[5],pt[4]);
    std::array<Plane, 5> Planes{FirstPlane, SecondPlane, ThirdPlane, FourthPlane, FifthPlane};
    // I calculate the intersection over four line (up, down, left, right) of the capillary
    ZLine PointZLine(point);
    ZLine UpPointZLine(Vector3D(point.x(),point.y()+2.0*mm,point.z()));
    ZLine DownPointZLine(Vector3D(point.x(),point.y()-2.0*mm,point.z()));
    ZLine LeftPointZLine(Vector3D(point.x()-2.0*mm,point.y(),point.z()));
    ZLine RightPointZLine(Vector3D(point.x()+2.0*mm,point.y(),point.z()));
    std::array<ZLine, 4> Lines{UpPointZLine, DownPointZLine, LeftPointZLine, RightPointZLine};
    // Intersection to be found
    Vector3D intersection;
    double FiberLength = 0.;
    int side = 99.; //which capillary side intersected first
    for(std::size_t k=0; k<Lines.size(); k++){
      for(std::size_t i=0; i<Planes.size(); i++){
        auto truth = intersectLinePlane(Lines.at(k), Planes.at(i), intersection);
        //auto truth = intersectLinePlane(PointZLine, Planes.at(i), intersection);
        double FiberLengthPlane = (Lines.at(k).origin-intersection).r();
        //double FiberLengthPlane = (point-intersection).r();
        //std::cout<<"Plane "<<i<<" FiberLengthPlane "<<FiberLengthPlane<<std::endl;
        if(k==0 && i==0) FiberLength = FiberLengthPlane;
        else if (FiberLength > FiberLengthPlane) FiberLength = FiberLengthPlane;
        else {;} 
        std::cout<<"k "<<k<<" i "<<i<<" "<<FiberLength<<std::endl;
      }
    }
    std::cout<<"point "<<point.x()<<" "<<point.y()<<" final fiber length "<<FiberLength<<std::endl;
    return FiberLength;
}

// A sphere a air inside a box of water
static Ref_t create_detector(Detector &description, xml_h e, SensitiveDetector /* sens */)  {

  std::cout<<"--> MyDetector::create_detector() start"<<std::endl;

  // Get MyDetector info coded in MyDetector.xml file
  xml_det_t   x_det    = e;
  std::string det_name = x_det.nameStr();
  std::cout<<"--> Getting detector name from xml file: "<<det_name<<", with ID: "<<x_det.id()<<std::endl;

  // Get xml info also for subdetectors: tank, stave and tower
  xml_det_t   x_tank   = x_det.child(_Unicode(tank));
  xml_det_t   x_stave  = x_det.child(_Unicode(stave));
  xml_det_t   x_tower  = x_det.child(_Unicode(tower));
  xml_det_t   x_capillary  = x_det.child(_Unicode(capillary));

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

  // Create helper class
  DREndcapTubeHelper Helper;
  Helper.SetInnerR(innerR);
  Helper.SetTower_height(tower_height);
  Helper.SetNumZRot(NbOfZRot);

  dd4hep::rec::Vector3D pt[8]; // dummy inizialization of pt

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
  phiERLog.setVisAttributes(description, x_stave.visStr());

  // Place logical volume containing the Endcap R slice multiple times
  //
  double deltatheta_barrel[40] = {0};
  // Assign delta theta of each barrel tower
  for(int i=0;i<NbOfBarrel;i++) deltatheta_barrel[i] = M_PI/4/(NbOfBarrel);
  double thetaB = 0; // add up max theta at the end of barrel
  for(int i=0;i<NbOfBarrel;i++) thetaB += deltatheta_barrel[i];
  // Rotate the endcap R slice around the Z axis
  for(int j=0;j<NbOfZRot;j++){
  //for(int j=0;j<1;j++){

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

  // Create a brass tube
  const double tubeRadius = 2.0*mm;
  Tube capillary(0.*mm, tubeRadius, tower_height/2., 2*M_PI);
  Volume capillaryLog("capillaryLog",capillary,description.material(x_tank.attr<std::string>(_U(material))));
  capillaryLog.setVisAttributes(description, x_capillary.visStr());

  // Build the towers inside and endcap R slice
  //
  Helper.Rbool(1);
  double thetaofcenter = 0.; // theta angle to center of tower being constructed
  double thetaofcenter2 = 0.; // theta angle to center of next tower
  double deltatheta_endcap[40] = {0.};
  for(int i=0;i<NbOfEndcap+1;i++) deltatheta_endcap[i] = M_PI/4/(NbOfEndcap+1);
  double thetaE = 0.;
  for(int i=0;i<NbOfEndcap+1;i++) thetaE += deltatheta_endcap[i];
  double fulltheta = thetaE; // 45 deg by construction

  //for(int i=0;i<NbOfEndcap-1;i++){
  for(int i=0;i<1;i++){

    // this is theta angle from beam pipe up
    // center of first (highest) tower is 45 deg - deltatheta_endcap[i]
    thetaofcenter = fulltheta-deltatheta_endcap[i]/2.;
    //G4cout<<"thetaofcenter "<<thetaofcenter<<" rad"<<G4endl;
    // this is theta of center of the next tower
    thetaofcenter2=thetaofcenter-deltatheta_endcap[i]/2.-deltatheta_endcap[i+1]/2.;
    //G4cout<<"thetaofcenter2 "<<thetaofcenter2<<" rad"<<G4endl;
    //std::cout<<i<<" deltatheta_endcap_i "<<deltatheta_endcap[i]<<" thetaofcenter "<<thetaofcenter<<
    //    " deltatheta_endcap_i+1 "<<deltatheta_endcap[i+1]<<" thetaofcenter2 "<<thetaofcenter2<<std::endl;
    Helper.SetDeltaTheta(deltatheta_endcap[i]);
    Helper.SetThetaOfCenter(thetaofcenter);
    Helper.SetDeltaTheta2(deltatheta_endcap[i+1]);
    Helper.SetThetaOfCenter2(thetaofcenter2);
    Helper.CalBasic();
    Helper.Getpt(pt);
    for(int i=0; i<8;i++){
      std::cout<<i<<" "<<pt[i].x()<<" "<<pt[i].y()<<" "<<pt[i].z()<<std::endl;
    }
    
    // Calculation of parameters to create Trap
    auto pDz = (pt[7]).z();
    auto pDy1     = ((pt[2]).y()-(pt[1]).y())*0.5;
    auto pDx1     = ((pt[1]).x()-(pt[0]).x())*0.5;
    auto pDx2     = ((pt[3]).x()-(pt[2]).x())*0.5;
    auto fTalpha1 = ((pt[2]).x()+(pt[3]).x()-(pt[1]).x()-(pt[0]).x())*0.25/pDy1;
    auto pAlp1 = std::atan(fTalpha1);
    auto pDy2     = ((pt[6]).y()-(pt[5]).y())*0.5;
    auto pDx3     = ((pt[5]).x()-(pt[4]).x())*0.5;
    auto pDx4     = ((pt[7]).x()-(pt[6]).x())*0.5;
    auto fTalpha2 = ((pt[6]).x()+(pt[7]).x()-(pt[5]).x()-(pt[4]).x())*0.25/pDy2;
    auto pAlp2 = std::atan(fTalpha2);
    auto fThetaCphi = ((pt[4]).x()+pDy2*fTalpha2+pDx3)/pDz;
    std::cout<<pt[4].x()<<" "<<pDy2<<" "<<fTalpha2<<" "<<pDx3<<" "<<pDz<<std::endl;
    auto fThetaSphi = ((pt[4]).y()+pDy2)/pDz;
    std::cout<<"theta c phi "<<fThetaCphi<<" "<<fThetaSphi<<std::endl;
    double pPhi = 0.;
    double pTheta = 0.;
    if(fThetaSphi==0. && fThetaCphi==0.){}
    else{
         pPhi = std::atan(fThetaSphi/fThetaCphi);
         pTheta = std::atan(fThetaCphi/std::cos(pPhi));
    }
    
    //sprintf(name,"tower%d",NbOfBarrel+i+1); 
    //std::cout<<name<<std::endl;

    Trap tower( "phiER", pDz, pTheta, pPhi, pDy1, pDx1, pDx2, pAlp1, pDy2, pDx3, pDx4, pAlp2);
    std::cout<<"tower dim "<<tower.dZ()<<" "<<tower.theta()<<" "<<tower.phi()<<std::endl;
    Volume towerLog("towerLog", tower, description.material(x_tank.attr<std::string>(_U(material))));
    towerLog.setVisAttributes(description, x_tower.visStr());

    //PlacedVolume towerPlaced = tank_vol.placeVolume(towerLog);
    //towerPlaced.addPhysVolID("towerPlaced",i);

    //G4RotationMatrix* rm = new G4RotationMatrix();
    //rm->rotateX(thetaofcenter);
    //RotationZYX tower_rot(0.,0.,-thetaofcenter); //phi theta psi
    //std::cout<<"around Z "<<tower_rot.Phi()<<" around X "<<tower_rot.Psi()<<" around Y "<<tower_rot.Theta()<<std::endl;

    RotationX rotX(-thetaofcenter);
    RotationY rotY(0.);
    RotationZ rotZ(0.);
    //G4ThreeVector c = dimE->GetOrigin(0);
    dd4hep::rec::Vector3D c = Helper.GetOrigin(0);
    std::cout<<"origin x "<<c.x()<<" "<<" y "<<c.y()<<" z "<<c.z()<<std::endl;
    //G4ThreeVector c_new(-c.getY(),c.getZ(),c.getX()-(innerR+0.5*length));
    dd4hep::rec::Vector3D c_new(-c.y(),c.z(),c.x()-(innerR+0.5*length));
    //if(i<35) new G4PVPlacement(rm,c_new,towerLV,name,phiERLog,false,NbOfBarrel+i+1,false);
    //if(i<35) {
    //    Transform3D tower_trnsform(rotX*rotY*rotZ, Position(c_new.x(),c_new.y(),c_new.z())); 
    //    phiERLog.placeVolume(towerLog,i,tower_trnsform);
    //}
    //if(i<35) tank_vol.placeVolume(towerLog,i,tower_rot);
    if(i==0) tank_vol.placeVolume(towerLog);
    //if(i<35) tank_vol.placeVolume(towerLog,i,rotX*rotY*rotZ);









    // Capillary placement inside tower(s)
    //
    // offset to be applied to capillary tubes
    // with length reduced by the intersection
    // with the tower's planes.
    // It is needed for the intersections with
    // the lateral sides for which the intersecting
    // point is not exactly the left or right edge of the tube.
    const double FiberLengthOffset{0.1*mm}; 
    const double tubeDiameter = tubeRadius*2.;
    const double y_pitch = tubeDiameter*sqrt(3)/2.; //y-distance from center to center
    //Fill a tower along y
    double y_backplane = pt[4].y();
    double x_backplane = pt[4].x();
    double x_start = 0;
    for(std::size_t k=14; k<15; k++){
      x_start = x_backplane;
      double y_tube = 0.;
      double delta_x = 0.;
      if(k==0) y_tube = y_backplane + tubeRadius + 1.0*mm; //add 1 mm otherwise the tube immediately intersect with a plane
      else {
        y_tube = y_backplane + tubeRadius + 1.0*mm + k*2.*y_pitch;

        // Adjust starting x given the opening angle of the Trap back face
        double hypotenuse = sqrt(pow((-1*y_backplane)*2,2) + pow(pt[6].x()-pt[4].x(),2));
        double angle = acos(((-1*y_backplane)*2) / hypotenuse);
        delta_x = ((-1.*y_backplane) - (-1.*y_tube)) * tan(angle);
        std::cout<<k<<" delta x "<<delta_x<<std::endl;
      }
      if(delta_x > tubeDiameter) {
          // Calculate how many fibers should be inserted in x
          // Caveat: casting to int does not round takes only the 
          // integer part (i.e. 6.9 -> 6)
          auto newFibersNo = static_cast<int>(delta_x/tubeDiameter);
          x_start = x_start - newFibersNo*tubeDiameter;
          std::cout<<k<<" longerrrrrrrrrrr delta_x "<<delta_x/cm<<" cm "<<newFibersNo<<" "<<newFibersNo*tubeDiameter<<std::endl;
          std::cout<<"x_start "<<x_start<<std::endl;
      }

      //Fill a tower row along x
      for(std::size_t j=0;j<5;j++){
         auto x_tube = x_start  + tubeRadius + j*(tubeDiameter); 
         std::cout<<"x_tube "<<x_tube<<" x_backplane "<<x_backplane<<" tubeRadius "<<tubeRadius<<std::endl;
         std::cout<<"y_tube "<<y_tube<<" y_backplane "<<y_backplane<<" tubeRadius "<<tubeRadius<<std::endl;
         Vector3D capillaryPos(x_tube, y_tube, length/2.);
         auto capillaryLength = GetFiberLength(pt,capillaryPos);
         //std::cout<<" my x y "<<x_tube<<" "<<y_tube<<" length "<<capillaryLength<<std::endl;
         if(capillaryLength == length){
           PlacedVolume capillaryPlaced = towerLog.placeVolume(capillaryLog, 1000*k+j, Position(x_tube, y_tube, 0.));
         }
         else{
           std::cout<<"short fiber"<<std::endl;
           std::cout<<x_tube<<" y_tube "<<y_tube<<" k "<<k<<std::endl;
           std::cout<<"length "<<capillaryLength<<std::endl;
           Tube capillaryShort(0.*mm, tubeRadius, (capillaryLength-FiberLengthOffset)/2., 2*M_PI); // reduced capillary length
                                                                                                   // by a fixed offset
           Volume capillaryShortLog("capillaryShortLog",capillaryShort,description.material(x_tank.attr<std::string>(_U(material))));
           capillaryShortLog.setVisAttributes(description, x_capillary.visStr());
           PlacedVolume capillaryShortPlaced = towerLog.placeVolume(capillaryShortLog, 1000*k+j, Position(x_tube, y_tube, length/2.-capillaryLength/2.+FiberLengthOffset/2.));
         }
         if((-1.*x_backplane)+delta_x-x_tube < (tubeDiameter+tubeRadius)) break;
      } // end x loop
 
      // y_backplane is equal up and down so I can keep the same for exiting loop
      if((-1.*y_backplane)-y_tube < (2.*y_pitch+tubeRadius)) break;
    } // end y loop




    // Update parameters
    Helper.Getpt(pt);
    fulltheta = fulltheta-deltatheta_endcap[i];
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
