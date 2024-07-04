#ifndef DREndcapTube_H
#define DREndcapTube_H

#include "DDRec/Vector3D.h" // Used for dd4hep::rec::Vector3D

using namespace dd4hep;

class DREndcapTubeHelper {

  public:
    // Constructor and de-constructor
    DREndcapTubeHelper() = default;
    ~DREndcapTubeHelper() = default;

  private:
    // Class fields
    double fPhiZRot;
    bool fRbool;
    bool fcalbasicbool;
    double finnerR;
    double ftower_height;
    double fnumzrot;
    double fdeltatheta;
    double fthetaofcenter;
    double finnerR_new;
    double fTrns_Length;
    dd4hep::rec::Vector3D fTrns_Vector;
    dd4hep::rec::Vector3D fV1;
    dd4hep::rec::Vector3D fV2;
    dd4hep::rec::Vector3D fV3;
    dd4hep::rec::Vector3D fV4;
    double fPMTT;
    double fthetaofcenter2;
    double fdeltatheta2;
    double finnerR_new2;
    double Ratio;
    double Ratio2;
 
  public: 
    // Methods to initialize class parameters
    void Rbool(bool Rbool) { fRbool = Rbool; } // set to True if building right endcap
    void SetInnerR(double innerR) { finnerR = innerR; }
    void SetTower_height(double tower_height) { ftower_height = tower_height; }
    void SetNumZRot(int num) { fnumzrot = num; fPhiZRot = 2*M_PI/(double)num; }
    void SetDeltaTheta(double theta) { fdeltatheta = theta; }
    void SetThetaOfCenter(double theta) { fthetaofcenter = theta; }
    void SetDeltaTheta2(double theta) { fdeltatheta2 = theta; }
    void SetThetaOfCenter2(double theta) { fthetaofcenter2 = theta; }
    void SetPMTT(double PMTT) { fPMTT = PMTT; }

    void CalBasic();
    /*
    // Methods to get class parameters
    G4double GetInnerR_new();
    G4double GetTrns_Length();
    G4ThreeVector GetTrns_Vector();
    G4ThreeVector GetV1();
    G4ThreeVector GetV2();
    G4ThreeVector GetV3();
    G4ThreeVector GetV4();
    */
    dd4hep::rec::Vector3D GetOrigin(int i);
    void Getpt(dd4hep::rec::Vector3D *pt);
		
    /*G4ThreeVector GetOrigin_PMTG(G4int i);
    void Getpt_PMTG(G4ThreeVector *pt);

    void Getpt_PMTCath(G4ThreeVector *pt);

    G4RotationMatrix* GetRM(G4int i);
    */
};

inline void DREndcapTubeHelper::CalBasic(){

    fcalbasicbool = 1;
    
    // distance from center to front face of first tower
    // radius is 2500 mm which is same as distance from center to front face of very
    // last tower.
    // To get distance to first tower I divide this radius by cosine of first tower.
    // To understand impact of sin*tan part
    finnerR_new = finnerR/(cos(fthetaofcenter)-sin(fthetaofcenter)*tan(fdeltatheta/2.));
    // distance from center to front face of second tower (same as above)
    finnerR_new2 = finnerR/(cos(fthetaofcenter2)-sin(fthetaofcenter2)*tan(fdeltatheta2/2.));
    
    // Half size at front face of first tower given by effective radius times tan of half delta
    double innerSide_half = finnerR_new*tan(fdeltatheta/2.);
    // Half size at back face
    double outerSide_half = (finnerR_new+ftower_height)*tan(fdeltatheta/2.);	
    
    // Half size at front face of second tower (same as above)
    double innerSide_half2 = finnerR_new2*tan(fdeltatheta2/2.);

    // Distance from origin to center of first tower
    fTrns_Length = ftower_height/2.+finnerR_new;
	
    // Vector from origin to center of first tower
    // Remember towers are placed starting from the one laying on the x-axis 
    fTrns_Vector = dd4hep::rec::Vector3D(cos(fthetaofcenter)*fTrns_Length,0,sin(fthetaofcenter)*fTrns_Length);
    
    double dx1=finnerR; // inner radius hardcoded in detector construction
    double dxi=sin(fthetaofcenter)*finnerR_new+innerSide_half*cos(fthetaofcenter);
    double dxi2=sin(fthetaofcenter2)*finnerR_new2+innerSide_half2*cos(fthetaofcenter2);
    Ratio=dxi/dx1;
    Ratio2=dxi2/dx1;
    std::cout<<"finnerR_new "<<finnerR_new/dd4hep::mm<<" innerSide_half "<<innerSide_half/dd4hep::mm<<" Ratio: "<<Ratio<<std::endl;

    // This vector are distributed over the plane of the tower closest to the barrel
    // The difference between fV1 and fV2 (fV3 and fV4) is the correction of finnerR_new+ftower_height that meake the vector start at the begin and end of the tower surface
    fV1 = dd4hep::rec::Vector3D(
			(Ratio2)*(cos(fthetaofcenter)*finnerR_new+sin(fthetaofcenter)*finnerR_new*tan(fdeltatheta/2.)),
			0,
			sin(fthetaofcenter)*finnerR_new-sin(fthetaofcenter)*finnerR_new*tan(fdeltatheta/2.)
			);
        
    std::cout<<fV1.x()/mm<<" "<<fV1.y()/mm<<" "<<fV1.z()/mm<<" fV1"<<std::endl;        
    fV2 = dd4hep::rec::Vector3D(
			(Ratio2)*(cos(fthetaofcenter)*(finnerR_new+ftower_height)+sin(fthetaofcenter)*(finnerR_new+ftower_height)*tan(fdeltatheta/2.)),
			0,
			sin(fthetaofcenter)*(finnerR_new+ftower_height)-sin(fthetaofcenter)*(finnerR_new+ftower_height)*tan(fdeltatheta/2.)
			);

    std::cout<<fV2.x()/mm<<" "<<fV2.y()/mm<<" "<<fV2.z()/mm<<" fV2"<<std::endl;        
    fV3 = dd4hep::rec::Vector3D(
			(Ratio)*(cos(fthetaofcenter)*finnerR_new-sin(fthetaofcenter)*finnerR_new*tan(fdeltatheta/2.)),
			0,
			sin(fthetaofcenter)*finnerR_new+sin(fthetaofcenter)*finnerR_new*tan(fdeltatheta/2.)
			);
	
    std::cout<<fV3.x()/mm<<" "<<fV3.y()/mm<<" "<<fV3.z()/mm<<" fV3"<<std::endl;        
    fV4 = dd4hep::rec::Vector3D(
			(Ratio)*(cos(fthetaofcenter)*(finnerR_new+ftower_height)-sin(fthetaofcenter)*(finnerR_new+ftower_height)*tan(fdeltatheta/2.)),
			0,
			sin(fthetaofcenter)*(finnerR_new+ftower_height)+sin(fthetaofcenter)*(finnerR_new+ftower_height)*tan(fdeltatheta/2.)
			);
    std::cout<<fV4.x()/mm<<" "<<fV4.y()/mm<<" "<<fV4.z()/mm<<" fV4"<<std::endl;        

}

inline void DREndcapTubeHelper::Getpt(dd4hep::rec::Vector3D *pt) {
	double innerSide_half = finnerR_new*tan(fdeltatheta/2.);
	double outerSide_half= (finnerR_new+ftower_height)*tan(fdeltatheta/2.);

	if(fRbool == 1){
		pt[0]=dd4hep::rec::Vector3D(-(fV1.x()*tan(fPhiZRot/2.)),-innerSide_half,-ftower_height/2.);
		pt[1]=dd4hep::rec::Vector3D((fV1.x()*tan(fPhiZRot/2.)),-innerSide_half,-ftower_height/2.);
		pt[2]=dd4hep::rec::Vector3D(-(fV3.x()*tan(fPhiZRot/2.)),innerSide_half,-ftower_height/2.);
		pt[3]=dd4hep::rec::Vector3D((fV3.x()*tan(fPhiZRot/2.)),innerSide_half,-ftower_height/2.);
		pt[4]=dd4hep::rec::Vector3D(-(fV2.x()*tan(fPhiZRot/2.)),-outerSide_half,ftower_height/2.);
		pt[5]=dd4hep::rec::Vector3D((fV2.x()*tan(fPhiZRot/2.)),-outerSide_half,ftower_height/2.);
		pt[6]=dd4hep::rec::Vector3D(-(fV4.x()*tan(fPhiZRot/2.)),outerSide_half,ftower_height/2.);
		pt[7]=dd4hep::rec::Vector3D((fV4.x()*tan(fPhiZRot/2.)),outerSide_half,ftower_height/2.);
	}

	else{
		pt[0]=dd4hep::rec::Vector3D(-(fV1.x()*tan(fPhiZRot/2.)),-innerSide_half,-ftower_height/2.);
		pt[1]=dd4hep::rec::Vector3D((fV1.x()*tan(fPhiZRot/2.)),-innerSide_half,-ftower_height/2.);
		pt[2]=dd4hep::rec::Vector3D(-(fV3.x()*tan(fPhiZRot/2.)),innerSide_half,-ftower_height/2.);
		pt[3]=dd4hep::rec::Vector3D((fV3.x()*tan(fPhiZRot/2.)),innerSide_half,-ftower_height/2.);
		pt[4]=dd4hep::rec::Vector3D(-(fV2.x()*tan(fPhiZRot/2.)),-outerSide_half,ftower_height/2.);
		pt[5]=dd4hep::rec::Vector3D((fV2.x()*tan(fPhiZRot/2.)),-outerSide_half,ftower_height/2.);
		pt[6]=dd4hep::rec::Vector3D(-(fV4.x()*tan(fPhiZRot/2.)),outerSide_half,ftower_height/2.);
		pt[7]=dd4hep::rec::Vector3D((fV4.x()*tan(fPhiZRot/2.)),outerSide_half,ftower_height/2.);
	}
	
        for(std::size_t i=0;i<8;i++){
            std::cout<<"pt at i "<<pt[i].x()/mm<<" "<<pt[i].y()/mm<<" "<<pt[i].z()/mm<<std::endl;
        }
	//cout<<"ENDCAP Y_int = "<<innerSide_half*2<<" Y_out = "<<outerSide_half*2<<" X_inn = "<<pt[3](0)-pt[2](0)<<" X_out = "<<pt[7](0)-pt[6](0)<<std::endl;
//cout<<std::endl;
}

inline dd4hep::rec::Vector3D DREndcapTubeHelper::GetOrigin(int i){

	if(fcalbasicbool==0) 
	{std::cout<<"fcalbasicbool = 0"<<std::endl; 
		return dd4hep::rec::Vector3D();
	}
	else{ 
		if(fRbool==1){
			double x,y,z;
			x=cos(i*fPhiZRot)*fTrns_Vector.x();
			y=sin(i*fPhiZRot)*fTrns_Vector.x();
			z=fTrns_Vector.z();

			return dd4hep::rec::Vector3D(x,y,z);
		}

		else
		{
			double x,y,z;
			x=cos(i*fPhiZRot)*fTrns_Vector.x();
			y=sin(i*fPhiZRot)*fTrns_Vector.x();
			z=-fTrns_Vector.z();

			return dd4hep::rec::Vector3D(x,y,z);
		}
	}
}

#endif // DREndcapTube_H
