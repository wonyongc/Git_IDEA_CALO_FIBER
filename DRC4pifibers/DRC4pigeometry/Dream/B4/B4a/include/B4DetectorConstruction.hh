
//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: B4DetectorConstruction.hh 75215 2013-10-29 16:07:06Z gcosmo $
//
/// \file B4DetectorConstruction.hh
/// \brief Definition of the B4DetectorConstruction class

#ifndef B4DetectorConstruction_h
#define B4DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4RotationMatrix.hh"
#include "G4FieldManager.hh"
#include <string.h>
#include <vector>
#include "dimensionB.hh"
#include "dimensionE.hh"
#include "G4Trap.hh"
#include "G4Tubs.hh"
#include "B4MyMaterials.hh"

#include "G4VSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Material.hh"
#include "G4VSensitiveDetector.hh"
#include "G4IntersectionSolid.hh"
#include "G4VPVParameterisation.hh"

class B4MagneticField;
class G4VPhysicalVolume;
class G4Material;
class G4VSensitiveDetector;
class G4VisAttributes;
class G4GenericMessenger;

class G4VPhysicalVolume;
class G4GlobalMagFieldMessenger;

class B4DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    B4DetectorConstruction();
    virtual ~B4DetectorConstruction();
    
public:
    virtual G4VPhysicalVolume* Construct();
    void fiberBR(G4int i,G4double deltatheta_);
    void fiberBL(G4int i,G4double deltatheta_);
    void fiberER(G4int i,G4double deltatheta_);
    void fiberEL(G4int i,G4double deltatheta_);
    
    void initializeMaterials();
    void placeSCEPCal_Barrel();
    
    
    virtual void ConstructSDandField();
    
    // get methods
    //
    const G4VPhysicalVolume* GetmodulePV() const;
    
private:
    // methods
    //
    
    G4GenericMessenger* fMessenger;
    
    static G4ThreadLocal B4MagneticField* fMagneticField;
    static G4ThreadLocal G4FieldManager* fFieldMgr;
    G4LogicalVolume* fMagneticLogical;
    
    void DefineMaterials();
    G4VPhysicalVolume* DefineVolumes();
    
    G4bool checkOverlaps;
    
    // data members
    //
    static G4ThreadLocal G4GlobalMagFieldMessenger*  fMagFieldMessenger;
    // magnetic field messenger
    
    G4VPhysicalVolume*   modulePV; // the module physical volume
    
    G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps
    
    G4VisAttributes* visAttrC;
    G4VisAttributes* visAttrS;
    std::vector<G4VisAttributes*> fVisAttributes;
    G4VisAttributes* visAttr;
    
    G4double innerR;
    G4double tower_height;
    
    G4int NbOfBarrel;
    G4int NbOfEndcap;
    G4int NbOfZRot;
    
    G4double theta_unit;
    G4double phi_unit;
    
    G4double deltatheta;
    G4double thetaofcenter;
    G4double fulltheta;
    G4double lastdeltatheta;
    G4ThreeVector pt[8];
    
    
    G4double SCEP_innerR ;
    G4double SCEP_xtal_L ;
    G4int SCEP_NbOfBarrel;
    G4int SCEP_NbOfEndcap;
    G4int SCEP_NbOfZRot  ;
    G4double SCEP_theta_unit;
    G4double SCEP_phi_unit;
    G4double FR_X0ratio;
    
    G4double SCEP_deltatheta;
    G4double SCEP_thetaofcenter;
    G4double SCEP_thetaofcenter2;
    G4double SCEP_fulltheta;
    G4double SCEP_lastdeltatheta;
    G4ThreeVector SCEP_pt[8];
    
    
    G4double solenoid_L ;
    G4double solenoid_IR;
    G4double solenoid_OR;
    
    G4double PMTT;
    
    dimensionB* dimB;
    dimensionE* dimE;
    
    char name[80];
    G4Trap* tower;
    G4Trap* pmtg;
    G4Trap* pmtcath;
    G4Tubs* fiber;
    G4Tubs* fiber_S;
    G4Tubs* fiber_C;
    G4Tubs* fiberS;
    G4Tubs* fiberC;
    G4VSolid* intersect;
    G4VSolid* intersect_;
    
    G4LogicalVolume* fiberCladCLog;
    G4LogicalVolume* fiberCladSLog;
    
    G4LogicalVolume* towerLogicalBR[40];
    G4LogicalVolume* towerLogicalBL[40];
    //     G4LogicalVolume* towerLogicalER[47];
    //     G4LogicalVolume* towerLogicalEL[47];
    G4LogicalVolume* towerLogicalER[40];
    G4LogicalVolume* towerLogicalEL[40];
    
    G4LogicalVolume* PMTGLogicalBR[40];
    G4LogicalVolume* PMTGLogicalBL[40];
    //     G4LogicalVolume* PMTGLogicalER[47];
    //     G4LogicalVolume* PMTGLogicalEL[47];
    G4LogicalVolume* PMTGLogicalER[40];
    G4LogicalVolume* PMTGLogicalEL[40];
    
    G4LogicalVolume* PMTCathLogicalBR[40];
    G4LogicalVolume* PMTCathLogicalBL[40];
    //     G4LogicalVolume* PMTCathLogicalER[47];
    //     G4LogicalVolume* PMTCathLogicalEL[47];
    G4LogicalVolume* PMTCathLogicalER[40];
    G4LogicalVolume* PMTCathLogicalEL[40];
    
    vector<G4LogicalVolume*> fiberLogical_BR[40];
    vector<G4LogicalVolume*> fiberLogical_BR_[40];
    vector<G4LogicalVolume*> fiberLogical_BL[40];
    vector<G4LogicalVolume*> fiberLogical_BL_[40];
    //     vector<G4LogicalVolume*> fiberLogical_ER[47];
    //     vector<G4LogicalVolume*> fiberLogical_ER_[47];
    //     vector<G4LogicalVolume*> fiberLogical_EL[47];
    //     vector<G4LogicalVolume*> fiberLogical_EL_[47];
    vector<G4LogicalVolume*> fiberLogical_ER[40];
    vector<G4LogicalVolume*> fiberLogical_ER_[40];
    vector<G4LogicalVolume*> fiberLogical_EL[40];
    vector<G4LogicalVolume*> fiberLogical_EL_[40];
    
    map<int,G4LogicalVolume*> fiberCLog;
    map<int,G4LogicalVolume*> fiberSLog;
    
    G4VSensitiveDetector* PMTSDBR[40];
    G4VSensitiveDetector* PMTSDBL[40];
    //     G4VSensitiveDetector* PMTSDER[47];
    //     G4VSensitiveDetector* PMTSDEL[47];
    G4VSensitiveDetector* PMTSDER[40];
    G4VSensitiveDetector* PMTSDEL[40];
    
    G4double clad_C_rMin;
    G4double clad_C_rMax;
    G4double clad_C_Dz  ;
    G4double clad_C_Sphi;
    G4double clad_C_Dphi;
    
    G4double core_C_rMin;
    G4double core_C_rMax;
    G4double core_C_Dz  ;
    G4double core_C_Sphi;
    G4double core_C_Dphi;
    
    G4double clad_S_rMin;
    G4double clad_S_rMax;
    G4double clad_S_Dz  ;
    G4double clad_S_Sphi;
    G4double clad_S_Dphi;
    
    G4double core_S_rMin;
    G4double core_S_rMax;
    G4double core_S_Dz  ;
    G4double core_S_Sphi;
    G4double core_S_Dphi;
    
    //---Materials for Cerenkov fiber---
    G4Material *clad_C_Material;
    G4Material *core_C_Material;
    
    //---Materials for Scintillation fiber---
    G4Material *clad_S_Material;
    G4Material *core_S_Material;
    
    //--Material for PMT glass---
    G4Material *Glass_Material;
    
    //--- Material for PMT Photocathod ---
    G4Material *PMTPC_Material;
    
    
    //SCEPCal 
    G4Material* SCEPCalMaterial;
    
    G4Trap* crystal_BF;    
    G4Trap* crystal_BR;    
    G4LogicalVolume* crystalLogicalF_B;
    G4LogicalVolume* crystalLogicalR_B;    
    
            
    
    G4Cons *SCEP_EndcapRing_F;
    G4Cons *SCEP_EndcapRing_R;
    
    G4Cons* crystal_EF;
    G4Cons* crystal_ER;   
    G4LogicalVolume* SCEP_endcapRingLog_F;
    G4LogicalVolume* SCEP_endcapRingLog_R;
    
    
    G4Cons *SCEP_EndcapRing_F_L;
    G4Cons *SCEP_EndcapRing_R_L;
    
    G4Cons* crystal_EF_L;
    G4Cons* crystal_ER_L;   
    G4LogicalVolume* SCEP_endcapRingLog_F_L;
    G4LogicalVolume* SCEP_endcapRingLog_R_L;
    
    
    
    
};

// inline functions

inline const G4VPhysicalVolume* B4DetectorConstruction::GetmodulePV() const {
    return modulePV;
}



class BarrelPhiParameterisation : public G4VPVParameterisation
{
    G4double fInnerR;
    G4double fPhiUnit;
    G4double fHeight;
    G4double fOffset;
        
    public:
        
        BarrelPhiParameterisation(G4double inner_R, G4double phi_unit, G4double height, G4double offset);
        
        ~BarrelPhiParameterisation(){};
        
        void ComputeTransformation (const G4int copyNo, G4VPhysicalVolume* physVol) const;                
        
};


class EndcapThetaParameterisation : public G4VPVParameterisation
{
    G4double fInnerR;
    G4double fThetaUnit;
    G4double fHeight;
    G4double fThetaB;
        
    public:
        
        EndcapThetaParameterisation(G4double inner_R, G4double theta_unit, G4double height, G4double thetaB);
        
        ~EndcapThetaParameterisation(){};
        
        void ComputeTransformation (const G4int copyNo, G4VPhysicalVolume* physVol) const;                
        
        void ComputeDimensions (G4Cons& endcapRing, const G4int copyNo, const G4VPhysicalVolume* physVol) const;
        
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

