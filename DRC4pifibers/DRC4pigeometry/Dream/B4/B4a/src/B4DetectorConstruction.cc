
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
// $Id: B4DetectorConstruction.cc 87359 2014-12-01 16:04:27Z gcosmo $
//
/// \file B4DetectorConstruction.cc
/// \brief Implementation of the B4DetectorConstruction class

#include "B4DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"
#include "G4Sphere.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4OpBoundaryProcess.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"

#include "G4Trap.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Polycone.hh"
#include "G4SubtractionSolid.hh"
#include "G4Trd.hh"
#include "G4GenericTrap.hh"

#include "dimensionB.hh"
#include "dimensionE.hh"
#include "B4MyMaterials.hh"


#include <string>

#include "B4MagneticField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4GenericMessenger.hh"
#include "G4PVParameterised.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal B4MagneticField* B4DetectorConstruction::fMagneticField = 0;
G4ThreadLocal G4FieldManager* B4DetectorConstruction::fFieldMgr = 0;

B4DetectorConstruction::B4DetectorConstruction()
: G4VUserDetectorConstruction(),
modulePV(0),
fCheckOverlaps(true),
fMagneticLogical(nullptr)
{
}

B4DetectorConstruction::~B4DetectorConstruction()
{
}

G4VPhysicalVolume* B4DetectorConstruction::Construct()
{
    // Define materials
	DefineMaterials();
    // Define volumes
	return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4DetectorConstruction::DefineMaterials()
{
    // Copper material defined using NIST Manager
    // I use Cu as default absorber material but you can switch to lead
    G4NistManager* CunistManager = G4NistManager::Instance();
    CunistManager->FindOrBuildMaterial("G4_Cu");

    // Lead material defined using NIST Manager
    // G4NistManager* PbnistManager = G4NistManager::Instance();
    // PbnistManager->FindOrBuildMaterial("G4_Pb");

    // Polystyrene material defined using NIST Manager
    // I use this material for the core of plastic scintillating fibers
    // cannot find any G4_Polystyrene, I build it later
    //G4NistManager* PynistManager = G4NistManager::Instance();
    //PynistManager->FindOrBuildMaterial("G4_Polystyrene");

    // PMMA material, there's no default G4_PMMA, I build it (C502H8)
    G4String name_e, symbol;    // a=mass of a mole;
    G4double a, z;            // z=mean number of protons;
    
    // create elements
    a = 1.01*g/mole;
    G4Element* elH  = new G4Element(name_e="Hydrogen",symbol="H" , z= 1., a); //Hidrogen
    
    a = 12.01*g/mole;
    G4Element* elC  = new G4Element(name_e="Carbon"  ,symbol="C" , z= 6., a); //Carbon
    
    a = 16.00*g/mole;
    G4Element* elO  = new G4Element(name_e="Oxygen"  ,symbol="O" , z= 8., a); //Oxygen
    
    a = 28.09*g/mole;
    G4Element* elSi = new G4Element(name_e="Silicon", symbol="Si", z=14., a); //Silicon
    
    a = 18.9984*g/mole;
    G4Element* elF  = new G4Element("Fluorine",symbol="F" , z= 9., a); //Fluorine
    
    a = 63.546*g/mole;
    G4Element* elCu = new G4Element("Copper", symbol="Cu", z=29., a); //Copper
    
    a = 65.38*g/mole;
    G4Element* elZn = new G4Element("Zinc", symbol="Zn", z=30., a); //Zinc
    
    // create PMMA
    G4Material* PMMA = new G4Material("PMMA", 1.19*g/cm3, 3); //name, density and number of elements
    PMMA -> AddElement(elC, 5);
    PMMA -> AddElement(elO, 2);
    PMMA -> AddElement(elH, 8); //PMMA building complete
    
    // create Polystyrene (C5H5)
    G4Material* Polystyrene = new G4Material("Polystyrene", 1.05*g/cm3, 2);
    Polystyrene -> AddElement(elC, 8);
    Polystyrene -> AddElement(elH, 8); //Polystyrene building complete
    
    // create Fluorinated Polymer (C2F2)
    // I use it for the cladding of the Cherenkov fibers
    G4Material* fluorinatedPolymer =
    new G4Material("Fluorinated_Polymer", 1.43*g/cm3, 2);
    fluorinatedPolymer->AddElement(elC,2);
    fluorinatedPolymer->AddElement(elF,2);
   //fluorinatedPolymer->AddElement(H,2); //Fluorinated Polymer building complete
    
    // create Glass (SiO2)
    G4Material* Glass = new G4Material("Glass", 2.4*g/cm3, 2);
    Glass -> AddElement(elSi, 1);
    Glass -> AddElement(elO, 2); //Glass building complete
    
    // Vacuum material defined using NIST Manager
    G4NistManager* VanistManager = G4NistManager::Instance();
    VanistManager->FindOrBuildMaterial("G4_Galactic");
    
    // Silicon material defined using NIST Manager
    G4NistManager* SinistManager = G4NistManager::Instance();
    SinistManager->FindOrBuildMaterial("G4_Si");
    
    // create Cu260 (Brass)
    // I use it for the absorber of the real small beam tested module
    double density = 8.53*g/cm3;
    int ncomponentsbrass = 2;
    G4Material* Cu260 = new G4Material(name_e="Brass", density, ncomponentsbrass);
    Cu260->AddElement(elCu, 70*perCent);
    Cu260->AddElement(elZn, 30*perCent);
    
    // Air material defined using NIST Manager
    // You can use Air instead of vacuum
    G4NistManager* AirnistManager = G4NistManager::Instance();
    AirnistManager->FindOrBuildMaterial("G4_AIR");
    
    // Print materials
    // I don't want to print materials all the times,
    // if you want uncomment it
    //G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B4DetectorConstruction::DefineVolumes()
{
    
    // define colors
    G4Colour white  (1.00, 1.00, 1.00);  // white
    G4Colour gray   (0.50, 0.50, 0.50);  // gray
    G4Colour black  (0.00, 0.00, 0.00);  // black
    G4Colour red    (1.00, 0.00, 0.00);  // red
    G4Colour green  (0.00, 1.00, 0.00);  // green
    G4Colour blue   (0.00, 0.00, 1.00);  // blue
    G4Colour cyan   (0.00, 1.00, 1.00);  // cyan
    G4Colour air    (0.90, 0.94, 1.00);  // cyan
    G4Colour magenta(1.00, 0.00, 1.00);  // magenta 
    G4Colour yellow (1.00, 1.00, 0.00);  // yellow
    G4Colour brass  (0.80, 0.60, 0.40);  // brass
    G4Colour brown  (0.70, 0.40, 0.10);  // brown
  
    // Geometry parameters of the world, world is a box
    G4double worldX = 14*m;
    G4double worldY = 14*m;
    G4double worldZ = 14*m;
    
    // Get materials for vacuum, absorber, scintillating and cherenkov fibers, SiPM
    G4Material* defaultMaterial = G4Material::GetMaterial("G4_AIR"); // G4_AIR or G4_Galactic
    
    
    
    // Building the calorimeter
    // Build the world
    G4VSolid* worldS
    = new G4Box("World",                        // its name
                worldX/2, worldY/2, worldZ/2); // its size
    
    G4LogicalVolume* worldLV
    = new G4LogicalVolume(
                          worldS,           // its solid
                          defaultMaterial,  // its material (Galactic or Air)
                          "World");         // its name
    
    // Set the world as invisible
    worldLV->SetVisAttributes(G4VisAttributes::Invisible);
    
    G4VPhysicalVolume* worldPV
    = new G4PVPlacement(
                        0,                // no rotation
                        G4ThreeVector(),  // at (0,0,0)
                        worldLV,          // its logical volume
                        "World",          // its name
                        0,                // its mother  volume
                        false,            // no boolean operation
                        0,                // copy number
                        fCheckOverlaps);  // checking overlaps
    //absorber to calculate leakage
    G4VSolid* leakageabsorber
    = new G4Sphere("leakageabsorber",                        // its name
                6000., 6500., 0.*deg, 360.*deg, 0.*deg, 180.*deg); // its size
    
    G4LogicalVolume* leakageabsorberLV
    = new G4LogicalVolume(
                          leakageabsorber,           // its solid
                          defaultMaterial,  // its material (Galactic or Air)
                          "leakageabsorber");         // its name
    
    leakageabsorberLV->SetVisAttributes(G4VisAttributes::Invisible);   
    G4VPhysicalVolume* leakageabsorberPV
    = new G4PVPlacement(
                        0,                // no rotation
                        G4ThreeVector(),  // at (0,0,0)
                        leakageabsorberLV,          // its logical volume
                        "leakageabsorber",          // its name
                        worldLV,                // its mother  volume
                        false,            // no boolean operation
                        0,                // copy number
                        fCheckOverlaps);  // checking overlaps

    
    G4NistManager* nistManager = G4NistManager::Instance();
    
    G4String symbol;             //a=mass of a mole;
    G4double a, z, density;      //z=mean number of protons;
//     G4int iz, n;                 //iz=number of protons  in an isotope;
    // n=number of nucleons in an isotope;
    
    G4int ncomponents, natoms;
    G4double abundance, fractionmass;
    G4Material* cu  =new G4Material("Copper"   , z=29., a=63.546*g/mole, density=8.96*g/cm3);
    G4Element* H  = nistManager->FindOrBuildElement(1);
    G4Element* C  = nistManager->FindOrBuildElement(6);
    G4Element* N  = nistManager->FindOrBuildElement(7);
    G4Element* O  = nistManager->FindOrBuildElement(8);
    G4Element* F  = nistManager->FindOrBuildElement(9);
    G4Element* Si = nistManager->FindOrBuildElement(14);
    
    
    
    
    //---for Solenoid COIL
    G4Material* Fe = new G4Material("Iron",z=26., a=55.845*g/mole, density=7.874*g/cm3);
    
    //---for PSD
    G4Material* Pb = new G4Material("Lead",z=82., a=207.2*g/mole, density=11.35*g/cm3);
    
    //--- for PMT Cathod ---
    G4Material* Al = new G4Material("Aluminium", z=13., a=26.98*g/mole, density=2.700*g/cm3);
    
    //--- for PMT Glass ---
    G4Material* Glass = new G4Material("Glass", density=1.032*g/cm3,2);
    Glass->AddElement(C,91.533*perCent);
    Glass->AddElement(H,8.467*perCent);
    
    ///--- for scintillation fiber core ---
    G4Material* polystyrene =
    new G4Material("Polystyrene",density= 1.05*g/cm3, ncomponents=2);
    polystyrene->AddElement(C, natoms=8);
    polystyrene->AddElement(H, natoms=8);
    
    ///--- for cladding (scintillation fibers) ---
    G4Material* pmma_clad =
    new G4Material("PMMA_Clad",density= 1.19*g/cm3, ncomponents=3);
    pmma_clad->AddElement(C, natoms=5);
    pmma_clad->AddElement(H, natoms=8);
    pmma_clad->AddElement(O, natoms=2);
    
    ///--- for Cerenkov fiber core ---
    G4Material* pmma =
    new G4Material("PMMA",density= 1.19*g/cm3, ncomponents=3);
    pmma->AddElement(C, natoms=5);
    pmma->AddElement(H, natoms=8);
    pmma->AddElement(O, natoms=2);
    
    ///--- for cladding (Cerenkov fibers) ---
    G4Material* fluorinatedPolymer =
    new G4Material("Fluorinated_Polymer", density= 1.43*g/cm3, ncomponents=2);
    fluorinatedPolymer->AddElement(C,2);
    fluorinatedPolymer->AddElement(F,2);
    
    G4Material* Air = nistManager->FindOrBuildMaterial("G4_AIR",false);
    
    ///--- Material property tables for fiber materials ---
    G4MaterialPropertiesTable* mpAir;
    G4MaterialPropertiesTable* mpPS;
    G4MaterialPropertiesTable* mpPMMA;
    G4MaterialPropertiesTable* mpFS;
    G4MaterialPropertiesTable* mpGlass;
    G4MaterialPropertiesTable* mpPMTPC;
    
    //--- Generate and add material properties table ---
    G4double PhotonEnergy[] = {2.00*eV,2.03*eV,2.06*eV,2.09*eV,2.12*eV,
    	2.15*eV,2.18*eV,2.21*eV,2.24*eV,2.27*eV,
    	2.30*eV,2.33*eV,2.36*eV,2.39*eV,2.42*eV,
    	2.45*eV,2.48*eV,2.51*eV,2.54*eV,2.57*eV,
    	2.60*eV,2.63*eV,2.66*eV,2.69*eV,2.72*eV,
    	2.75*eV,2.78*eV,2.81*eV,2.84*eV,2.87*eV,
    	2.90*eV,2.93*eV,2.96*eV,2.99*eV,3.02*eV,
    	3.05*eV,3.08*eV,3.11*eV,3.14*eV,3.17*eV,
    	3.20*eV,3.23*eV,3.26*eV,3.29*eV,3.32*eV,
    	3.35*eV,3.38*eV,3.41*eV,3.44*eV,3.47*eV};

    	
    const G4int nEntries = sizeof(PhotonEnergy) / sizeof(G4double);
    //--- PMMA ---
    G4double RefractiveIndex_PMMA[nEntries] =
    {
            1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
            1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
            1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
            1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
            1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49
    };
    mpPMMA = new G4MaterialPropertiesTable();
    mpPMMA->AddProperty("RINDEX",PhotonEnergy,RefractiveIndex_PMMA,nEntries);
    pmma->SetMaterialPropertiesTable(mpPMMA);

    //--- Fluorinated Polymer (FS) ---
    G4double RefractiveIndex_FluorinatedPolymer[nEntries] =
    {
            1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42,
            1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42,
            1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42,
            1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42,
            1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42
    };
    mpFS = new G4MaterialPropertiesTable();
    mpFS->AddProperty("RINDEX",PhotonEnergy,RefractiveIndex_FluorinatedPolymer,nEntries);
    fluorinatedPolymer->SetMaterialPropertiesTable(mpFS);

    //Glass
    G4double RefractiveIndex_Glass[nEntries] =
    {  1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
            1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
            1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
            1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
            1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49
    };
    mpGlass = new G4MaterialPropertiesTable();
    mpGlass->AddProperty("RINDEX",PhotonEnergy,RefractiveIndex_Glass,nEntries);
    //mpGlass->AddProperty("ABSLENGTH",PhotonEnergy,Glass_AbsLength,nEntries);
    //Glass->SetMaterialPropertiesTable(mpGlass);
    /*if you want air refractive index
    	G4double RefractiveIndex_Air[nEntries] =
    	{
    		1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
    		1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
    		1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
    		1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
    		1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00
    	};

    	mpAir = new G4MaterialPropertiesTable();
    	mpAir->AddProperty("RINDEX", PhotonEnergy, RefractiveIndex_Air, nEntries);
    	Air->SetMaterialPropertiesTable(mpAir);*/
    //default materials of the World
    //---Materials for Cerenkov fiber---
    	clad_C_Material = fluorinatedPolymer;
    	core_C_Material = pmma;
    //---Materials for Scintillation fiber---
    	clad_S_Material = pmma_clad;
    	core_S_Material = polystyrene;
    //--Material for PMT glass---
    	Glass_Material = Glass;
    //--- Material for PMT Photocathod ---
    	PMTPC_Material = Al;
    //--- Photocathod property ---
    /* uncomment if you want photocayhod properties
    	G4double p_mppc[2] = {2.00*eV, 3.47*eV};
    	G4double refl_mppc[2] = {0.0, 0.0};
    G4double effi_mppc[2] = {0.11, 0.11}; // mimic Quantum Efficiency
    G4double photocath_ReR[] = {1.92, 1.92};
    G4double photocath_ImR[] = {1.69, 1.69};
    
    mpPMTPC = new G4MaterialPropertiesTable();
    mpPMTPC->AddProperty("REFLECTIVITY",p_mppc,refl_mppc,2);
    mpPMTPC->AddProperty("EFFICIENCY",p_mppc,effi_mppc,2);
    
    G4OpticalSurface* photocath_opsurf = new G4OpticalSurface("photocath_opsurf",glisur,polished,dielectric_metal);
    photocath_opsurf->SetMaterialPropertiesTable(mpPMTPC);*/
    
    bool placeFIBERS  = true;
    bool placeHCAL    = true;
    bool placeSCEPCAL = true;
    bool placeTiming  = true;
    
    bool displaySlice = false;    
    bool displayECALHalfSlice = false;
    bool displayTimingModule  = false;
    
    
    ////////////// Calorimeter parameters
    innerR = 2500; //inner radius /1800
    tower_height = 2000; //tower height 2500
    NbOfBarrel = 40; //(it was 52 before) number of towers in barrel right (left)
    NbOfEndcap = NbOfBarrel-1; //number of towers in endcap
    NbOfZRot = 252; //number of Z to round around the center
    //PMTT = 1*mm;
    PMTT = 0*mm;
    fulltheta = 0;

    
    //2*pi/number of tower to complete a rotation around the center
    phi_unit = 2*M_PI/(G4double)NbOfZRot;
    
    
    // Segmented Crystal Electromagnetic Section parameters (SCEPCal)
    SCEPCalMaterial     = MyMaterials::PWO();
    SCEP_TimingMaterial = MyMaterials::LSO();
//     SCEPCalMaterial = Fe;
    
    G4double crystal_size   = 1*cm; 
//     G4double crystal_size   = 1.15*cm; // to have <1Mch --> ~ 922 906 channels
    SCEP_innerR = 1800;
    SCEP_xtal_L = 200;
    FR_X0ratio  = 6./22.;
    
    SCEP_NbOfBarrel = SCEP_innerR/crystal_size;;               
    SCEP_NbOfEndcap = SCEP_NbOfBarrel-1; 
//     SCEP_NbOfEndcap = 40; 
    SCEP_NbOfZRot = SCEP_innerR*2*M_PI/crystal_size;                 
    
//     SCEP_NbOfBarrel = 100;
//     SCEP_NbOfEndcap = SCEP_NbOfBarrel-1; 
//     SCEP_NbOfZRot = 2*M_PI*100;                 
    
//     SCEP_NbOfBarrel = 50;
//     SCEP_NbOfEndcap = SCEP_NbOfBarrel-1;//   
//     SCEP_NbOfZRot = 2*M_PI*50;   
    
    SCEP_phi_unit = 2.*M_PI/(G4double)SCEP_NbOfZRot;
    SCEP_fulltheta = 0;
    
    
    //solenoid parameters
    solenoid_L  = 2.1*m;
    solenoid_IR = 2.11757*m;
    solenoid_OR = 2.5*m;
    
    
    //////////////
    
    //Parameters for fibers
    clad_C_rMin = 0.49*mm; //cladding cherenkov minimum radius
    clad_C_rMax = 0.50*mm; //cladding cherenkov max radius
    clad_C_Dz   = 2.5*m;   //cladding cherenkov lenght
    clad_C_Sphi = 0.;      //cladding cherenkov min rotation
    clad_C_Dphi = 2.*M_PI; //cladding chrenkov max rotation
    
    core_C_rMin = 0.*mm;
    core_C_rMax = 0.49*mm;
    core_C_Dz   = 2.5*m;
    core_C_Sphi = 0.;
    core_C_Dphi = 2.*M_PI;
    
    clad_S_rMin = 0.485*mm;
    clad_S_rMax = 0.50*mm;
    clad_S_Dz   = 2.5*m;
    clad_S_Sphi = 0.;
    clad_S_Dphi = 2.*M_PI;
    
    core_S_rMin = 0.*mm;
    core_S_rMax = 0.485*mm;
    core_S_Dz   = 2.5*m;
    core_S_Sphi = 0.;
    core_S_Dphi = 2.*M_PI;
    
    //Inizialise for Barrel R
    theta_unit=0; 
    deltatheta=0; 
    thetaofcenter=0;
    
    //creating fibers solids
    //G4cout << "r_clad= " << clad_C_rMax << " r_coreC=" << core_C_rMax << " r_coreS=" << core_S_rMax << G4endl;
    fiber  = new G4Tubs("Fiber", 0,clad_C_rMax,tower_height/2.,0*deg,360.*deg);// S is the same
    fiberC = new G4Tubs("fiberC",0,core_C_rMax,tower_height/2.,0*deg,360.*deg);
    fiberS = new G4Tubs("fiberS",0,core_S_rMax,tower_height/2.,0*deg,360.*deg);
    
    //vector for logical volumes of fibers
    //G4LogicalVolume* fiberCLog[2500];
    //G4LogicalVolume* fiberSLog[2500];
	
    ////////////// Tube with Local Magnetic field //////////////////////// 

    auto magneticSolid = new G4Tubs("magneticTubs",0.,2.5*m,(250-2*0.5612)*cm,0.,360.*deg);    
//     auto magneticSolid = new G4Tubs("magneticTubs",0,1.79*m,1.98*m,0.,360.*deg);
//     auto magneticSolid = new G4Tubs("magneticTubs",0.,solenoid_L,,0.,360.*deg);
    fMagneticLogical = new G4LogicalVolume(magneticSolid, Air, "magneticLogical");
    G4RotationMatrix* fieldRot = new G4RotationMatrix();
    new G4PVPlacement(fieldRot,G4ThreeVector(),fMagneticLogical, 
                      "magneticPhysical",worldLV,
                       false,0,checkOverlaps);
                      
  
    // set step limit in tube with magnetic field  
    //G4UserLimits* userLimits = new G4UserLimits(1*m);
    //fMagneticLogical->SetUserLimits(userLimits);
    G4VisAttributes* MagnetAttributes = new G4VisAttributes(G4Colour(0.9,0.9,0.9));   // LightGray
    MagnetAttributes->SetForceLineSegmentsPerCircle(100);
    MagnetAttributes->SetVisibility(false);
    fMagneticLogical->SetVisAttributes(MagnetAttributes);	
	
    ////////////// Tube COIL //////////////////////// 
    // I did it of 0.74 X0 IRON: X0 = 1.757 cm -> here 1X0 to take into account preshower
//     auto Solenoid = new G4Tubs("Solenoid",solenoid_L,solenoid_IR,solenoid_OR,0.,360.*deg);
//     auto Solenoid = new G4Tubs("Solenoid",solenoid_IR,solenoid_OR,solenoid_L, 0.,360.*deg);
//     auto Solenoid = new G4Tubs("Solenoid",2.1*m,2.11757*m,2.5*m,0.,360.*deg);
    G4double solenoid_X0 = 0.7*1.757*cm;
    auto Solenoid = new G4Tubs("Solenoid",2.1*m, 2.1*m+solenoid_X0, 2.5*m,0.,360.*deg);
    if (displaySlice) Solenoid = new G4Tubs("Solenoid",solenoid_IR, solenoid_OR, 2.5*m,90.*deg,1.*deg);
    G4LogicalVolume* SolenoidLV = new G4LogicalVolume(Solenoid, Fe, "SolenoidLV");
    G4RotationMatrix* SolenoidRot = new G4RotationMatrix();
    new G4PVPlacement(SolenoidRot,G4ThreeVector(),SolenoidLV,
                        "SolenoidPV",worldLV,
                        false,0,checkOverlaps);
    G4VisAttributes* SolenoidAttributes = new G4VisAttributes(G4Colour(1.,0.,0.));   // red
    SolenoidAttributes->SetForceLineSegmentsPerCircle(100);
    SolenoidAttributes->SetVisibility(false);        
    if (displaySlice) SolenoidAttributes->SetVisibility(true);

    SolenoidLV->SetVisAttributes(SolenoidAttributes);	
	
    ////////////// PSD //////////////////////// 
    // I did it of 1 X0 LEAD: X0 = 0.5612 cm
    auto PSDdx = new G4Tubs("PSDdx",250.8367*mm,solenoid_L,(0.5612/2.0)*cm,0.,360.*deg);
    G4LogicalVolume* PSDdxLV = new G4LogicalVolume(PSDdx, Pb, "PSDdxLV");
    G4RotationMatrix* PSDdxRot = new G4RotationMatrix();
    /*new G4PVPlacement(PSDdxRot,G4ThreeVector(0.,0.,(250-0.5612/2.0)*cm),PSDdxLV,
                    "PSDdxPV",worldLV,
                    false,0,checkOverlaps);*/
    
    G4VisAttributes* PSDAttributes = new G4VisAttributes(G4Colour(0.,1.,0.));   // green
    PSDdxLV->SetVisAttributes(PSDAttributes);	
    // I did it of 1 X0 LEAD: X0 = 0.5612 cm
    auto PSDsx = new G4Tubs("PSDsx",250.8367*mm,solenoid_L,(0.5612/2.0)*cm,0.,360.*deg);
    G4LogicalVolume* PSDsxLV = new G4LogicalVolume(PSDsx, Pb, "PSDsxLV");
    G4RotationMatrix* PSDsxRot = new G4RotationMatrix();
    /*new G4PVPlacement(PSDsxRot,G4ThreeVector(0.,0.,(-250+0.5612/2)*cm),PSDsxLV,
                    "PSDsxPV",worldLV,
                    false,0,checkOverlaps);*/
        
    PSDsxLV->SetVisAttributes(PSDAttributes);	
    //////////////////////////////////////////////////////////
	
    // Prepare for logical volume of fiber tower_height=2000
    for(int length=1;length<=tower_height;length++){ //from 1 to 20000
        double half=0.5*length; //half to build objects with proper dimensions
//         char name[80];
        sprintf(name,"Fiber%d",length);
        fiber = new G4Tubs(name,0,clad_C_rMax,half,0*deg,360.*deg); //creating fibers G4Tubs
        sprintf(name,"fiberC%d",length);
        fiberC = new G4Tubs(name,0,core_C_rMax,half,0*deg,360.*deg);
        sprintf(name,"fiberS%d",length);
        fiberS = new G4Tubs(name,0,core_S_rMax,half,0*deg,360.*deg);
        visAttrC = new G4VisAttributes(G4Colour(0.,0.,1.0));
        visAttrC->SetVisibility(true);
        visAttrC->SetDaughtersInvisible(true);
        visAttrC->SetForceWireframe(true);
        visAttrC->SetForceSolid(true);
        visAttrS = new G4VisAttributes(G4Colour(1.,0.,0.));
        visAttrS->SetVisibility(true);
        visAttrS->SetDaughtersInvisible(true);
        visAttrS->SetForceWireframe(true);
        visAttrS->SetForceSolid(true);
        
        fiberCLog[length] = new G4LogicalVolume(fiber,clad_C_Material,"fiberCladC");
        fiberSLog[length] = new G4LogicalVolume(fiber,clad_S_Material,"fiberCladS");
        fiberCLog[length]->SetVisAttributes(visAttrC);
        fiberSLog[length]->SetVisAttributes(visAttrS);
        G4LogicalVolume* fiberCoreCLog = new G4LogicalVolume(fiberC,core_C_Material,"fiberCoreC");
        G4LogicalVolume* fiberCoreSLog = new G4LogicalVolume(fiberS,core_S_Material,"fiberCoreS");
        new G4PVPlacement(0,G4ThreeVector(0,0,0),fiberCoreCLog,"fiberCoreCherePhys",fiberCLog[length],false,0);
        new G4PVPlacement(0,G4ThreeVector(0,0,0),fiberCoreSLog,"fiberCoreScintPhys",fiberSLog[length],false,0);
        /*if(sd){
         fiberCoreCLog->SetSensitiveDetector(sd);
         fiberCoreSLog->SetSensitiveDetector(sd);
         }*/
    }
    //Final logical volumes of fibers
    //fiberCladCLog = fiberCLog[2500];
    //fiberCladSLog = fiberSLog[2500];
    
    //Counter on number of volumes
    G4int volnum=0;
    // vector length has to be the same of NbOfBarrel
    G4double deltatheta_barrel[40] = {0};
    for(int i=0;i<NbOfBarrel;i++) deltatheta_barrel [i] = M_PI/4/(NbOfBarrel);
    G4double deltatheta_endcap[40] = {0};
    for(int i=0;i<NbOfEndcap+1;i++) deltatheta_endcap [i] = M_PI/4/(NbOfEndcap+1);
    double thetaB = 0;
    for(int i=0;i<NbOfBarrel;i++) thetaB += deltatheta_barrel[i];
    double thetaE = 0;
    for(int i=0;i<NbOfEndcap+1;i++) thetaE += deltatheta_endcap[i];
    
    double length = tower_height; //(was tower height + 100)length of physical volumes
    double innerR_Endcap = innerR;
    
    //BARREL
    G4Trd* phiBarrel = new G4Trd("phiBarrel",(innerR)*tan(0.5*phi_unit),(innerR+length)*tan(0.5*phi_unit),(innerR)*tan(thetaB),(innerR+length)*tan(thetaB),0.5*length);
    G4LogicalVolume* phiBLog = new G4LogicalVolume(phiBarrel,Air,"phiBLog");
    
    // ENDCAP R
    // I use the G4Generictrap class, so I define the 8 points needed 
    vector<G4TwoVector> vertices; 
    vertices.push_back( G4TwoVector(0,0));
    vertices.push_back( G4TwoVector(0,0));
    vertices.push_back( G4TwoVector(-innerR*tan(0.5*phi_unit),innerR));
    vertices.push_back( G4TwoVector(innerR*tan(0.5*phi_unit),innerR));
    vertices.push_back( G4TwoVector(0,0));
    vertices.push_back( G4TwoVector(0,0));
    vertices.push_back( G4TwoVector(-(innerR+tower_height)*tan(0.5*phi_unit),innerR+tower_height));
    vertices.push_back( G4TwoVector((innerR+tower_height)*tan(0.5*phi_unit),innerR+tower_height));
    
    G4GenericTrap* phiER = new G4GenericTrap( "phiER", tower_height/2., vertices);
    G4LogicalVolume* phiERLog = new G4LogicalVolume(phiER,Air,"phiERLog");
    // ENDCAP L
    // I use the G4Generictrap class, so I define the 8 points needed 
    vector<G4TwoVector> vertices2; 
    vertices2.push_back( G4TwoVector(0,0));
    vertices2.push_back( G4TwoVector(0,0));
    vertices2.push_back( G4TwoVector(-innerR*tan(0.5*phi_unit),innerR));
    vertices2.push_back( G4TwoVector(innerR*tan(0.5*phi_unit),innerR));
    vertices2.push_back( G4TwoVector(0,0));
    vertices2.push_back( G4TwoVector(0,0));
    vertices2.push_back( G4TwoVector(-(innerR+tower_height)*tan(0.5*phi_unit),innerR+tower_height));
    vertices2.push_back( G4TwoVector((innerR+tower_height)*tan(0.5*phi_unit),innerR+tower_height));
    
    G4GenericTrap* phiEL = new G4GenericTrap( "phiEL", tower_height/2., vertices2);
    G4LogicalVolume* phiELLog = new G4LogicalVolume(phiEL,Air,"phiELLog");
    
    //Mother volumes and Z ROTATION
    for(int j=0;j<NbOfZRot;j++){ //j<NbOfZRot
        
    //place physical spacing of BARREL
    G4RotationMatrix* rmB = new G4RotationMatrix();
    rmB->rotateZ(M_PI/2.);
    rmB->rotateZ(-j*phi_unit);
    rmB->rotateX(M_PI/2.);
    
    //ER
    G4RotationMatrix* rmER = new G4RotationMatrix();
    rmER->rotateZ(-M_PI/2.);
    rmER->rotateZ(-j*phi_unit);
    rmER->rotateX(M_PI);
    rmER->rotateY(M_PI);
    	
    //EL
    G4RotationMatrix* rmEL = new G4RotationMatrix();
    rmEL->rotateZ(-M_PI/2.);
    rmEL->rotateZ(-j*phi_unit);
    rmEL->rotateX(M_PI);
    	        
//         if (j == NbOfZRot/4 || j == NbOfZRot/4*3) 
//     if (j < NbOfZRot/4  || j > NbOfZRot/4*3) 
//     if (j < NbOfZRot/3  || j > NbOfZRot/3*1.5) 
//         if (j < NbOfZRot/3*0.5  || j > NbOfZRot/3*1.5) 
    if (placeHCAL)
    {
        new G4PVPlacement(rmB,G4ThreeVector((innerR+0.5*length)*cos(j*phi_unit),(innerR+0.5*length)*sin(j*phi_unit),0),phiBLog,"phiDivPhys",worldLV,false,j,false);
        new G4PVPlacement(rmER,G4ThreeVector(0,0,(innerR)*tan(thetaB)+length/2.),phiERLog,"phiERPhys",worldLV,false,j,false);
        new G4PVPlacement(rmEL,G4ThreeVector(0,0,-(innerR)*tan(thetaB)-length/2.),phiELLog,"phiELPhys",worldLV,false,j,false);  
    }

    if (displaySlice)
    {
        if (j == NbOfZRot/4 || j == NbOfZRot/4*3) 
        {
            new G4PVPlacement(rmB,G4ThreeVector((innerR+0.5*length)*cos(j*phi_unit),(innerR+0.5*length)*sin(j*phi_unit),0),phiBLog,"phiDivPhys",worldLV,false,j,false);
            new G4PVPlacement(rmER,G4ThreeVector(0,0,(innerR)*tan(thetaB)+length/2.),phiERLog,"phiERPhys",worldLV,false,j,false);
            new G4PVPlacement(rmEL,G4ThreeVector(0,0,-(innerR)*tan(thetaB)-length/2.),phiELLog,"phiELPhys",worldLV,false,j,false);  
        }
    }
    
        
            
    }
    
    G4double detOpacity = 1;
    G4bool visHCAL = true;
    G4bool wireFrame = false;
    G4bool solidView = true;
    
    // Vis Attributes
    G4VisAttributes* phiVisAttr = new G4VisAttributes(G4Colour(0.8,0.8,0.8, detOpacity));
    phiVisAttr->SetVisibility(false);
    phiBLog->SetVisAttributes(phiVisAttr);
    phiERLog->SetVisAttributes(phiVisAttr);
    phiELLog->SetVisAttributes(phiVisAttr);
    
    
    G4VisAttributes* towerVisAttr = new G4VisAttributes(G4Colour(0,0,1, detOpacity));
//     G4VisAttributes* towerVisAttr = new G4VisAttributes(G4Colour(0.9, 0.4, 0.1, detOpacity));
    towerVisAttr->SetVisibility(visHCAL);
    towerVisAttr->SetDaughtersInvisible(false);
    towerVisAttr->SetForceWireframe(wireFrame);
//     towerVisAttr->SetForceSolid(solidView);
    
    G4VisAttributes* towerVisAttr2 = new G4VisAttributes(G4Colour(0.9, 0.4, 0.1, detOpacity));
    towerVisAttr2->SetVisibility(visHCAL);
    towerVisAttr2->SetDaughtersInvisible(false);
    towerVisAttr2->SetForceWireframe(wireFrame);
//     towerVisAttr2->SetForceSolid(solidView);

    
    G4VisAttributes* PMTVisAttr = new G4VisAttributes(G4Colour(0.3,0.6,0.0, detOpacity));
    PMTVisAttr->SetVisibility(true);
    PMTVisAttr->SetDaughtersInvisible(false);
    PMTVisAttr->SetForceWireframe(wireFrame);
    
    //BARREL....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    dimB = new dimensionB();
    dimB->SetInnerR(innerR);
    dimB->SetTower_height(tower_height);
    dimB->SetNumZRot(NbOfZRot);
    dimB->SetPMTT(PMTT);
    
    // barrel R
    G4cout << "Barrel R..." << G4endl;
    dimB->Rbool(1);
    
    for(int i=0;i<NbOfBarrel;i++){    //i<NbOfBarrel
		
		//cout<<i<<" ";
    	thetaofcenter=fulltheta+deltatheta_barrel[i]/2.;
    	dimB->SetDeltaTheta(deltatheta_barrel[i]);
    	dimB->SetThetaOfCenter(thetaofcenter);
    	dimB->CalBasic();
    	dimB->Getpt(pt);

    	sprintf(name,"tower%d",i+1);
    	tower = new G4Trap("TowerBR",pt);
    	towerLogicalBR[i] = new G4LogicalVolume(tower,cu,name);
    	towerLogicalBR[i]->SetVisAttributes(towerVisAttr2);
    	dimB->Getpt_PMTG(pt);

    	G4RotationMatrix* rm = new G4RotationMatrix();
    	rm->rotateX(-thetaofcenter);
    	G4ThreeVector c = dimB->GetOrigin(0);
    	G4ThreeVector c_new(c.getY(),-c.getZ(),c.getX()-(innerR+0.5*length));
        //placing towers in barrel R
     	new G4PVPlacement(rm,c_new,towerLogicalBR[i],name,phiBLog,false,i+1,checkOverlaps);
        sprintf(name,"PMT%d",volnum);
        
    	dimB->Getpt(pt);
    	sprintf(name,"fiber%d",volnum);
        //VERY IMPORTANT TO PLACE FIBERS
        if (placeFIBERS) fiberBR(i,deltatheta_barrel[i]);

    	fulltheta = fulltheta+deltatheta_barrel[i];
    	volnum++;
    }
    
    // barrel L
    G4cout << "Barrel L..." << G4endl;
    dimB->Rbool(0);
    thetaofcenter=0;
    fulltheta=0;
    
    for(int i=0;i<NbOfBarrel;i++){  //i<NbOfBarrel
    	thetaofcenter=fulltheta+deltatheta_barrel[i]/2.;
    	dimB->SetDeltaTheta(deltatheta_barrel[i]);
    	dimB->SetThetaOfCenter(thetaofcenter);
    	dimB->CalBasic();
    	dimB->Getpt(pt);

    	sprintf(name,"tower%d",-i-1);
    	tower = new G4Trap("TowerBL",pt);
    	towerLogicalBL[i] = new G4LogicalVolume(tower,cu,name);
    	towerLogicalBL[i]->SetVisAttributes(towerVisAttr2);

    	G4RotationMatrix* rm = new G4RotationMatrix();
    	rm->rotateX(thetaofcenter);
    	G4ThreeVector c = dimB->GetOrigin(0);
    	G4ThreeVector c_new(c.getY(),-c.getZ(),c.getX()-(innerR+0.5*length));
    	new G4PVPlacement(rm,c_new,towerLogicalBL[i],name,phiBLog,false,-i-1,checkOverlaps);
       
    	dimB->Getpt(pt);
    	sprintf(name,"fiber%d",volnum);
    	if (placeFIBERS) fiberBL(i,deltatheta_barrel[i]);

    	fulltheta = fulltheta+deltatheta_barrel[i];
    	volnum++;
    }
    
    // ENDCAP....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    dimE = new dimensionE();
    dimE->SetInnerR(innerR);
    dimE->SetTower_height(tower_height);
    dimE->SetNumZRot(NbOfZRot);
    //dimE->SetDeltaTheta(lastdeltatheta);
    dimE->SetPMTT(PMTT);
    
    // endcap R
    G4cout << "Endcap R..." << G4endl;
    dimE->Rbool(1);
    thetaofcenter=0;
    G4double thetaofcenter2=0;
    fulltheta = thetaE;
    
    for(int i=0;i<NbOfEndcap-1;i++){//NbofEndcap-1
    	thetaofcenter=fulltheta-deltatheta_endcap[i]/2.;
        thetaofcenter2=thetaofcenter-deltatheta_endcap[i]/2.-deltatheta_endcap[i+1]/2.;
        
        dimE->SetDeltaTheta(deltatheta_endcap[i]);
    	dimE->SetThetaOfCenter(thetaofcenter);
        dimE->SetDeltaTheta2(deltatheta_endcap[i+1]);
    	dimE->SetThetaOfCenter2(thetaofcenter2);
    	dimE->CalBasic();
    	dimE->Getpt(pt);

    	sprintf(name,"tower%d",NbOfBarrel+i+1);
        
    	tower = new G4Trap("TowerER",pt);
    	towerLogicalER[i] = new G4LogicalVolume(tower,cu,name);
    	towerLogicalER[i]->SetVisAttributes(towerVisAttr);
    	G4RotationMatrix* rm = new G4RotationMatrix();
    	rm->rotateX(thetaofcenter);
    	G4ThreeVector c = dimE->GetOrigin(0);
    	G4ThreeVector c_new(-c.getY(),c.getZ(),c.getX()-(innerR+0.5*length));
        if(i<35)new G4PVPlacement(rm,c_new,towerLogicalER[i],name,phiERLog,false,NbOfBarrel+i+1,checkOverlaps);
        
        dimE->Getpt(pt);
        if (i<35 && placeFIBERS) fiberER(i,deltatheta_endcap[i]);
        fulltheta = fulltheta-deltatheta_endcap[i];
        volnum++;
    }
    
    // endcap L
    G4cout << "Endcap L..." << G4endl;
    dimE->Rbool(0);
    thetaofcenter=0;
    thetaofcenter2=0;
    fulltheta = thetaE;
    
    for(int i=0;i<NbOfEndcap-1;i++){
    	thetaofcenter=fulltheta-deltatheta_endcap[i]/2.;
        thetaofcenter2=thetaofcenter-deltatheta_endcap[i]/2.-deltatheta_endcap[i+1]/2.;
        dimE->SetDeltaTheta(deltatheta_endcap[i]);
    	dimE->SetThetaOfCenter(thetaofcenter);
        dimE->SetDeltaTheta2(deltatheta_endcap[i+1]);
    	dimE->SetThetaOfCenter2(thetaofcenter2);
    	dimE->CalBasic();
    	dimE->Getpt(pt);

    	sprintf(name,"tower%d",volnum);
        tower = new G4Trap("TowerEL",pt);
    	towerLogicalEL[i] = new G4LogicalVolume(tower,cu,name);
    	towerLogicalEL[i]->SetVisAttributes(towerVisAttr);
    	G4RotationMatrix* rm = new G4RotationMatrix();
    	rm->rotateX(thetaofcenter);
    	G4ThreeVector c = dimE->GetOrigin(0);
    	G4ThreeVector c_new(c.getY(),-c.getZ(),c.getX()-(innerR+0.5*length));
        if(i<35)new G4PVPlacement(rm,c_new,towerLogicalEL[i],name,phiELLog,false,-NbOfBarrel-i-1,checkOverlaps);
        
        dimE->Getpt(pt);
        if(i<35 && placeFIBERS) fiberEL(i,deltatheta_endcap[i]);
        fulltheta = fulltheta-deltatheta_endcap[i];
        volnum++;
    }

    
    
    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    //                      Define and place the SCEPCal    
    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    
    
    std::cout << "*****************************************************************************" << std::endl;    
    std::cout << "SCEPCal crystal cross section: " << crystal_size/cm << " x " << crystal_size/cm << " cm2"  << std::endl;
    std::cout << std::endl;
    std::cout << "*****************************************************************************" << std::endl;    
    std::cout << "building SCEPCal barrel with: " << std::endl;
    std::cout << std::endl;
    std::cout << "  --> eta segmentation (nChannels) = " << 2*SCEP_NbOfBarrel <<std::endl;    
    std::cout << "  --> phi segmentation (nChannels) = " << SCEP_NbOfZRot <<std::endl;    
    std::cout << "  --> total barrel channel count   = " << 2*(2*SCEP_NbOfBarrel * SCEP_NbOfZRot) <<std::endl;
    std::cout << std::endl;
    std::cout << "*****************************************************************************" << std::endl;    
    std::cout << std::endl;
    
    
    
    
    G4int SCEP_volnum=0;            
    SCEP_theta_unit=0; 
    SCEP_deltatheta=0; 
    
    G4double SCEP_deltatheta_barrel[180] = {0};
    for(int i=0;i<SCEP_NbOfBarrel;i++)   SCEP_deltatheta_barrel [i] = M_PI/4/(SCEP_NbOfBarrel);
    
    G4double SCEP_deltatheta_endcap[180] = {0};
    for(int i=0;i<SCEP_NbOfEndcap+1;i++) SCEP_deltatheta_endcap [i] = M_PI/4/(SCEP_NbOfEndcap+1);
    
    double SCEP_thetaB = 0;
    for(int i=0;i<SCEP_NbOfBarrel;i++)   SCEP_thetaB += SCEP_deltatheta_barrel[i];
    
    double SCEP_thetaE = 0;    
    for(int i=0;i<SCEP_NbOfEndcap+1;i++) SCEP_thetaE += SCEP_deltatheta_endcap[i];
        
    
    G4VisAttributes* crystalFVisAttr = new G4VisAttributes(G4Colour(0,1,0, detOpacity));
    crystalFVisAttr->SetVisibility(true);
    crystalFVisAttr->SetForceSolid(solidView);
    crystalFVisAttr->SetForceWireframe(wireFrame);

    G4VisAttributes* crystalRVisAttr = new G4VisAttributes(G4Colour(0,0.7,0.3, detOpacity));    
    crystalRVisAttr->SetVisibility(true);
    crystalRVisAttr->SetForceSolid(solidView);
    crystalRVisAttr->SetForceWireframe(wireFrame);
    
    G4VisAttributes* crystalFVisAttrEnd = new G4VisAttributes(G4Colour(0,1,1, detOpacity));    
    crystalFVisAttrEnd->SetVisibility(true);
    crystalFVisAttrEnd->SetForceSolid(solidView);
    crystalFVisAttrEnd->SetForceWireframe(wireFrame);
    
    G4VisAttributes* crystalRVisAttrEnd = new G4VisAttributes(G4Colour(0,0.5,1, detOpacity));    
    crystalRVisAttrEnd->SetVisibility(true);
    crystalRVisAttrEnd->SetForceSolid(solidView);
    crystalRVisAttrEnd->SetForceWireframe(wireFrame);
    
    
    G4VisAttributes* barrelEnvAttr = new G4VisAttributes(cyan);    
    barrelEnvAttr->SetVisibility(false);    
    barrelEnvAttr->SetForceSolid(true);
    barrelEnvAttr->SetForceLineSegmentsPerCircle(50);
    
    G4VisAttributes* endcapEnvAttr = new G4VisAttributes(yellow);    
    endcapEnvAttr->SetVisibility(false);
    endcapEnvAttr->SetForceSolid(true);
    endcapEnvAttr->SetForceLineSegmentsPerCircle(50);

    
    G4VisAttributes* motherEnvelopesVisAttr = new G4VisAttributes(G4Colour(1,1,1));
    motherEnvelopesVisAttr->SetVisibility(false);
//     motherEnvelopesVisAttr->SetForceSolid(true);
    motherEnvelopesVisAttr->SetForceLineSegmentsPerCircle(50);
    
    
    G4VisAttributes* phiSegVisAttr = new G4VisAttributes(G4Colour(0,1,1));
    phiSegVisAttr->SetVisibility(true);
    phiSegVisAttr->SetForceSolid(true);
    phiSegVisAttr->SetForceLineSegmentsPerCircle(50);
    
    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
        
    G4cout << "SCEPCal Barrel ..." << G4endl;
    //SCEPCal BARREL mother volumes for phi slices
    G4double barrel_dPhi = 2.* M_PI * rad;    
    G4Tubs* SCEP_b_temp = new G4Tubs("SCEP_b_temp", SCEP_innerR, SCEP_innerR+SCEP_xtal_L, (SCEP_innerR+SCEP_xtal_L)*tan(SCEP_thetaB), 0., barrel_dPhi);
    
    G4double cone_height = (SCEP_innerR+SCEP_xtal_L)/SCEP_xtal_L*(tan(SCEP_thetaB)*SCEP_xtal_L/2) ;
    
    G4Cons *SCEP_e_temp = new G4Cons("SCEP_e_temp", 0 , 0, 0, SCEP_innerR+SCEP_xtal_L, cone_height, 0., barrel_dPhi);
    G4RotationMatrix* rotSub = new G4RotationMatrix();
    rotSub->rotateX(M_PI);
    
        
//     G4LogicalVolume* SCEP_BLog = new G4LogicalVolume(SCEP_b_temp,Air,"SCEP_BLog");
//     G4LogicalVolume* SCEP_ELog = new G4LogicalVolume(SCEP_e_temp,Air,"SCEP_ELog");
//     new G4PVPlacement(0,G4ThreeVector(0,0,0),SCEP_BLog,"SCEP_BLog_Envelope",worldLV,false,0,checkOverlaps);        
//     new G4PVPlacement(rotSub,G4ThreeVector(0, 0, -(SCEP_innerR+SCEP_xtal_L)*tan(SCEP_thetaB)+cone_height),SCEP_ELog,"SCEP_ELog_EnvelopeL",worldLV,false,0,checkOverlaps);
//     new G4PVPlacement(0,G4ThreeVector(0, 0, (SCEP_innerR+SCEP_xtal_L)*tan(SCEP_thetaB)-cone_height),SCEP_ELog,"SCEP_ELog_EnvelopeR",worldLV,false,0,checkOverlaps);
//     SCEP_ELog->SetVisAttributes(endcapEnvAttr);
//     
    G4VSolid* SCEP_Barrel_temp = new G4SubtractionSolid("SCEP_Barrel_temp", SCEP_b_temp, SCEP_e_temp, rotSub, G4ThreeVector(0, 0, -(SCEP_innerR+SCEP_xtal_L)*tan(SCEP_thetaB)+cone_height));
    G4VSolid* SCEP_Barrel      = new G4SubtractionSolid("SCEP_Barrel",      SCEP_Barrel_temp, SCEP_e_temp, 0, G4ThreeVector(0, 0, (SCEP_innerR+SCEP_xtal_L)*tan(SCEP_thetaB)-cone_height));
        
    
    G4LogicalVolume* SCEP_BLog = new G4LogicalVolume(SCEP_Barrel, Air, "SCEP_BLog");
    new G4PVPlacement(0,G4ThreeVector(0,0,0),SCEP_BLog,"SCEP_BLog_Envelope",fMagneticLogical,false,0,checkOverlaps);                            
    SCEP_BLog->SetVisAttributes(barrelEnvAttr);
    

    G4Trd* SCEP_phiBarrel = new G4Trd("SCEP_phiBarrel",(SCEP_innerR)*tan(0.5*SCEP_phi_unit),(SCEP_innerR+SCEP_xtal_L)*tan(0.5*SCEP_phi_unit),(SCEP_innerR)*tan(SCEP_thetaB),(SCEP_innerR+SCEP_xtal_L)*tan(SCEP_thetaB),0.5*SCEP_xtal_L);
    G4LogicalVolume* SCEP_phiBLog = new G4LogicalVolume(SCEP_phiBarrel,Air,"SCEP_phiBLog");
    SCEP_phiBLog->SetVisAttributes(motherEnvelopesVisAttr);
        
    
   int offset_phi = 0;
   
   if (displaySlice) offset_phi = SCEP_NbOfZRot/4*1;
//     int offset = 100/4*3;

    G4VPVParameterisation* barrelPhiParam = new BarrelPhiParameterisation( SCEP_innerR, SCEP_phi_unit, SCEP_xtal_L, offset_phi);        
    if (placeSCEPCAL) G4VPhysicalVolume* SCEP_phiDivPhys  = new G4PVParameterised( "SCEP_phiDivPhys", SCEP_phiBLog, SCEP_BLog, kUndefined, SCEP_NbOfZRot, barrelPhiParam);

    
    //to display a single slice uncomment the line below
    if (displaySlice) G4VPhysicalVolume* SCEP_phiDivPhys  = new G4PVParameterised( "SCEP_phiDivPhys", SCEP_phiBLog, SCEP_BLog, kUndefined, 1, barrelPhiParam);
    if (displayECALHalfSlice) G4VPhysicalVolume* SCEP_phiDivPhys  = new G4PVParameterised( "SCEP_phiDivPhys", SCEP_phiBLog, SCEP_BLog, kUndefined, 1, barrelPhiParam);
//     if (displaySlice) G4VPhysicalVolume* SCEP_phiDivPhys  = new G4PVParameterised( "SCEP_phiDivPhys", SCEP_phiBLog, SCEP_BLog, kUndefined, SCEP_NbOfZRot/2, barrelPhiParam);
    
    
    
    
//     G4LogicalVolume* SCEP_BLog2 = new G4LogicalVolume(SCEP_Barrel, Air, "SCEP_BLog2");
//     new G4PVPlacement(0,G4ThreeVector(0,0,0),SCEP_BLog2,"SCEP_BLog_Envelope2",fMagneticLogical,false,0,checkOverlaps);                            
//     SCEP_BLog2->SetVisAttributes(barrelEnvAttr);
//     
//     G4LogicalVolume* SCEP_phiBLog2 = new G4LogicalVolume(SCEP_phiBarrel,Air,"SCEP_phiBLog2");
//     offset_phi = SCEP_NbOfBarrel/4*3;
//     G4VPVParameterisation* barrelPhiParam2 = new BarrelPhiParameterisation( SCEP_innerR, SCEP_phi_unit, SCEP_xtal_L, offset_phi);        
//     G4VPhysicalVolume* SCEP_phiDivPhys2  = new G4PVParameterised( "SCEP_phiDivPhys2", SCEP_phiBLog2, SCEP_BLog2, kUndefined, 1, barrelPhiParam2);
//     SCEP_phiBLog2->SetVisAttributes(motherEnvelopesVisAttr);
    
    dimB = new dimensionB();    
    dimB->SetNumZRot(SCEP_NbOfZRot);
        
    //fill crystals along eta
    for (int iSide = -1; iSide < 2; iSide+=2)
    {
        if (iSide == -1 && (displayECALHalfSlice) ) continue;
        char side_name[20];
        if (iSide == -1 ) 
        {
            sprintf(side_name, "L");
            dimB->Rbool(0);
        }
        else
        {
            sprintf(side_name, "R");
            dimB->Rbool(1);
        }
        
        int id_0 = SCEP_NbOfBarrel*(iSide+1)/2;
        SCEP_thetaofcenter=0;
        SCEP_fulltheta = 0;
        
//         std::cout << "doing side: " << iSide << " :: id_0 = " << id_0 << std::endl;
        
        for(int i=0;i<SCEP_NbOfBarrel;i++)
        {   
            
//          std::cout << " SCEP_deltatheta_barrel[ " << i << "] = " << SCEP_deltatheta_barrel[i] << " :: SCEP_thetaofcenter = " << SCEP_thetaofcenter << std::endl;
//             std::cout << "placing barrel phi slice: " << id_0+i << std::endl;
            SCEP_thetaofcenter=SCEP_fulltheta+SCEP_deltatheta_barrel[i]/2.;
            dimB->SetDeltaTheta(SCEP_deltatheta_barrel[i]);
            dimB->SetThetaOfCenter(SCEP_thetaofcenter);
            G4RotationMatrix* rm = new G4RotationMatrix();
            rm->rotateX(-iSide*SCEP_thetaofcenter);
        
            //front segment
            dimB->SetInnerR(SCEP_innerR);
            dimB->SetTower_height(SCEP_xtal_L*FR_X0ratio);
            dimB->CalBasic();
            dimB->Getpt(SCEP_pt);        
        
            G4ThreeVector c = dimB->GetOrigin(0);            	        
//             sprintf(name,"crystalFront%d", iSide*( i+1));
            
//             crystal_F = new G4Trap(name,SCEP_pt);                
            crystal_BF = new G4Trap("crystalFront_B_S",SCEP_pt);                
            sprintf(name,"crystalECALBarrelFront_B%s_L_%d", side_name, i);
            crystalLogicalF_B = new G4LogicalVolume(crystal_BF, SCEPCalMaterial, name);        
            crystalLogicalF_B->SetVisAttributes(crystalFVisAttr);        
    	
//             std::cout << " c.getZ() = " <<  c.getZ() << std::endl;
            G4ThreeVector c_front(c.getY(),-c.getZ(),c.getX()-(SCEP_innerR+0.5*SCEP_xtal_L));
            sprintf(name,"crystalECALBarrelFront_B%s_P_%d", side_name, i);             	
            new G4PVPlacement(rm,c_front,crystalLogicalF_B,name,SCEP_phiBLog,false,iSide*(i+1),checkOverlaps);
//             new G4PVPlacement(rm,c_front,crystalLogicalF_B,name,SCEP_phiBLog2,false,iSide*(i+1),false);
        
            //rear segment            
            dimB->SetInnerR(SCEP_innerR+SCEP_xtal_L*FR_X0ratio*cos(SCEP_thetaofcenter));
            dimB->SetTower_height(SCEP_xtal_L*(1-FR_X0ratio));
            dimB->CalBasic();
            dimB->Getpt(SCEP_pt);
            
            c = dimB->GetOrigin(0);            	
            crystal_BR = new G4Trap("crystalRear_B_S",SCEP_pt);                
            sprintf(name,"crystalECALBarrelRear_B%s_L_%d",side_name,i);
//             std::cout << "name is = " << name << std::endl;
            crystalLogicalR_B = new G4LogicalVolume(crystal_BR, SCEPCalMaterial, name);        
            crystalLogicalR_B->SetVisAttributes(crystalRVisAttr);
    	
            G4ThreeVector c_rear(c.getY(),-c.getZ(),c.getX()-(SCEP_innerR+0.5*SCEP_xtal_L));
            sprintf(name,"crystalECALBarrelRear_B%s_P_%d", side_name,i);             	
            new G4PVPlacement(rm,c_rear,crystalLogicalR_B,name,SCEP_phiBLog,false,iSide*(i+1),checkOverlaps);
//             new G4PVPlacement(rm,c_rear,crystalLogicalR_B,name,SCEP_phiBLog2,false,iSide*(i+1),false);
            
            SCEP_fulltheta = SCEP_fulltheta+SCEP_deltatheta_barrel[i];
            SCEP_volnum++;
            
        }
    }
    


    
    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//     G4cout << "SCEPCal Endcap..." << G4endl;
    
    std::cout << "*****************************************************************************" << std::endl;    
    std::cout << "building SCEPCal endcaps with: " << std::endl;
    std::cout << std::endl;
    std::cout << "  --> eta segmentation (Nb of Rings) = " << SCEP_NbOfEndcap <<std::endl;
                
    int nTotEndcapChannels = 0;
    int crystal_0IDR [SCEP_NbOfEndcap];
    
    G4RotationMatrix* rmEL = new G4RotationMatrix();
    rmEL->rotateX(M_PI/1.);
    
    G4RotationMatrix * rotSlice = new G4RotationMatrix();
    rotSlice->rotateZ(-M_PI/2); 
            
    for (int iRing = 0; iRing < SCEP_NbOfEndcap; iRing++)
    {
     
        G4double this_theta  = M_PI/4+SCEP_deltatheta_endcap[0]/2+SCEP_deltatheta_endcap[iRing]*iRing;
        
        G4double rmin1_F = SCEP_innerR/tan(this_theta+SCEP_deltatheta_endcap[iRing]/2);
        G4double rmax1_F = SCEP_innerR/tan(this_theta-SCEP_deltatheta_endcap[iRing]/2);
        G4double rmin2_F = rmin1_F + SCEP_xtal_L*FR_X0ratio*cos(this_theta+SCEP_deltatheta_endcap[iRing]/2);
        G4double rmax2_F = rmax1_F + SCEP_xtal_L*FR_X0ratio*cos(this_theta-SCEP_deltatheta_endcap[iRing]/2);
        G4double crystal_z_proj = SCEP_xtal_L*FR_X0ratio*sin(this_theta);
        
//         G4double rmin1_R = (SCEP_innerR+SCEP_xtal_L*FR_X0ratio)/tan(this_theta+SCEP_deltatheta_endcap[iRing]/2);
//         G4double rmax1_R = (SCEP_innerR+SCEP_xtal_L*FR_X0ratio)/tan(this_theta-SCEP_deltatheta_endcap[iRing]/2);
//         G4double rmin1_R = rmin1_F + SCEP_xtal_L*FR_X0ratio*cos(this_theta+SCEP_deltatheta_endcap[iRing]/2);
        G4double rmin1_R = rmin2_F;
        G4double rmax1_R = rmax2_F;
        G4double rmin2_R = rmin1_R + SCEP_xtal_L*(1-FR_X0ratio)*cos(this_theta+SCEP_deltatheta_endcap[iRing]/2);
        G4double rmax2_R = rmax1_R + SCEP_xtal_L*(1-FR_X0ratio)*cos(this_theta-SCEP_deltatheta_endcap[iRing]/2);
        G4double crystal_z_proj_R = SCEP_xtal_L*(1-FR_X0ratio)*sin(this_theta);
        
        G4int    n_divisionRing   = rmin1_F*2*M_PI/crystal_size;        
        barrel_dPhi = barrel_dPhi;
        G4double divided_con_dPhi = barrel_dPhi/n_divisionRing;
        
        
     
//         std::cout << "Ring[" << iRing <<  "]: this_theta = "  << this_theta/deg << " :: ring0_rmin1 = " << ring0_rmin1 << " :: ring0_rmax1 = " << ring0_rmax1 << " :: ring0_rmin2 = " << ring0_rmin2 << " :: ring0_rmax2 = " << ring0_rmax2 << " :: crystal_z_proj = " << crystal_z_proj << std::endl;        
              	
//         if (iRing == 0 || iRing == 0 || iRing == SCEP_NbOfEndcap-10) 
        if (iRing < SCEP_NbOfEndcap*0.9-1) 
        {
//             std::cout << "crystal_z_proj[" << iRing << "] = " << crystal_z_proj << std::endl;
            
            std::cout << "  xx --> Nb of channels in ring [" << iRing << "]: " << n_divisionRing << " , starting crystal_0IDR[" << iRing << "] = " << nTotEndcapChannels << std::endl;
            crystal_0IDR[iRing] = nTotEndcapChannels;
            nTotEndcapChannels +=n_divisionRing;
            
            //front crystals
            
            
            
            //right side front
            SCEP_EndcapRing_F = new G4Cons("SCEP_EndcapRing_F", rmin1_F , rmax1_F, rmin2_F, rmax2_F, crystal_z_proj/2., 0., barrel_dPhi);
            SCEP_endcapRingLog_F = new G4LogicalVolume(SCEP_EndcapRing_F,SCEPCalMaterial,"SCEP_ringELog_F");
            SCEP_endcapRingLog_F->SetVisAttributes(motherEnvelopesVisAttr);
                                
            sprintf(name,"SCEP_endcapR_Front_P_ring%d", iRing+1);                   
            if (placeSCEPCAL) new G4PVPlacement(0,G4ThreeVector(0, 0, (SCEP_innerR+SCEP_xtal_L*FR_X0ratio)*tan(SCEP_thetaB)-SCEP_xtal_L*FR_X0ratio+crystal_z_proj/2),SCEP_endcapRingLog_F,name,fMagneticLogical,false,iRing+1,checkOverlaps);                          
            if (displaySlice) new G4PVPlacement(rotSlice,G4ThreeVector(0, 0, (SCEP_innerR+SCEP_xtal_L*FR_X0ratio)*tan(SCEP_thetaB)-SCEP_xtal_L*FR_X0ratio+crystal_z_proj/2),SCEP_endcapRingLog_F,name,fMagneticLogical,false,iRing+1,checkOverlaps);                          
            
            crystal_EF =            new G4Cons("crystal_EF", rmin1_F , rmax1_F, rmin2_F, rmax2_F, crystal_z_proj/2.,        -divided_con_dPhi/2., divided_con_dPhi);
            G4LogicalVolume* crystalF_log = new G4LogicalVolume(crystal_EF, SCEPCalMaterial ,"crystalF_log", 0, 0, 0);        
            sprintf(name,"crystalECALEndcapFront_P_ring%d", iRing+1);
            if (placeSCEPCAL) G4VPhysicalVolume* crystalF_P = new G4PVReplica(name, crystalF_log, SCEP_endcapRingLog_F, kPhi, n_divisionRing, divided_con_dPhi);
            if (displaySlice) G4VPhysicalVolume* crystalF_P = new G4PVReplica(name, crystalF_log, SCEP_endcapRingLog_F, kPhi, 1, divided_con_dPhi);
            crystalF_log->SetVisAttributes(crystalFVisAttr);
//             crystalF_log->SetVisAttributes(crystalFVisAttrEnd);
            
            //left side front
            SCEP_EndcapRing_F_L = new G4Cons("SCEP_EndcapRing_F_L", rmin2_F , rmax2_F, rmin1_F, rmax1_F, crystal_z_proj/2., 0., barrel_dPhi);            
            SCEP_endcapRingLog_F_L = new G4LogicalVolume(SCEP_EndcapRing_F_L,SCEPCalMaterial,"SCEP_ringELog_F_L");
            SCEP_endcapRingLog_F_L->SetVisAttributes(motherEnvelopesVisAttr);
            
            sprintf(name,"SCEP_endcapL_Front_P_ring%d", iRing+1);     
            if (placeSCEPCAL) new G4PVPlacement(0,G4ThreeVector(0, 0, -((SCEP_innerR+SCEP_xtal_L*FR_X0ratio)*tan(SCEP_thetaB)-SCEP_xtal_L*FR_X0ratio+crystal_z_proj/2)),SCEP_endcapRingLog_F_L,name,fMagneticLogical,false,-(iRing+1),checkOverlaps); 
            if (displaySlice) new G4PVPlacement(rotSlice,G4ThreeVector(0, 0, -((SCEP_innerR+SCEP_xtal_L*FR_X0ratio)*tan(SCEP_thetaB)-SCEP_xtal_L*FR_X0ratio+crystal_z_proj/2)),SCEP_endcapRingLog_F_L,name,fMagneticLogical,false,-(iRing+1),checkOverlaps); 
                                                                        
            crystal_EF_L =  new G4Cons("crystal_EF_L", rmin2_F , rmax2_F, rmin1_F, rmax1_F, crystal_z_proj/2.,        -divided_con_dPhi/2., divided_con_dPhi);
            G4LogicalVolume* crystalF_L_log = new G4LogicalVolume(crystal_EF_L, SCEPCalMaterial ,"crystalF_L_log", 0, 0, 0);        
            sprintf(name,"crystalECALEndcapFront_P_ring%d", iRing+1);
            if (placeSCEPCAL) G4VPhysicalVolume* crystalF_L_P = new G4PVReplica(name, crystalF_L_log, SCEP_endcapRingLog_F_L, kPhi, n_divisionRing, divided_con_dPhi);
            if (displaySlice) G4VPhysicalVolume* crystalF_L_P = new G4PVReplica(name, crystalF_L_log, SCEP_endcapRingLog_F_L, kPhi, 1, divided_con_dPhi);
            crystalF_L_log->SetVisAttributes(crystalFVisAttr);
//             crystalF_L_log->SetVisAttributes(crystalFVisAttrEnd);
            
            
            //rear crystals
            
            //right side rear
            SCEP_EndcapRing_R = new G4Cons("SCEP_EndcapRing_R", rmin1_R , rmax1_R, rmin2_R, rmax2_R, crystal_z_proj_R/2., 0., barrel_dPhi);
            SCEP_endcapRingLog_R = new G4LogicalVolume(SCEP_EndcapRing_R,SCEPCalMaterial,"SCEP_ringELog_R");
            SCEP_endcapRingLog_R->SetVisAttributes(motherEnvelopesVisAttr);                    
            
            sprintf(name,"SCEP_endcapR_Rear_P_ring%d", iRing+1);       
            if (placeSCEPCAL) new G4PVPlacement(0,G4ThreeVector(0, 0, (SCEP_innerR+SCEP_xtal_L*FR_X0ratio)*tan(SCEP_thetaB)-SCEP_xtal_L*FR_X0ratio+crystal_z_proj+crystal_z_proj_R/2),SCEP_endcapRingLog_R,name,fMagneticLogical,false,iRing+1,checkOverlaps);                          
            if (displaySlice) new G4PVPlacement(rotSlice,G4ThreeVector(0, 0, (SCEP_innerR+SCEP_xtal_L*FR_X0ratio)*tan(SCEP_thetaB)-SCEP_xtal_L*FR_X0ratio+crystal_z_proj+crystal_z_proj_R/2),SCEP_endcapRingLog_R,name,fMagneticLogical,false,iRing+1,checkOverlaps);                          
            
            crystal_ER =            new G4Cons("crystalR", rmin1_R , rmax1_R, rmin2_R, rmax2_R, crystal_z_proj_R/2.,        -divided_con_dPhi/2., divided_con_dPhi);
            G4LogicalVolume* crystalR_log = new G4LogicalVolume(crystal_ER, SCEPCalMaterial ,"crystalR_log", 0, 0, 0);        
            sprintf(name,"crystalECALEndcapRear_P_ring%d", iRing+1);             	            
            if (placeSCEPCAL) G4VPhysicalVolume* crystalR_P = new G4PVReplica(name, crystalR_log, SCEP_endcapRingLog_R, kPhi, n_divisionRing, divided_con_dPhi);
            if (displaySlice) G4VPhysicalVolume* crystalR_P = new G4PVReplica(name, crystalR_log, SCEP_endcapRingLog_R, kPhi, 1, divided_con_dPhi);
            crystalR_log->SetVisAttributes(crystalRVisAttr);   
//             crystalR_log->SetVisAttributes(crystalRVisAttrEnd);   
            
            //left side rear
            SCEP_EndcapRing_R_L = new G4Cons("SCEP_EndcapRing_R_L", rmin2_R , rmax2_R, rmin1_R, rmax1_R, crystal_z_proj_R/2., 0., barrel_dPhi);
            SCEP_endcapRingLog_R_L = new G4LogicalVolume(SCEP_EndcapRing_R_L,SCEPCalMaterial,"SCEP_ringELog_R_L");
            SCEP_endcapRingLog_R_L->SetVisAttributes(motherEnvelopesVisAttr);                    
            
            sprintf(name,"SCEP_endcapL_Rear_L_P_ring%d", iRing+1);       
            if (placeSCEPCAL) new G4PVPlacement(0,G4ThreeVector(0, 0, -((SCEP_innerR+SCEP_xtal_L*FR_X0ratio)*tan(SCEP_thetaB)-SCEP_xtal_L*FR_X0ratio+crystal_z_proj+crystal_z_proj_R/2)),SCEP_endcapRingLog_R_L,name,fMagneticLogical,false,-(iRing+1),checkOverlaps);                          
            if (displaySlice) new G4PVPlacement(rotSlice,G4ThreeVector(0, 0, -((SCEP_innerR+SCEP_xtal_L*FR_X0ratio)*tan(SCEP_thetaB)-SCEP_xtal_L*FR_X0ratio+crystal_z_proj+crystal_z_proj_R/2)),SCEP_endcapRingLog_R_L,name,fMagneticLogical,false,-(iRing+1),checkOverlaps);                          
                                
            crystal_ER_L =            new G4Cons("crystalR_L", rmin2_R , rmax2_R, rmin1_R, rmax1_R, crystal_z_proj_R/2.,        -divided_con_dPhi/2., divided_con_dPhi);
            G4LogicalVolume* crystalR_L_log = new G4LogicalVolume(crystal_ER_L, SCEPCalMaterial ,"crystalR_L_log", 0, 0, 0);        
            sprintf(name,"crystalECALEndcapRear_P_ring%d", iRing+1);             	            
            if (placeSCEPCAL) G4VPhysicalVolume* crystalR_L_P = new G4PVReplica(name, crystalR_L_log, SCEP_endcapRingLog_R_L, kPhi, n_divisionRing, divided_con_dPhi);
            if (displaySlice) G4VPhysicalVolume* crystalR_L_P = new G4PVReplica(name, crystalR_L_log, SCEP_endcapRingLog_R_L, kPhi, 1, divided_con_dPhi);
            crystalR_L_log->SetVisAttributes(crystalRVisAttr);
//             crystalR_L_log->SetVisAttributes(crystalRVisAttrEnd);
        }
    }
    
    std::cout << "  --> total endcaps (x2) channel count   = " << 2*(2*nTotEndcapChannels) <<std::endl;
    std::cout << std::endl;
    std::cout << "*****************************************************************************" << std::endl;    
    std::cout << std::endl;
    
    
    //print out the array of ring 0ID
    std::cout << "crystal_0IDR[" << int(SCEP_NbOfEndcap*0.9) << "] = {";
    for (int iRing = 0; iRing < SCEP_NbOfEndcap*0.9; iRing++)
    {
        std::cout << crystal_0IDR[iRing] << ", ";
    }
    std::cout << "};";

    
    
    
    
    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    //                Define and place the SCEPCal timing layers 
    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    
    G4double SCEP_Timing_InnerR = 1775.*mm;
    G4double SCEP_Timing_OuterR = 1795.*mm;
    G4double SCEP_Timing_Length = SCEP_Timing_InnerR;
    
        
    //to get a 60x60 mm module
    G4int    nBars      = 20;
    G4double bar_length = 60*mm;
    G4double bar_width  = 3*mm;
    G4double bar_thick  = 3*mm;
    
    
    G4int SCEP_Timing_NbOfPhiRot    = std::floor(2.*M_PI*(SCEP_Timing_InnerR+SCEP_Timing_OuterR)/2./bar_length);
    G4double SCEP_Timing_phi_unit = 2.*M_PI/(G4double)SCEP_Timing_NbOfPhiRot;
    G4int nBarrelTiming_Z         = std::floor(SCEP_Timing_Length/bar_length);
    
    std::cout << std::endl;
    std::cout << "*****************************************************************************" << std::endl;    
    std::cout << "building SCEPCal timing barrel with: " << std::endl;
    std::cout << std::endl;
    std::cout << "  --> Nb of modules along Z = " << nBarrelTiming_Z      << std::endl;
    std::cout << "  --> Nb of rotation in phi = " << SCEP_Timing_NbOfPhiRot << std::endl;
    std::cout << std::endl;
    std::cout << "*****************************************************************************" << std::endl;
    std::cout << std::endl;
    
    //************************************************************************************************************************************
    //barrel envelopes (non pointing, a flat layer)
    
    //full timing barrel mother envelope
    G4Tubs         * Timing_Barrel_S = new G4Tubs("TimingBarrel_S", SCEP_Timing_InnerR, SCEP_Timing_OuterR, SCEP_Timing_Length, 0., 2.* M_PI * rad);
    G4LogicalVolume* Timing_Barrel_L = new G4LogicalVolume(Timing_Barrel_S, Air, "TimingBarrel_L");
    new G4PVPlacement(0,G4ThreeVector(0,0,0),Timing_Barrel_L,"TimingBarrel_P",fMagneticLogical,false,0,checkOverlaps);
//     Timing_Barrel_L->SetVisAttributes(barrelEnvAttr);    
    Timing_Barrel_L->SetVisAttributes(motherEnvelopesVisAttr);

    //barrel timing phi slice to be parameterized and then filled with modules
    G4Tubs* TimingBarrel_PhiSlice_S = new G4Tubs("TimingBarrel_PhiSlice_S",SCEP_Timing_InnerR, SCEP_Timing_OuterR, SCEP_Timing_Length, -SCEP_Timing_phi_unit/2., SCEP_Timing_phi_unit);
    
    G4LogicalVolume* TimingBarrel_PhiSlice_L = new G4LogicalVolume(TimingBarrel_PhiSlice_S,Air,"TimingBarrel_PhiSlice_L");
    TimingBarrel_PhiSlice_L->SetVisAttributes(motherEnvelopesVisAttr);
    
    if (displaySlice) offset_phi = SCEP_Timing_NbOfPhiRot/4*2+1;
    G4VPVParameterisation* barrelTimingPhiParam = new BarrelTimingPhiParameterisation( SCEP_Timing_InnerR, SCEP_Timing_phi_unit, SCEP_Timing_OuterR-SCEP_Timing_InnerR, offset_phi);            
    if (placeTiming)  G4VPhysicalVolume* SCEP_Timing_phiDivPhys  = new G4PVParameterised( "SCEP_TimingBarrel_phiDivPhys", TimingBarrel_PhiSlice_L, Timing_Barrel_L, kUndefined, SCEP_Timing_NbOfPhiRot, barrelTimingPhiParam);
    if (displaySlice) G4VPhysicalVolume* SCEP_Timing_phiDivPhys  = new G4PVParameterised( "SCEP_TimingBarrel_phiDivPhys", TimingBarrel_PhiSlice_L, Timing_Barrel_L, kUndefined, 1, barrelTimingPhiParam);
    
    
    
    //endcap envelopes (non pointing, flat disk, filled as lines of modules)
    G4Tubs         * Timing_Endcap_S = new G4Tubs("TimingEndcap_S", 0, SCEP_Timing_OuterR, (SCEP_Timing_OuterR-SCEP_Timing_InnerR)/2., 0., 2.* M_PI * rad);            
    G4LogicalVolume* Timing_Endcap_L = new G4LogicalVolume(Timing_Endcap_S, Air, "TimingEndcap_L");

    G4RotationMatrix * rotEndcap = new G4RotationMatrix();
    rotEndcap->rotateX(M_PI);   
    new G4PVPlacement(rotEndcap,G4ThreeVector(0,0, (SCEP_Timing_OuterR+SCEP_Timing_InnerR)/2.),Timing_Endcap_L,"TimingEndcap_P_R",fMagneticLogical,false, 1,checkOverlaps);
    new G4PVPlacement(0,G4ThreeVector(0,0,-(SCEP_Timing_OuterR+SCEP_Timing_InnerR)/2.),Timing_Endcap_L,"TimingEndcap_P_L",fMagneticLogical,false,-1,checkOverlaps);           
//     Timing_Endcap_L->SetVisAttributes(barrelEnvAttr);    
    Timing_Endcap_L->SetVisAttributes(motherEnvelopesVisAttr);    
    

            
    //************************************************************************************************************************************
    //mother module envelope of 60x60x6mm^3 containing the two timing layers
    G4Box           * Timing_Module_Env_S = new G4Box          ("Timing_Module_Env_S", nBars*bar_width/2., nBars*bar_width/2., bar_thick);
    G4LogicalVolume * Timing_Module_Env_L = new G4LogicalVolume(Timing_Module_Env_S, Air, "Timing_Module_Env_L");
    G4RotationMatrix * rotModule = new G4RotationMatrix();
    rotModule->rotateY(M_PI/2.);        
    Timing_Module_Env_L->SetVisAttributes(crystalRVisAttrEnd);
//     rotModule->rotateZ(M_PI);            

    
    //placing modules in barrel
    sprintf(name,"TimingBarrelModule");             	            
    for (int iZ = 0; iZ < nBarrelTiming_Z; iZ++)
//         for (int iZ = 0; iZ < 1; iZ++)
//     for (int iZ = 0; iZ < 1; iZ++)        
    {        
//       if (placeTiming)      
          new G4PVPlacement(rotModule, G4ThreeVector( (SCEP_Timing_InnerR+SCEP_Timing_OuterR)/2., 0,  bar_length*(iZ+0.5)), Timing_Module_Env_L, name, TimingBarrel_PhiSlice_L, false,    iZ+1, checkOverlaps);
//       if (placeTiming)     
          new G4PVPlacement(rotModule, G4ThreeVector( (SCEP_Timing_InnerR+SCEP_Timing_OuterR)/2., 0, -bar_length*(iZ+0.5)), Timing_Module_Env_L, name, TimingBarrel_PhiSlice_L, false,  -(iZ+1), checkOverlaps);        
    }
    
    
    //placing modules in endcaps
    G4double SCEP_Timing_Endcap_InnerR = 200.*mm;
    G4int nEndcapModulePerLine = std::floor(SCEP_Timing_OuterR/bar_length*2);
    
    std::cout << "building SCEPCal timing endcaps with: " << std::endl;
    std::cout << std::endl;
    std::cout << "  --> Nb of modules per line = " << nEndcapModulePerLine      << std::endl;
    std::cout << std::endl;
    std::cout << "*****************************************************************************" << std::endl;
    std::cout << std::endl;
    
    G4RotationMatrix * rotModuleEndcap = new G4RotationMatrix();
    rotModuleEndcap->rotateZ(M_PI/2.);        
//     rotModuleEndcap->rotateZ(M_PI);        
    
    sprintf(name,"TimingEndcapModule");             	            
    for (int iX = 0; iX < nEndcapModulePerLine; iX++)
//     for (int iX = 0; iX < nEndcapModulePerLine/2; iX++)
    {
        for (int iY = 0; iY < nEndcapModulePerLine; iY++)
//         for (int iY = 0; iY < nEndcapModulePerLine/6; iY++)
        {
                        
            double posX   = (iX+0.5-nEndcapModulePerLine/2)*bar_length;
            double posY   = (iY+0.5-nEndcapModulePerLine/2)*bar_length;
            double radius = sqrt(posX*posX+posY*posY);
            
            double radius_min = radius-bar_length*sqrt(2)/2.;
            double radius_max = radius+bar_length*sqrt(2)/2.;
            
            if (radius_min>SCEP_Timing_Endcap_InnerR && radius_max< SCEP_Timing_OuterR)
            {
	      if (placeTiming)    new G4PVPlacement(rotModuleEndcap, G4ThreeVector(posX, posY, 0.), Timing_Module_Env_L, name, Timing_Endcap_L, false, iX*nEndcapModulePerLine+iY, checkOverlaps);
              if (displaySlice && iX ==nEndcapModulePerLine/2) new G4PVPlacement(rotModuleEndcap, G4ThreeVector(posX, posY, 0.), Timing_Module_Env_L, name, Timing_Endcap_L, false, iX*nEndcapModulePerLine+iY, checkOverlaps);
              
            }
            
        }
    }
    if (displayTimingModule) new G4PVPlacement(rotModuleEndcap, G4ThreeVector(0, 0, 0.), Timing_Module_Env_L, name, worldLV, false, 0, checkOverlaps);
    
    

    //************************************************************************************************************************************    
    //individual layers
    G4Box           * Timing_Layer_Env_S = new G4Box          ("Timing_Layer_Env_S", nBars*bar_width/2., nBars*bar_width/2., bar_thick/2.);
    G4LogicalVolume * Timing_Layer_Env_L = new G4LogicalVolume(Timing_Layer_Env_S, SCEP_TimingMaterial, "Timing_Module_Env_L");
    
    //front layer
    sprintf(name,"crystalTimingLayer_%d", 1);             	            
    new G4PVPlacement(0, G4ThreeVector(0, 0, bar_thick/2.), Timing_Layer_Env_L, name, Timing_Module_Env_L, false, 1, checkOverlaps);
    
    //rear layer
    sprintf(name,"crystalTimingLayer_%d", 2);             	            
    G4RotationMatrix * rmTiming = new G4RotationMatrix();
    rmTiming->rotateZ(M_PI/2.);        
    new G4PVPlacement(rmTiming, G4ThreeVector(0, 0, -bar_thick/2.), Timing_Layer_Env_L, name, Timing_Module_Env_L, false, 2, checkOverlaps);                          
    Timing_Layer_Env_L->SetVisAttributes(crystalRVisAttrEnd);
    
    
    //************************************************************************************************************************************
    //create replicas of single bars within layers
    G4Box * crystalBar_S = new G4Box("crystalBar_S", bar_width/2., bar_length/2., bar_thick/2.);
    G4LogicalVolume* crystalBar_L = new G4LogicalVolume(crystalBar_S, SCEP_TimingMaterial ,"crystalBar_L", 0, 0, 0);        
        
    G4VPhysicalVolume* crystalBar_P = new G4PVReplica("crystalTimingBar_P", crystalBar_L, Timing_Layer_Env_L, kXAxis, nBars, bar_width);
    crystalBar_L->SetVisAttributes(crystalRVisAttrEnd);                
    
    return worldPV;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4DetectorConstruction::ConstructSDandField()
{
    // Create global magnetic field messenger,
    // Uniform magnetic field is then created automatically if
    // the field value is not zero
	/*G4ThreeVector fieldValue = G4ThreeVector();
	fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
	fMagFieldMessenger->SetVerboseLevel(1);

    // Register the field messenger for deleting
	G4AutoDelete::Register(fMagFieldMessenger);*/
	
	
	// magnetic field ----------------------------------------------------------
  fMagneticField = new B4MagneticField();
  fFieldMgr = new G4FieldManager();
  fFieldMgr->SetDetectorField(fMagneticField);
  fFieldMgr->CreateChordFinder(fMagneticField);
  G4bool forceToAllDaughters = true;
  fMagneticLogical->SetFieldManager(fFieldMgr, forceToAllDaughters);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......




void B4DetectorConstruction::fiberBR(G4int i,G4double deltatheta_){

    vector<G4double> temp_x;
    vector<G4double> temp_y;// vector for grid calculation
    
    temp_x.clear();
    temp_y.clear();
    temp_x.push_back(0.*mm);
    temp_y.push_back(0.*mm);
	
	int T_index=i+1;
    
    int fiber_N = 1500;
    
    for(int j = 0 ; j<fiber_N;j++){
    	temp_x.push_back(0.*mm);
    	temp_y.push_back(1.5*(j+1)*mm);
    	temp_x.push_back(0.*mm);
    	temp_y.push_back(-1.5*(j+1)*mm);
    }
    for(int i=0; i<fiber_N;i++){
    	temp_x.push_back(1.5*(i+1)*mm);
    	temp_y.push_back(0.*mm);
    	temp_x.push_back(-1.5*(i+1)*mm);
    	temp_y.push_back(0.*mm);
    }
    for(int j = 0;j<fiber_N;j++){
    	for(int i = 0;i<fiber_N;i++){
    		temp_x.push_back(1.5*(i+1)*mm);
    		temp_y.push_back(1.5*(j+1)*mm);
    		temp_x.push_back(-1.5*(i+1)*mm);
    		temp_y.push_back(-1.5*(j+1)*mm);
    		temp_x.push_back(-1.5*(i+1)*mm);
    		temp_y.push_back(1.5*(j+1)*mm);
    		temp_x.push_back(1.5*(i+1)*mm);
    		temp_y.push_back(-1.5*(j+1)*mm);
    	}
    }
    
    G4double center_x;
    G4double center_y;
    
    int hi1;
    int hi2;
    int hi1reminder; // 0 -> even 1-> odd
    int hi2reminder;
    
    int numx;
    int numy;
    int reminderx;
    int remindery;
    
    G4ThreeVector v1 = dimB->GetV1();
    G4ThreeVector v2 = dimB->GetV2();
    G4ThreeVector v3 = dimB->GetV3();
    G4ThreeVector v4 = dimB->GetV4();
    
    G4double innerSide_half = dimB->GetInnerR_new()*tan(deltatheta_/2.);
    G4double outerSide_half = (dimB->GetInnerR_new()+tower_height)*tan(deltatheta_/2.);
    
    G4double theta_bc=0;
    
    
    //1. 4 types of grid coordinate(x,y) =  (o,o),(e,o),(o,e),(e,e) { e= 0, o=1 }
    int type_x_BR;
    int type_y_BR;
    
    numx = (int)(((v4.getX()*tan(phi_unit/2.)*2)-1.*mm)/1.5*mm);
    numy = (int)((outerSide_half*2-1.*mm)/(1.5*mm));
    reminderx = numx%2;
    remindery = numy%2;
    if(reminderx == 1) type_x_BR=0;
    if(reminderx == 0) type_x_BR=1;
    if(remindery == 1) type_y_BR=0;
    if(remindery == 0) type_y_BR=1;
    
    ////2. aplying the boundary conditions, reject the fibre which near the boundary of tower surface, get coord of cetres of fibre
    //prepare the grid vector : v2 ~ 114 ( #76 ) , outerSide ~ 132 (#88)
    
    ////////////////////////////////////////////
    vector<G4double> center_x_BR;
    vector<G4double> center_y_BR;
    vector<G4int> bool_cfiber_BR;// 0 -> c fiber 1 s fiber
    
    // select the type of grid coord
    theta_bc=atan(2*outerSide_half/((v2.getX()-v4.getX())*tan(phi_unit/2.)));
    
    if(type_x_BR==1&&type_y_BR==1){//(o,o)
    	for(int num = 0;num<temp_x.size();num++){
    		center_x = temp_x.at(num);
    		center_y = temp_y.at(num);
    		hi1 = fabs(temp_x.at(num))/(1.5*mm);
    		hi2 = fabs(temp_y.at(num))/(1.5*mm);
            hi1reminder = hi1%2; // 0 -> even 1-> odd
            hi2reminder = hi2%2;
            
            
            if(((center_y-cos(theta_bc)*clad_C_rMax)>(2*2*outerSide_half/(2*v2.getX()*tan(phi_unit/2.)-2*v4.getX()*tan(phi_unit/2.))*(center_x+sin(theta_bc)*clad_C_rMax-v2.getX()*tan(phi_unit/2.))+outerSide_half))
            	&&
            	((center_y+clad_C_rMax)<(outerSide_half))
            	&&
            	((center_y-clad_C_rMax)>(-outerSide_half))
            	&&
            	((center_y-cos(theta_bc)*clad_C_rMax)>(2*2*outerSide_half/(2*v4.getX()*tan(phi_unit/2.)-2*v2.getX()*tan(phi_unit/2.))*(center_x-sin(theta_bc)*clad_C_rMax+v2.getX()*tan(phi_unit/2.))+outerSide_half
            		)))
            {
            	center_x_BR.push_back(center_x);
            	center_y_BR.push_back(center_y);
            	if(((hi1reminder==0)&&(hi2reminder==0))||((hi1reminder==1)&&(hi2reminder==1)))
            		bool_cfiber_BR.push_back(0);
            	else bool_cfiber_BR.push_back(1);
            }
            
        }
    }
    
    if(type_x_BR==0&&type_y_BR==0){//(e,e)

    	for(int num = 0;num<temp_x.size();num++){
    		center_x = temp_x.at(num)+0.75*mm;
    		center_y = temp_y.at(num)-0.75*mm;

    		hi1 = fabs(temp_x.at(num))/(1.5*mm);
    		hi2 = fabs(temp_y.at(num))/(1.5*mm);
            hi1reminder = hi1%2; // 0 -> even 1-> odd
            hi2reminder = hi2%2;
            
            if(((center_y-cos(theta_bc)*clad_C_rMax)>(2*2*outerSide_half/(2*v2.getX()*tan(phi_unit/2.)-2*v4.getX()*tan(phi_unit/2.))*(center_x+sin(theta_bc)*clad_C_rMax-v2.getX()*tan(phi_unit/2.))+outerSide_half))
            	&&
            	((center_y+clad_C_rMax)<(outerSide_half))
            	&&
            	((center_y-clad_C_rMax)>(-outerSide_half))
            	&&
            	((center_y-cos(theta_bc)*clad_C_rMax)>(2*2*outerSide_half/(2*v4.getX()*tan(phi_unit/2.)-2*v2.getX()*tan(phi_unit/2.))*(center_x-sin(theta_bc)*clad_C_rMax+v2.getX()*tan(phi_unit/2.))+outerSide_half
            		)))
            {
            	center_x_BR.push_back(center_x);
            	center_y_BR.push_back(center_y);
            	if(((hi1reminder==0)&&(hi2reminder==0))||((hi1reminder==1)&&(hi2reminder==1)))
            		bool_cfiber_BR.push_back(0);
            	else bool_cfiber_BR.push_back(1);
            }
        }
    }
    
    
    if(type_x_BR==0&&type_y_BR==1){//(e,o)

    	for(int num = 0;num<temp_x.size();num++){
    		center_x = temp_x.at(num)+0.75*mm;
    		center_y = temp_y.at(num);

    		hi1 = fabs(temp_x.at(num))/(1.5*mm);
    		hi2 = fabs(temp_y.at(num))/(1.5*mm);
            hi1reminder = hi1%2; // 0 -> even 1-> odd
            hi2reminder = hi2%2;
            
            if(((center_y-cos(theta_bc)*clad_C_rMax)>(2*2*outerSide_half/(2*v2.getX()*tan(phi_unit/2.)-2*v4.getX()*tan(phi_unit/2.))*(center_x+sin(theta_bc)*clad_C_rMax-v2.getX()*tan(phi_unit/2.))+outerSide_half))
            	&&
            	((center_y+clad_C_rMax)<(outerSide_half))
            	&&
            	((center_y-clad_C_rMax)>(-outerSide_half))
            	&&
            	((center_y-cos(theta_bc)*clad_C_rMax)>(2*2*outerSide_half/(2*v4.getX()*tan(phi_unit/2.)-2*v2.getX()*tan(phi_unit/2.))*(center_x-sin(theta_bc)*clad_C_rMax+v2.getX()*tan(phi_unit/2.))+outerSide_half
            		)))
            {
            	center_x_BR.push_back(center_x);
            	center_y_BR.push_back(center_y);
            	if(((hi1reminder==0)&&(hi2reminder==0))||((hi1reminder==1)&&(hi2reminder==1)))
            		bool_cfiber_BR.push_back(0);
            	else bool_cfiber_BR.push_back(1);
            }
            
        }
    }
    
    if(type_x_BR==1&&type_y_BR==0){//(o,e)

    	for(int num = 0;num<temp_x.size();num++){
    		center_x = temp_x.at(num);
    		center_y = temp_y.at(num)-0.75*mm;

    		hi1 = fabs(temp_x.at(num))/(1.5*mm);
    		hi2 = fabs(temp_y.at(num))/(1.5*mm);
            hi1reminder = hi1%2; // 0 -> even 1-> odd
            hi2reminder = hi2%2;
            
            if(((center_y-cos(theta_bc)*clad_C_rMax)>(2*2*outerSide_half/(2*v2.getX()*tan(phi_unit/2.)-2*v4.getX()*tan(phi_unit/2.))*(center_x+sin(theta_bc)*clad_C_rMax-v2.getX()*tan(phi_unit/2.))+outerSide_half))
            	&&
            	((center_y+clad_C_rMax)<(outerSide_half))
            	&&
            	((center_y-clad_C_rMax)>(-outerSide_half))
            	&&
            	((center_y-cos(theta_bc)*clad_C_rMax)>(2*2*outerSide_half/(2*v4.getX()*tan(phi_unit/2.)-2*v2.getX()*tan(phi_unit/2.))*(center_x-sin(theta_bc)*clad_C_rMax+v2.getX()*tan(phi_unit/2.))+outerSide_half
            		)))
            {
            	center_x_BR.push_back(center_x);
            	center_y_BR.push_back(center_y);
            	if(((hi1reminder==0)&&(hi2reminder==0))||((hi1reminder==1)&&(hi2reminder==1)))
            		bool_cfiber_BR.push_back(0);
            	else bool_cfiber_BR.push_back(1);
            }
            
        }
    }
	
	//int counter=0;
    for(int j = 0; j<center_x_BR.size();j++){
    	double z = tower->GetZHalfLength();
    	bool outside = false;
    	for(int ip = 0;ip<4;ip++){
    		TrapSidePlane plane = tower->GetSidePlane(ip);
    		double zpoint = (-plane.a*center_x_BR.at(j)-plane.b*center_y_BR.at(j)-plane.d)/plane.c;
    		outside = (tower->Inside(G4ThreeVector(center_x_BR.at(j),center_y_BR.at(j),zpoint))==kOutside);
    		//G4cout << ip << ": " << plane.a << " " << plane.b << " " << plane.c << " " << plane.d << " " << outside << G4endl;
    		if(!outside){
    			G4ThreeVector normal = tower->SurfaceNormal(G4ThreeVector(center_x_BR.at(j),center_y_BR.at(j),zpoint));
    			double angle = normal.angle(G4ThreeVector(0,0,-1));
    			double shift = fabs(clad_C_rMax/tan(0.5*M_PI-angle));
    			int length = z - zpoint - shift;
    			//G4cout<<" "<< length <<" ";
				sprintf(name,"%d",length);
    			if(length>=1&&length<=2*z){
					int f=(j+1)+fiber_N*fiber_N*(T_index-1);
    				new G4PVPlacement(0,G4ThreeVector(center_x_BR.at(j),center_y_BR.at(j),z-0.5*length),
    					(bool_cfiber_BR.at(j)==0)?fiberCLog[length]:fiberSLog[length],name,towerLogicalBR[i],false,f,false);
					//counter++;

    			}
    			break;
    		}
    	}
    	if(outside) {
			sprintf(name,"%d",2000);
			//cout<<" 2000 ";
			int f=(j+1)+fiber_N*fiber_N*(T_index-1);
			new G4PVPlacement(0,G4ThreeVector(center_x_BR.at(j),center_y_BR.at(j),0),
    		(bool_cfiber_BR.at(j)==0)?fiberCLog[2*z]:fiberSLog[2*z],name,towerLogicalBR[i],false,f,false);
			//counter++;
		
        }   
    }//G4cout<<counter<<std::endl;
}

void B4DetectorConstruction::fiberBL(G4int i, G4double deltatheta_){
	vector<G4double> temp_x;
    vector<G4double> temp_y;// vector for grid calculation
    
    temp_x.clear();
    temp_y.clear();
    temp_x.push_back(0.*mm);
    temp_y.push_back(0.*mm);
    
	int T_index=i+1;
	
    int fiber_N=1500;
    
    for(int j = 0 ; j<fiber_N;j++){
    	temp_x.push_back(0.*mm);
    	temp_y.push_back(1.5*(j+1)*mm);
    	temp_x.push_back(0.*mm);
    	temp_y.push_back(-1.5*(j+1)*mm);
    }
    for(int i=0; i<fiber_N;i++){
    	temp_x.push_back(1.5*(i+1)*mm);
    	temp_y.push_back(0.*mm);
    	temp_x.push_back(-1.5*(i+1)*mm);
    	temp_y.push_back(0.*mm);
    }
    for(int j = 0;j<fiber_N;j++){
    	for(int i = 0;i<fiber_N;i++){
    		temp_x.push_back(1.5*(i+1)*mm);
    		temp_y.push_back(1.5*(j+1)*mm);
    		temp_x.push_back(-1.5*(i+1)*mm);
    		temp_y.push_back(-1.5*(j+1)*mm);
    		temp_x.push_back(-1.5*(i+1)*mm);
    		temp_y.push_back(1.5*(j+1)*mm);
    		temp_x.push_back(1.5*(i+1)*mm);
    		temp_y.push_back(-1.5*(j+1)*mm);
    	}
    }
    
    G4double center_x;
    G4double center_y;
    int hi1;
    int hi2;
    int hi1reminder; // 0 -> even 1-> odd
    int hi2reminder;
    
    int numx;
    int numy;
    int reminderx;
    int remindery;
    
    G4ThreeVector v1 = dimB->GetV1();
    G4ThreeVector v2 = dimB->GetV2();
    G4ThreeVector v3 = dimB->GetV3();
    G4ThreeVector v4 = dimB->GetV4();
    
    G4double innerSide_half = dimB->GetInnerR_new()*tan(deltatheta_/2.);
    G4double outerSide_half = (dimB->GetInnerR_new()+tower_height)*tan(deltatheta_/2.);
    
    G4double theta_bc=0;
    /////////////////////
    //1. 4 types of grid coordinate(x,y) =  (o,o),(e,o),(o,e),(e,e) { e= 0, o=1 }
    
    int type_x_BL;
    int type_y_BL;
    
    numx = (int)(((v4.getX()*tan(phi_unit/2.)*2)-1.*mm)/1.5*mm);
    numy = (int)((outerSide_half*2-1.*mm)/(1.5*mm));
    reminderx = numx%2;
    remindery = numy%2;
    
    if(reminderx == 1) type_x_BL=0;
    if(reminderx == 0) type_x_BL=1;
    if(remindery == 1) type_y_BL=0;
    if(remindery == 0) type_y_BL=1;
    
    
    ////2. aplying the boundary conditions, reject the fibre which near the boundary of tower surface, get coord of cetres of fibre
    vector<G4double> center_x_BL;
    vector<G4double> center_y_BL;
    vector<G4int> bool_cfiber_BL;// 0 -> c fiber 1 s fiber
    
    
    // select the type of grid coord
    theta_bc=atan(2*outerSide_half/((v2.getX()-v4.getX())*tan(phi_unit/2.)));
    if(type_x_BL==1&&type_y_BL==1){//(o,o)
    	for(int num = 0;num<temp_x.size();num++){
    		center_x = temp_x.at(num);
    		center_y = temp_y.at(num);
    		hi1 = fabs(temp_x.at(num))/(1.5*mm);
    		hi2 = fabs(temp_y.at(num))/(1.5*mm);
            hi1reminder = hi1%2; // 0 -> even 1-> odd
            hi2reminder = hi2%2;
            
            if(((center_y+cos(theta_bc)*clad_C_rMax)<(2*2*outerSide_half/(2*v4.getX()*tan(phi_unit/2.)-2*v2.getX()*tan(phi_unit/2.))*(center_x+sin(theta_bc)*clad_C_rMax-v4.getX()*tan(phi_unit/2.))+outerSide_half))
            	&&
            	((center_y+clad_C_rMax)<(outerSide_half))
            	&&
            	((center_y-clad_C_rMax)>(-outerSide_half))
            	&&
            	((center_y+cos(theta_bc)*clad_C_rMax)<(2*2*outerSide_half/(2*v2.getX()*tan(phi_unit/2.)-2*v4.getX()*tan(phi_unit/2.))*(center_x-sin(theta_bc)*clad_C_rMax+v4.getX()*tan(phi_unit/2.))+outerSide_half
            		)))
            {
            	center_x_BL.push_back(center_x);
            	center_y_BL.push_back(center_y);
            	if(((hi1reminder==0)&&(hi2reminder==0))||((hi1reminder==1)&&(hi2reminder==1)))
            		bool_cfiber_BL.push_back(0);
            	else bool_cfiber_BL.push_back(1);
            }
            
        }
    }
    
    if(type_x_BL==0&&type_y_BL==0){//(e,e)

    	for(int num = 0;num<temp_x.size();num++){
    		center_x = temp_x.at(num)+0.75*mm;
    		center_y = temp_y.at(num)-0.75*mm;

    		hi1 = fabs(temp_x.at(num))/(1.5*mm);
    		hi2 = fabs(temp_y.at(num))/(1.5*mm);
            hi1reminder = hi1%2; // 0 -> even 1-> odd
            hi2reminder = hi2%2;
            if(((center_y+cos(theta_bc)*clad_C_rMax)<(2*2*outerSide_half/(2*v4.getX()*tan(phi_unit/2.)-2*v2.getX()*tan(phi_unit/2.))*(center_x+sin(theta_bc)*clad_C_rMax-v4.getX()*tan(phi_unit/2.))+outerSide_half))
            	&&
            	((center_y+clad_C_rMax)<(outerSide_half))
            	&&
            	((center_y-clad_C_rMax)>(-outerSide_half))
            	&&
            	((center_y+cos(theta_bc)*clad_C_rMax)<(2*2*outerSide_half/(2*v2.getX()*tan(phi_unit/2.)-2*v4.getX()*tan(phi_unit/2.))*(center_x-sin(theta_bc)*clad_C_rMax+v4.getX()*tan(phi_unit/2.))+outerSide_half
            		)))
            {
            	center_x_BL.push_back(center_x);
            	center_y_BL.push_back(center_y);
            	if(((hi1reminder==0)&&(hi2reminder==0))||((hi1reminder==1)&&(hi2reminder==1)))
            		bool_cfiber_BL.push_back(0);
            	else bool_cfiber_BL.push_back(1);
            }
        }
    }
    
    
    if(type_x_BL==0&&type_y_BL==1){//(e,o)

    	for(int num = 0;num<temp_x.size();num++){
    		center_x = temp_x.at(num)+0.75*mm;
    		center_y = temp_y.at(num);

    		hi1 = fabs(temp_x.at(num))/(1.5*mm);
    		hi2 = fabs(temp_y.at(num))/(1.5*mm);
            hi1reminder = hi1%2; // 0 -> even 1-> odd
            hi2reminder = hi2%2;
            
            if(((center_y+cos(theta_bc)*clad_C_rMax)<(2*2*outerSide_half/(2*v4.getX()*tan(phi_unit/2.)-2*v2.getX()*tan(phi_unit/2.))*(center_x+sin(theta_bc)*clad_C_rMax-v4.getX()*tan(phi_unit/2.))+outerSide_half))
            	&&
            	((center_y+clad_C_rMax)<(outerSide_half))
            	&&
            	((center_y-clad_C_rMax)>(-outerSide_half))
            	&&
            	((center_y+cos(theta_bc)*clad_C_rMax)<(2*2*outerSide_half/(2*v2.getX()*tan(phi_unit/2.)-2*v4.getX()*tan(phi_unit/2.))*(center_x-sin(theta_bc)*clad_C_rMax+v4.getX()*tan(phi_unit/2.))+outerSide_half
            		)))
            {
            	center_x_BL.push_back(center_x);
            	center_y_BL.push_back(center_y);
            	if(((hi1reminder==0)&&(hi2reminder==0))||((hi1reminder==1)&&(hi2reminder==1)))
            		bool_cfiber_BL.push_back(0);
            	else bool_cfiber_BL.push_back(1);
            }
            
        }
    }
    if(type_x_BL==1&&type_y_BL==0){//(o,e)

    	for(int num = 0;num<temp_x.size();num++){
    		center_x = temp_x.at(num);
    		center_y = temp_y.at(num)-0.75*mm;

    		hi1 = fabs(temp_x.at(num))/(1.5*mm);
    		hi2 = fabs(temp_y.at(num))/(1.5*mm);
            hi1reminder = hi1%2; // 0 -> even 1-> odd
            hi2reminder = hi2%2;
            
            if(((center_y+cos(theta_bc)*clad_C_rMax)<(2*2*outerSide_half/(2*v4.getX()*tan(phi_unit/2.)-2*v2.getX()*tan(phi_unit/2.))*(center_x+sin(theta_bc)*clad_C_rMax-v4.getX()*tan(phi_unit/2.))+outerSide_half))
            	&&
            	((center_y+clad_C_rMax)<(outerSide_half))
            	&&
            	((center_y-clad_C_rMax)>(-outerSide_half))
            	&&
            	((center_y+cos(theta_bc)*clad_C_rMax)<(2*2*outerSide_half/(2*v2.getX()*tan(phi_unit/2.)-2*v4.getX()*tan(phi_unit/2.))*(center_x-sin(theta_bc)*clad_C_rMax+v4.getX()*tan(phi_unit/2.))+outerSide_half
            		)))
            {
            	center_x_BL.push_back(center_x);
            	center_y_BL.push_back(center_y);
            	if(((hi1reminder==0)&&(hi2reminder==0))||((hi1reminder==1)&&(hi2reminder==1)))
            		bool_cfiber_BL.push_back(0);
            	else bool_cfiber_BL.push_back(1);
            }
            
        }
    }
    
    for(int j = 0; j<center_x_BL.size();j++){
        // determine z value for center of fibre
    	double z = tower->GetZHalfLength();
    	bool outside = false;
    	for(int ip = 0;ip<4;ip++){
    		TrapSidePlane plane = tower->GetSidePlane(ip);
    		double zpoint = (-plane.a*center_x_BL.at(j)-plane.b*center_y_BL.at(j)-plane.d)/plane.c;
    		outside = (tower->Inside(G4ThreeVector(center_x_BL.at(j),center_y_BL.at(j),zpoint))==kOutside);
    		if(!outside){
    			G4ThreeVector normal = tower->SurfaceNormal(G4ThreeVector(center_x_BL.at(j),center_y_BL.at(j),zpoint));
    			double angle = normal.angle(G4ThreeVector(0,0,-1));
    			double shift = fabs(clad_C_rMax/tan(0.5*M_PI-angle));
    			int length = z - zpoint - shift;
				sprintf(name,"%d",length);
    			if(length>=1&&length<=2*z){
					int f=-((j+1)+fiber_N*fiber_N*(T_index-1));
    				new G4PVPlacement(0,G4ThreeVector(center_x_BL.at(j),center_y_BL.at(j),z-0.5*length),
    					(bool_cfiber_BL.at(j)==0)?fiberCLog[length]:fiberSLog[length],name,towerLogicalBL[i],false,f,false);
    			}
    			break;
    		}
    	}
    	if(outside){
			sprintf(name,"%d",2000);
			int f=-((j+1)+fiber_N*fiber_N*(T_index-1));
			new G4PVPlacement(0,G4ThreeVector(center_x_BL.at(j),center_y_BL.at(j),0),
    		(bool_cfiber_BL.at(j)==0)?fiberCLog[2*z]:fiberSLog[2*z],name,towerLogicalBL[i],false,f,false);
			}
    }
}

void B4DetectorConstruction::fiberER(G4int i,G4double deltatheta_){
	vector<G4double> temp_x;
    vector<G4double> temp_y;// vector for grid calculation
    
    temp_x.clear();
    temp_y.clear();
    temp_x.push_back(0.*mm);
    temp_y.push_back(0.*mm);
    
    int fiber_N=1500;
    int T_index=NbOfBarrel+i+1;
	
    for(int j = 0 ; j<fiber_N;j++){
    	temp_x.push_back(0.*mm);
    	temp_y.push_back(1.5*(j+1)*mm);
    	temp_x.push_back(0.*mm);
    	temp_y.push_back(-1.5*(j+1)*mm);
    }
    for(int i=0; i<fiber_N;i++){
    	temp_x.push_back(1.5*(i+1)*mm);
    	temp_y.push_back(0.*mm);
    	temp_x.push_back(-1.5*(i+1)*mm);
    	temp_y.push_back(0.*mm);
    }
    for(int j = 0;j<fiber_N;j++){
    	for(int i = 0;i<fiber_N;i++){
    		temp_x.push_back(1.5*(i+1)*mm);
    		temp_y.push_back(1.5*(j+1)*mm);
    		temp_x.push_back(-1.5*(i+1)*mm);
    		temp_y.push_back(-1.5*(j+1)*mm);
    		temp_x.push_back(-1.5*(i+1)*mm);
    		temp_y.push_back(1.5*(j+1)*mm);
    		temp_x.push_back(1.5*(i+1)*mm);
    		temp_y.push_back(-1.5*(j+1)*mm);
    	}
    }
    
    G4double center_x;
    G4double center_y;
    int hi1;
    int hi2;
    int hi1reminder; // 0 -> even 1-> odd
    int hi2reminder;
    
    int numx;
    int numy;
    int reminderx;
    int remindery;
    
    G4ThreeVector v1_ = dimE->GetV3();
    G4ThreeVector v2_ = dimE->GetV4();
    G4ThreeVector v3_ = dimE->GetV1();
    G4ThreeVector v4_ = dimE->GetV2();
    
    G4double innerSide_half_ = dimE->GetInnerR_new()*tan(deltatheta_/2.);
    G4double outerSide_half_ = (dimE->GetInnerR_new()+tower_height)*tan(deltatheta_/2.);
    
    G4double theta_bc=0;
    /////////////////////
    //1. 4 types of grid coordinate(x,y) =  (o,o),(e,o),(o,e),(e,e) { e= 0, o=1 }
    
    
    int type_x_ER;
    int type_y_ER;
    
    numx = (int)(((v4_.getX()*tan(phi_unit/2.)*2)-1.*mm)/1.5*mm);
    numy = (int)((outerSide_half_*2-1.*mm)/(1.5*mm));
    reminderx = numx%2;
    remindery = numy%2;
    
    if(reminderx == 1) type_x_ER=0;
    if(reminderx == 0) type_x_ER=1;
    if(remindery == 1) type_y_ER=0;
    if(remindery == 0) type_y_ER=1;
    
    ////2. aplying the boundary conditions, reject the fibre which near the boundary of tower surface, get coord of cetres of fibre
    
    vector<G4double> center_x_ER;
    vector<G4double> center_y_ER;
    vector<G4int> bool_cfiber_ER;// 0 -> c fiber 1 s fiber
    
    theta_bc=atan(2*outerSide_half_/((v2_.getX()-v4_.getX())*tan(phi_unit/2.)));
    
    if(type_x_ER==1&&type_y_ER==1){//(o,o)
    	for(int num = 0;num<temp_x.size();num++){
    		center_x = temp_x.at(num);
    		center_y = temp_y.at(num);
    		hi1 = fabs(temp_x.at(num))/(1.5*mm);
    		hi2 = fabs(temp_y.at(num))/(1.5*mm);
            hi1reminder = hi1%2; // 0 -> even 1-> odd
            hi2reminder = hi2%2;
            
            if(((center_y-cos(theta_bc)*clad_C_rMax)>(2*2*outerSide_half_/(2*v2_.getX()*tan(phi_unit/2.)-2*v4_.getX()*tan(phi_unit/2.))*(center_x+sin(theta_bc)*clad_C_rMax-v2_.getX()*tan(phi_unit/2.))+outerSide_half_))
            	&&
            	((center_y+clad_C_rMax)<(outerSide_half_))
            	&&
            	((center_y-clad_C_rMax)>(-outerSide_half_))
            	&&
            	((center_y-cos(theta_bc)*clad_C_rMax)>(2*2*outerSide_half_/(2*v4_.getX()*tan(phi_unit/2.)-2*v2_.getX()*tan(phi_unit/2.))*(center_x-sin(theta_bc)*clad_C_rMax+v2_.getX()*tan(phi_unit/2.))+outerSide_half_
            		)))
            {
            	center_x_ER.push_back(center_x);
            	center_y_ER.push_back(center_y);
            	if(((hi1reminder==0)&&(hi2reminder==0))||((hi1reminder==1)&&(hi2reminder==1)))
            		bool_cfiber_ER.push_back(0);
            	else bool_cfiber_ER.push_back(1);
            }
            
        }
    }
    
    if(type_x_ER==0&&type_y_ER==0){//(e,e)

    	for(int num = 0;num<temp_x.size();num++){
    		center_x = temp_x.at(num)+0.75*mm;
    		center_y = temp_y.at(num)-0.75*mm;

    		hi1 = fabs(temp_x.at(num))/(1.5*mm);
    		hi2 = fabs(temp_y.at(num))/(1.5*mm);
            hi1reminder = hi1%2; // 0 -> even 1-> odd
            hi2reminder = hi2%2;
            
            if(((center_y-cos(theta_bc)*clad_C_rMax)>(2*2*outerSide_half_/(2*v2_.getX()*tan(phi_unit/2.)-2*v4_.getX()*tan(phi_unit/2.))*(center_x+sin(theta_bc)*clad_C_rMax-v2_.getX()*tan(phi_unit/2.))+outerSide_half_))
            	&&
            	((center_y+clad_C_rMax)<(outerSide_half_))
            	&&
            	((center_y-clad_C_rMax)>(-outerSide_half_))
            	&&
            	((center_y-cos(theta_bc)*clad_C_rMax)>(2*2*outerSide_half_/(2*v4_.getX()*tan(phi_unit/2.)-2*v2_.getX()*tan(phi_unit/2.))*(center_x-sin(theta_bc)*clad_C_rMax+v2_.getX()*tan(phi_unit/2.))+outerSide_half_
            		)))
            {
            	center_x_ER.push_back(center_x);
            	center_y_ER.push_back(center_y);
            	if(((hi1reminder==0)&&(hi2reminder==0))||((hi1reminder==1)&&(hi2reminder==1)))
            		bool_cfiber_ER.push_back(0);
            	else bool_cfiber_ER.push_back(1);
            }
        }
    }
    
    
    if(type_x_ER==0&&type_y_ER==1){//(e,o)

    	for(int num = 0;num<temp_x.size();num++){
    		center_x = temp_x.at(num)+0.75*mm;
    		center_y = temp_y.at(num);

    		hi1 = fabs(temp_x.at(num))/(1.5*mm);
    		hi2 = fabs(temp_y.at(num))/(1.5*mm);
            hi1reminder = hi1%2; // 0 -> even 1-> odd
            hi2reminder = hi2%2;
            
            if(((center_y-cos(theta_bc)*clad_C_rMax)>(2*2*outerSide_half_/(2*v2_.getX()*tan(phi_unit/2.)-2*v4_.getX()*tan(phi_unit/2.))*(center_x+sin(theta_bc)*clad_C_rMax-v2_.getX()*tan(phi_unit/2.))+outerSide_half_))

            	&&
            	((center_y+clad_C_rMax)<(outerSide_half_))
            	&&
            	((center_y-clad_C_rMax)>(-outerSide_half_))
            	&&
            	((center_y-cos(theta_bc)*clad_C_rMax)>(2*2*outerSide_half_/(2*v4_.getX()*tan(phi_unit/2.)-2*v2_.getX()*tan(phi_unit/2.))*(center_x-sin(theta_bc)*clad_C_rMax+v2_.getX()*tan(phi_unit/2.))+outerSide_half_
            		)))
            {
            	center_x_ER.push_back(center_x);
            	center_y_ER.push_back(center_y);
            	if(((hi1reminder==0)&&(hi2reminder==0))||((hi1reminder==1)&&(hi2reminder==1)))
            		bool_cfiber_ER.push_back(0);
            	else bool_cfiber_ER.push_back(1);
            }
            
        }
    }
    if(type_x_ER==1&&type_y_ER==0){//(o,e)

    	for(int num = 0;num<temp_x.size();num++){
    		center_x = temp_x.at(num);
    		center_y = temp_y.at(num)-0.75*mm;

    		hi1 = fabs(temp_x.at(num))/(1.5*mm);
    		hi2 = fabs(temp_y.at(num))/(1.5*mm);
            hi1reminder = hi1%2; // 0 -> even 1-> odd
            hi2reminder = hi2%2;
            
            if(((center_y-cos(theta_bc)*clad_C_rMax)>(2*2*outerSide_half_/(2*v2_.getX()*tan(phi_unit/2.)-2*v4_.getX()*tan(phi_unit/2.))*(center_x+sin(theta_bc)*clad_C_rMax-v2_.getX()*tan(phi_unit/2.))+outerSide_half_))
            	&&
            	((center_y+clad_C_rMax)<(outerSide_half_))
            	&&
            	((center_y-clad_C_rMax)>(-outerSide_half_))
            	&&
            	((center_y-cos(theta_bc)*clad_C_rMax)>(2*2*outerSide_half_/(2*v4_.getX()*tan(phi_unit/2.)-2*v2_.getX()*tan(phi_unit/2.))*(center_x-sin(theta_bc)*clad_C_rMax+v2_.getX()*tan(phi_unit/2.))+outerSide_half_
            		)))
            {
            	center_x_ER.push_back(center_x);
            	center_y_ER.push_back(center_y);
            	if(((hi1reminder==0)&&(hi2reminder==0))||((hi1reminder==1)&&(hi2reminder==1)))
            		bool_cfiber_ER.push_back(0);
            	else bool_cfiber_ER.push_back(1);
            }
            
        }
    }
    
	//int counter=0;
    for(int j=0;j<center_x_ER.size();j++){
    	double z = tower->GetZHalfLength();
    	bool outside = false;
    	for(int ip = 0;ip<4;ip++){
    		TrapSidePlane plane = tower->GetSidePlane(ip);
    		double zpoint = (-plane.a*center_x_ER.at(j)-plane.b*center_y_ER.at(j)-plane.d)/plane.c;
    		outside = (tower->Inside(G4ThreeVector(center_x_ER.at(j),center_y_ER.at(j),zpoint))==kOutside);
    		if(!outside){
				
    			G4ThreeVector normal = tower->SurfaceNormal(G4ThreeVector(center_x_ER.at(j),center_y_ER.at(j),zpoint));
    			double angle = normal.angle(G4ThreeVector(0,0,-1));
    			double shift = fabs(clad_C_rMax/tan(0.5*M_PI-angle));
    			int length = z - zpoint - shift;
				//cout<<length<<std::endl;
				sprintf(name,"%d",length);
    			if(length>=1&&length<=2*z){
					//counter++; //sprintf(name,"%.1f_%.1f",center_x_ER.at(j),center_y_ER.at(j));
					int f=(j+1)+fiber_N*fiber_N*(T_index-1);
    				new G4PVPlacement(0,G4ThreeVector(center_x_ER.at(j),center_y_ER.at(j),z-0.5*length),
    					(bool_cfiber_ER.at(j)==0)?fiberCLog[length]:fiberSLog[length],name,towerLogicalER[i],false,f,false);
    			}
    			break;
    		}
    	}
    	if(outside){
			//cout<<" 2000 ";
			sprintf(name,"%d",2000);
			int f=(j+1)+fiber_N*fiber_N*(T_index-1);
			
			new G4PVPlacement(0,G4ThreeVector(center_x_ER.at(j),center_y_ER.at(j),0),
    		(bool_cfiber_ER.at(j)==0)?fiberCLog[2*z]:fiberSLog[2*z],name,towerLogicalER[i],false,f,false);
			//counter++; 
		}
    }
	//G4cout<<counter<<std::endl;
}

void B4DetectorConstruction::fiberEL(G4int i,G4double deltatheta_)
{
    vector<G4double> temp_x;
    vector<G4double> temp_y;// vector for grid calculation
    
    temp_x.clear();
    temp_y.clear();
    temp_x.push_back(0.*mm);
    temp_y.push_back(0.*mm);
    
    int T_index=NbOfBarrel+i+1;
    int fiber_N = 1500;
    
    for(int j = 0 ; j<fiber_N;j++){
    	temp_x.push_back(0.*mm);
    	temp_y.push_back(1.5*(j+1)*mm);
    	temp_x.push_back(0.*mm);
    	temp_y.push_back(-1.5*(j+1)*mm);
    }
    for(int i=0; i<fiber_N;i++){
    	temp_x.push_back(1.5*(i+1)*mm);
    	temp_y.push_back(0.*mm);
    	temp_x.push_back(-1.5*(i+1)*mm);
    	temp_y.push_back(0.*mm);
    }
    for(int j = 0;j<fiber_N;j++){
    	for(int i = 0;i<fiber_N;i++){
    		temp_x.push_back(1.5*(i+1)*mm);
    		temp_y.push_back(1.5*(j+1)*mm);
    		temp_x.push_back(-1.5*(i+1)*mm);
    		temp_y.push_back(-1.5*(j+1)*mm);
    		temp_x.push_back(-1.5*(i+1)*mm);
    		temp_y.push_back(1.5*(j+1)*mm);
    		temp_x.push_back(1.5*(i+1)*mm);
    		temp_y.push_back(-1.5*(j+1)*mm);
    	}
    }
    
    G4double center_x;
    G4double center_y;
    int hi1;
    int hi2;
    int hi1reminder; // 0 -> even 1-> odd
    int hi2reminder;
    
    int numx;
    int numy;
    int reminderx;
    int remindery;
    
    G4ThreeVector v1_ = dimE->GetV3();
    G4ThreeVector v2_ = dimE->GetV4();
    G4ThreeVector v3_ = dimE->GetV1();
    G4ThreeVector v4_ = dimE->GetV2();
    G4double innerSide_half_ = dimE->GetInnerR_new()*tan(deltatheta_/2.);
    G4double outerSide_half_ = (dimE->GetInnerR_new()+tower_height)*tan(deltatheta_/2.);
    
    G4double theta_bc=0;
    /////////////////////
    //1. 4 types of grid coordinate(x,y) =  (o,o),(e,o),(o,e),(e,e) { e= 0, o=1 }
    
    int type_x_EL;
    int type_y_EL;
    
    numx = (int)(((v4_.getX()*tan(phi_unit/2.)*2)-1.*mm)/1.5*mm);
    numy = (int)((outerSide_half_*2-1.*mm)/(1.5*mm));
    reminderx = numx%2;
    remindery = numy%2;
    
    if(reminderx == 1) type_x_EL=0;
    if(reminderx == 0) type_x_EL=1;
    if(remindery == 1) type_y_EL=0;
    if(remindery == 0) type_y_EL=1;
    
    ////2. aplying the boundary conditions, reject the fibre which near the boundary of tower surface, get coord of centres of fibre
    
    
    vector<G4double> center_x_EL;
    vector<G4double> center_y_EL;
    vector<G4int> bool_cfiber_EL;// 0 -> c fiber 1 s fiber
    
    theta_bc=atan(2*outerSide_half_/((v2_.getX()-v4_.getX())*tan(phi_unit/2.)));
    
    if(type_x_EL==1&&type_y_EL==1){//(o,o)
    	for(int num = 0;num<temp_x.size();num++){
    		center_x = temp_x.at(num);
    		center_y = temp_y.at(num);
    		hi1 = fabs(temp_x.at(num))/(1.5*mm);
    		hi2 = fabs(temp_y.at(num))/(1.5*mm);
            hi1reminder = hi1%2; // 0 -> even 1-> odd
            hi2reminder = hi2%2;
            
            if(((center_y+cos(theta_bc)*clad_C_rMax)<(2*2*outerSide_half_/(2*v4_.getX()*tan(phi_unit/2.)-2*v2_.getX()*tan(phi_unit/2.))*(center_x+sin(theta_bc)*clad_C_rMax-v4_.getX()*tan(phi_unit/2.))+outerSide_half_))
            	&&
            	((center_y+clad_C_rMax)<(outerSide_half_))
            	&&
            	((center_y-clad_C_rMax)>(-outerSide_half_))
            	&&
            	((center_y+cos(theta_bc)*clad_C_rMax)<(2*2*outerSide_half_/(2*v2_.getX()*tan(phi_unit/2.)-2*v4_.getX()*tan(phi_unit/2.))*(center_x-sin(theta_bc)*clad_C_rMax+v4_.getX()*tan(phi_unit/2.))+outerSide_half_
            		)))
            {
            	center_x_EL.push_back(-center_x);
            	center_y_EL.push_back(-center_y);
            	if(((hi1reminder==0)&&(hi2reminder==0))||((hi1reminder==1)&&(hi2reminder==1)))
            		bool_cfiber_EL.push_back(0);
            	else bool_cfiber_EL.push_back(1);
            }
            
        }
    }
    if(type_x_EL==0&&type_y_EL==0){//(e,e)

    	for(int num = 0;num<temp_x.size();num++){
    		center_x = temp_x.at(num)+0.75*mm;
    		center_y = temp_y.at(num)-0.75*mm;

    		hi1 = fabs(temp_x.at(num))/(1.5*mm);
    		hi2 = fabs(temp_y.at(num))/(1.5*mm);
            hi1reminder = hi1%2; // 0 -> even 1-> odd
            hi2reminder = hi2%2;
            
            if(((center_y+cos(theta_bc)*clad_C_rMax)<(2*2*outerSide_half_/(2*v4_.getX()*tan(phi_unit/2.)-2*v2_.getX()*tan(phi_unit/2.))*(center_x+sin(theta_bc)*clad_C_rMax-v4_.getX()*tan(phi_unit/2.))+outerSide_half_))
            	&&
            	((center_y+clad_C_rMax)<(outerSide_half_))
            	&&
            	((center_y-clad_C_rMax)>(-outerSide_half_))
            	&&
            	((center_y+cos(theta_bc)*clad_C_rMax)<(2*2*outerSide_half_/(2*v2_.getX()*tan(phi_unit/2.)-2*v4_.getX()*tan(phi_unit/2.))*(center_x-sin(theta_bc)*clad_C_rMax+v4_.getX()*tan(phi_unit/2.))+outerSide_half_
            		)))
            {
            	center_x_EL.push_back(-center_x);
            	center_y_EL.push_back(-center_y);
            	if(((hi1reminder==0)&&(hi2reminder==0))||((hi1reminder==1)&&(hi2reminder==1)))
            		bool_cfiber_EL.push_back(0);
            	else bool_cfiber_EL.push_back(1);
            }
        }
    }
    if(type_x_EL==0&&type_y_EL==1){//(e,o)

    	for(int num = 0;num<temp_x.size();num++){
    		center_x = temp_x.at(num)+0.75*mm;
    		center_y = temp_y.at(num);

    		hi1 = fabs(temp_x.at(num))/(1.5*mm);
    		hi2 = fabs(temp_y.at(num))/(1.5*mm);
            hi1reminder = hi1%2; // 0 -> even 1-> odd
            hi2reminder = hi2%2;
            
            if(((center_y+cos(theta_bc)*clad_C_rMax)<(2*2*outerSide_half_/(2*v4_.getX()*tan(phi_unit/2.)-2*v2_.getX()*tan(phi_unit/2.))*(center_x+sin(theta_bc)*clad_C_rMax-v4_.getX()*tan(phi_unit/2.))+outerSide_half_))
            	&&
            	((center_y+clad_C_rMax)<(outerSide_half_))
            	&&
            	((center_y-clad_C_rMax)>(-outerSide_half_))
            	&&
            	((center_y+cos(theta_bc)*clad_C_rMax)<(2*2*outerSide_half_/(2*v2_.getX()*tan(phi_unit/2.)-2*v4_.getX()*tan(phi_unit/2.))*(center_x-sin(theta_bc)*clad_C_rMax+v4_.getX()*tan(phi_unit/2.))+outerSide_half_
            		)))
            {
            	center_x_EL.push_back(-center_x);
            	center_y_EL.push_back(-center_y);
            	if(((hi1reminder==0)&&(hi2reminder==0))||((hi1reminder==1)&&(hi2reminder==1)))
            		bool_cfiber_EL.push_back(0);
            	else bool_cfiber_EL.push_back(1);
            }
            
        }
    }
    if(type_x_EL==1&&type_y_EL==0){//(o,e)

    	for(int num = 0;num<temp_x.size();num++){
    		center_x = temp_x.at(num);
    		center_y = temp_y.at(num)-0.75*mm;

    		hi1 = fabs(temp_x.at(num))/(1.5*mm);
    		hi2 = fabs(temp_y.at(num))/(1.5*mm);
            hi1reminder = hi1%2; // 0 -> even 1-> odd
            hi2reminder = hi2%2;
            
            if(((center_y+cos(theta_bc)*clad_C_rMax)<(2*2*outerSide_half_/(2*v4_.getX()*tan(phi_unit/2.)-2*v2_.getX()*tan(phi_unit/2.))*(center_x+sin(theta_bc)*clad_C_rMax-v4_.getX()*tan(phi_unit/2.))+outerSide_half_))
            	&&
            	((center_y+clad_C_rMax)<(outerSide_half_))
            	&&
            	((center_y-clad_C_rMax)>(-outerSide_half_))
            	&&
            	((center_y+cos(theta_bc)*clad_C_rMax)<(2*2*outerSide_half_/(2*v2_.getX()*tan(phi_unit/2.)-2*v4_.getX()*tan(phi_unit/2.))*(center_x-sin(theta_bc)*clad_C_rMax+v4_.getX()*tan(phi_unit/2.))+outerSide_half_
            		)))
            {
            	center_x_EL.push_back(-center_x);
            	center_y_EL.push_back(-center_y);
            	if(((hi1reminder==0)&&(hi2reminder==0))||((hi1reminder==1)&&(hi2reminder==1)))
            		bool_cfiber_EL.push_back(0);
            	else bool_cfiber_EL.push_back(1);
            }
            
        }
    }
    
    for(int j = 0; j<center_x_EL.size();j++){
    	double z = tower->GetZHalfLength();
    	bool outside = false;
    	for(int ip = 0;ip<4;ip++){
    		TrapSidePlane plane = tower->GetSidePlane(ip);
    		double zpoint = (-plane.a*center_x_EL.at(j)-plane.b*center_y_EL.at(j)-plane.d)/plane.c;
    		outside = (tower->Inside(G4ThreeVector(center_x_EL.at(j),center_y_EL.at(j),zpoint))==kOutside);
    		if(!outside){
    			G4ThreeVector normal = tower->SurfaceNormal(G4ThreeVector(center_x_EL.at(j),center_y_EL.at(j),zpoint));
    			double angle = normal.angle(G4ThreeVector(0,0,-1));
    			double shift = fabs(clad_C_rMax/tan(0.5*M_PI-angle));
    			int length = z - zpoint - shift;
				sprintf(name,"%d",length);
    			if(length>=1&&length<=2*z){
					int f=-((j+1)+fiber_N*fiber_N*(T_index-1));
    				new G4PVPlacement(0,G4ThreeVector(center_x_EL.at(j),center_y_EL.at(j),z-0.5*length),
    					(bool_cfiber_EL.at(j)==0)?fiberCLog[length]:fiberSLog[length],name,towerLogicalEL[i],false,f,false);
    			}
    			break;
    		}
    	}
    	if(outside){
			sprintf(name,"%d",2000);
			int f=-((j+1)+fiber_N*fiber_N*(T_index-1));
			new G4PVPlacement(0,G4ThreeVector(center_x_EL.at(j),center_y_EL.at(j),0),
    		(bool_cfiber_EL.at(j)==0)?fiberCLog[2*z]:fiberSLog[2*z],name,towerLogicalEL[i],false,f,false);
			}
}}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BarrelPhiParameterisation::BarrelPhiParameterisation(G4double inner_R, G4double phi_unit, G4double height, G4double offset)
{
    fInnerR  = inner_R;
    fPhiUnit = phi_unit;
    fHeight  = height;
    fOffset  = offset;
}

void BarrelPhiParameterisation::ComputeTransformation (const G4int copyNo, G4VPhysicalVolume* physVol) const
{

    G4ThreeVector origin((fInnerR+0.5*fHeight)*cos((copyNo+fOffset)*fPhiUnit),(fInnerR+0.5*fHeight)*sin((copyNo+fOffset)*fPhiUnit),0);    
    physVol->SetTranslation(origin);
        
    G4RotationMatrix* rmB = new G4RotationMatrix();
    rmB->rotateZ(M_PI/2.);
    rmB->rotateZ(-(copyNo+fOffset)*fPhiUnit);
    rmB->rotateX(M_PI/2.);
    
    physVol->SetRotation(rmB);
}


BarrelTimingPhiParameterisation::BarrelTimingPhiParameterisation(G4double inner_R, G4double phi_unit, G4double height, G4double offset)
{
    fInnerR  = inner_R;
    fPhiUnit = phi_unit;
    fHeight  = height;
    fOffset  = offset;
}

void BarrelTimingPhiParameterisation::ComputeTransformation (const G4int copyNo, G4VPhysicalVolume* physVol) const
{
       
    G4RotationMatrix* rmB = new G4RotationMatrix();
    rmB->rotateZ(M_PI/2.);
    rmB->rotateZ(-(copyNo+fOffset)*fPhiUnit);    
    physVol->SetRotation(rmB);
}



EndcapThetaParameterisation::EndcapThetaParameterisation(G4double inner_R, G4double theta_unit, G4double height, G4double thetaB)
{
    fInnerR    = inner_R;
    fThetaUnit = theta_unit;
    fHeight    = height;
    fThetaB    = thetaB;
}

void EndcapThetaParameterisation::ComputeTransformation (const G4int copyNo, G4VPhysicalVolume* physVol) const
{

    
//     int offset = 100/4*3;
//     int offset = 100/4*3;
//     int offset = 0;
//     G4ThreeVector origin(0, 0, (fInnerR+fHeight)*tan(fThetaB)-fHeight/2);    
//     physVol->SetTranslation(origin);
        
//     G4RotationMatrix* rmB = new G4RotationMatrix();
//     rmB->rotateZ(M_PI/2.);
//     rmB->rotateZ(-(copyNo+offset)*fThetaUnit);
//     rmB->rotateX(M_PI/2.);
//     
//     physVol->SetRotation(rmB);
    
    
}

void EndcapThetaParameterisation::ComputeDimensions (G4Cons& endcapRing, const G4int copyNo, const G4VPhysicalVolume* physVol) const
{
    
    G4double ring_rmin1  = fInnerR*tan(M_PI/4-fThetaUnit*copyNo);
    G4double ring_rmax1  = fInnerR;
    G4double ring_rmin2  = (fInnerR+fHeight)*tan(M_PI/4-fThetaUnit*copyNo);
    G4double ring_rmax2  = fInnerR+fHeight;
    G4double ring_length = fHeight;
    
    endcapRing.SetInnerRadiusMinusZ (ring_rmin1);
    endcapRing.SetOuterRadiusMinusZ (ring_rmax1);
    endcapRing.SetInnerRadiusPlusZ (ring_rmin2);
    endcapRing.SetOuterRadiusPlusZ (ring_rmax2);
    endcapRing.SetZHalfLength (ring_length/2.);
    endcapRing.SetStartPhiAngle (0.);
    endcapRing.SetDeltaPhiAngle (360.);

    
}


