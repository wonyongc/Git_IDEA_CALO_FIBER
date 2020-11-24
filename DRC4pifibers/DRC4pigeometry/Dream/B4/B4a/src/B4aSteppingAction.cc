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
// $Id: B4aSteppingAction.cc 68058 2013-03-13 14:47:43Z gcosmo $
// 
/// \file B4aSteppingAction.cc
/// \brief Implementation of the B4aSteppingAction class

#include "B4aSteppingAction.hh"
#include "B4aEventAction.hh"
#include "B4DetectorConstruction.hh"
#include "G4Material.hh"
#include "G4UnitsTable.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"
#include <stdlib.h>

#include "G4OpBoundaryProcess.hh"
//#include "Fiber_Info.hh"

#include <chrono>
#include <random>

#include "TLorentzVector.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4aSteppingAction::B4aSteppingAction(
                      const B4DetectorConstruction* detectorConstruction,
                      B4aEventAction* eventAction)
  : G4UserSteppingAction(),
    fDetConstruction(detectorConstruction),
    fEventAction(eventAction)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4aSteppingAction::~B4aSteppingAction()
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4aSteppingAction::UserSteppingAction(const G4Step* step)
{
	
  // get volume of the current pre-step
  G4VPhysicalVolume* PreStepVolume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  G4double energydeposited = step->GetTotalEnergyDeposit();
  G4double steplength = step->GetStepLength();
  
  G4StepPoint * thePrePoint  = step->GetPreStepPoint () ;  
  const G4ThreeVector & thePrePosition  = thePrePoint->GetPosition () ;
  G4double global_x = thePrePosition.x()/mm;
  G4double global_y = thePrePosition.y()/mm;
  G4double global_z = thePrePosition.z()/mm;
//   if (!strstr(PreStepVolume->GetName(),"tower"))   std::cout << " :: x = " << global_x << " :: y = " << global_y << " :: z = " << global_z << " --> step is in: " << PreStepVolume->GetName() << std::endl;
  
  //define Birk's constant
  double k_B = 0.126; 
  G4double saturatedenergydeposited = 0.;

  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);

  std::poisson_distribution<int> cher_distribution(1*0.2723);

  std::poisson_distribution<int> scin_distribution(1.);
  
  if (PreStepVolume->GetName() != "World")
  {
    fEventAction->Addenergy(energydeposited);
    if (PreStepVolume->GetLogicalVolume()->GetMaterial()->GetName() == "Copper")
    {
        G4double copynumbertower = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(); 
    	G4double copynumberslice = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(1); 
    	if (copynumbertower > 0){ //im in barrel right or endcap right
            fEventAction->AddVectorR(energydeposited,copynumbertower, copynumberslice);
    	}
        if (copynumbertower < 0){ //im in barrel left or endcap left
            fEventAction->AddVectorL(energydeposited, copynumbertower, copynumberslice);
        }
    }
  }
  
  //looper finding
  if (step->GetTrack()->GetParentID() < 10 && PreStepVolume->GetLogicalVolume()->GetMaterial()->GetName() == "Copper")
  {
  	auto p = step->GetTrack()->GetVertexMomentumDirection();
  	auto e = step->GetTrack()->GetVertexKineticEnergy();
  	auto m = step->GetTrack()->GetDynamicParticle()->GetMass();
  	G4double pmodule = std::sqrt(e*e+2*e*m);
  	G4double pt = std::sqrt(p[0]*pmodule*p[0]*pmodule+p[1]*pmodule*p[1]*pmodule);
  	//G4cout<< p << " " << e <<" " << pt <<G4endl;
  	//G4cout<<step->GetTrack()->GetMomentum()<<" "<<step->GetTrack()->GetKineticEnergy()<<" "<<step->GetTrack()->GetMomentumDirection()<<" "<<step->GetTrack()->GetVertexMomentumDirection()<<" "<<step->GetTrack()->GetDefinition()->GetPDGCharge()<<" "<<step->GetTrack()->GetDefinition()->GetParticleName()<<G4endl;
  	if (pt<630. && step->GetTrack()->GetDefinition()->GetPDGCharge() != 0){	
            //G4cout<<"one looper found! "<<step->GetTrack()->GetVertexKineticEnergy()<<G4endl;
            //G4cout<<step->GetTrack()->GetDefinition()->GetParticleName()<<" "<<step->GetTrack()->GetTrackID()<<" "<<pt<<" "<<step->GetTrack()->GetDefinition()->GetParticleName()<<G4endl;
            G4double copynumbertower = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(); 
            G4double copynumberslice = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(1); 
            if (copynumbertower > 0){ //im in barrel right or endcap right
     		fEventAction->AddVectorR_loop(pt,copynumbertower, copynumberslice);
            }
            if (copynumbertower < 0){ //im in barrel left or endcap left
                fEventAction->AddVectorL_loop(pt, copynumbertower, copynumberslice);
            }
  	}
  }
  //end of looper finding

  //compute leak
  if (PreStepVolume->GetName() == "leakageabsorber"){
    auto name = step->GetTrack()->GetDefinition()->GetParticleName();
    if (name=="nu_mu" || name=="nu_e" || name=="anti_nu_e" || name=="anti_nu_mu"){
      //G4cout<<step->GetTrack()->GetCreatorProcess()->GetProcessName()<<G4endl;
      fEventAction->Addneutrinoleakage(step->GetTrack()->GetKineticEnergy());
      step->GetTrack()->SetTrackStatus(fStopAndKill);}
    else{
      fEventAction->Addleakage(step->GetTrack()->GetKineticEnergy());
      step->GetTrack()->SetTrackStatus(fStopAndKill);}
  }
  //end compute leak

  //cm fraction
  if (PreStepVolume->GetName() != "World" ) {
    if (step->GetTrack()->GetDefinition()->GetParticleName() == "e-" || step->GetTrack()->GetDefinition()->GetParticleName() == "e+"){
      //Function to add up energy deposited by em component
      fEventAction->Addem(energydeposited);
    }
  }
  //end of em fraction
 
 //primary particle energy
//    if ( step->GetTrack()->GetTrackID() == 1 && step->GetTrack()->GetCurrentStepNumber() == 1){
//     // Function to save primary particle energy and name
//     fEventAction->SavePrimaryParticle(step->GetTrack()->GetDefinition()->GetParticleName());
//     fEventAction->SavePrimaryEnergy(step->GetTrack()->GetVertexKineticEnergy());
//     fEventAction->SavePrimaryMomentum(step->GetTrack()->GetVertexKineticEnergy());
//   }
  //end of primary particle energy
		
		

  
  
  //Here I compute and save all informations about scintillating and Cherenkov fibers
  std::string Fiber;
  std::string S_fiber = "fiberCoreScint";
  std::string C_fiber = "fiberCoreChere";
  Fiber = PreStepVolume->GetName(); //name of current step fiber
  
  if ( strstr(Fiber.c_str(),S_fiber.c_str())){ //it's a scintillating fiber
    //Function to add up energy deposited in scintillating fibers:
    //- as signal saturated by Birk's law in VectorSignals
    //- as regular energy deposition in all scintillating fibers in EnergyScin

    if(step->GetTrack()->GetDefinition()->GetPDGCharge() != 0.)
    {
        if (steplength != 0)
                {
                    saturatedenergydeposited = (energydeposited/steplength) / ( 1+k_B*(energydeposited/steplength) ) * steplength;
                }
    }
    else
    {
    	saturatedenergydeposited = 0.;
    }
	  
    fEventAction->AddScin(energydeposited); //All energy deposited in scin fibers (not saturated)
    G4double copynumbertower = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(2); 
    G4double copynumberslice = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(3); 
    G4int Sfibercopynumber = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(1);
	  
    std::string LengthFibr =  step->GetPreStepPoint()->GetTouchableHandle()->GetVolume(1)->GetName(); 
	  
    G4double S_fiber_ID = 0;

    scin_distribution = std::poisson_distribution<int> (saturatedenergydeposited*12.5);
    int s_signal = scin_distribution(generator);
	
    if (copynumbertower > 0){ //im in barrel right or endcap right
        fEventAction->AddVectorScinEnergyR(s_signal,copynumbertower, copynumberslice); //energy deposited in any scintillating fiber (saturated)
        fEventAction->AddVectorR(energydeposited, copynumbertower, copynumberslice);
	 //I want unique Fiber ID: 168750000 is the max of Sfibercopynumber
	 S_fiber_ID = Sfibercopynumber+(168750000*copynumberslice);
    }  
	  
    if (copynumbertower < 0){ //im in barrel left or endcap left
        fEventAction->AddVectorScinEnergyL(s_signal, copynumbertower, copynumberslice);
        fEventAction->AddVectorL(energydeposited, copynumbertower, copynumberslice);
        //I want unique Fiber ID: 168750000 is the max of Sfibercopynumber
        S_fiber_ID = Sfibercopynumber-(168750000*copynumberslice);
    }
	
	// Fibers routine: fill the S fibres info 
    if (saturatedenergydeposited>0.)
    {
        G4VPhysicalVolume* physVol=step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
		
        //G4ThreeVector vectPos(0,0,0);
        G4int k=0;
		
        //LOCAL TO GLOBAL TRANSFORMATIONS
        G4TouchableHandle theTouchable = step->GetPreStepPoint()->GetTouchableHandle();
        G4ThreeVector origin(0.,0.,0.);
        G4ThreeVector zdir(0.,0.,1.);
        G4ThreeVector vectPos = theTouchable->GetHistory()->
    	GetTopTransform().Inverse().TransformPoint(origin);
        G4ThreeVector direction = theTouchable->GetHistory()->
    	GetTopTransform().Inverse().TransformAxis(zdir);
        G4double lengthfiber = atof(LengthFibr.c_str());
        G4ThreeVector Halffibervect = direction*lengthfiber/2;
        // Fibre tip position
        G4ThreeVector vectPostip = vectPos-Halffibervect;
        // SiPM position
        G4ThreeVector SiPMvecPos = vectPos+Halffibervect;
        if (s_signal>0.0)
        {
            fEventAction->WriteFiber_Info(S_fiber_ID,s_signal,1,vectPostip,copynumberslice,copynumbertower);// 1 == S 0 == C
            // Extract info for z time
            std::ofstream TimeFile;
            TimeFile.open("Time.txt", std::ios_base::app);
            TimeFile<<"Scin "<< std::fixed << std::setprecision(3) <<S_fiber_ID<<" "<<vectPostip.getX()<<" "<<vectPostip.getY()<<" "<<vectPostip.getZ()<<" "<<s_signal<<" "<<sqrt((SiPMvecPos[0]-step->GetTrack()->GetPosition().getX())*(SiPMvecPos[0]-step->GetTrack()->GetPosition().getX())+(SiPMvecPos[1]-step->GetTrack()->GetPosition().getY())*(SiPMvecPos[1]-step->GetTrack()->GetPosition().getY())+(SiPMvecPos[2]-step->GetTrack()->GetPosition().getZ())*(SiPMvecPos[2]-step->GetTrack()->GetPosition().getZ()))<<" "<<step->GetTrack()->GetGlobalTime()<<G4endl;
		//TimeFile.close();
        }
    }
  }

  if ( strstr(Fiber.c_str(),C_fiber.c_str())){//it's a Cherenkov fiber
    //Function to add up energy deposited in Cherenkov fibres
	fEventAction->AddCher(step->GetTotalEnergyDeposit());
  
  	G4double copynumbertower = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(2); 
  	G4double copynumberslice = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(3);

  	if (copynumbertower>0) {fEventAction->AddVectorR(energydeposited, copynumbertower, copynumberslice);}
  	if (copynumbertower<0) {fEventAction->AddVectorL(energydeposited, copynumbertower, copynumberslice);}
  	}
 //part for cherenkov photons
 G4OpBoundaryProcessStatus theStatus = Undefined;

 G4ProcessManager* OpManager =
                    G4OpticalPhoton::OpticalPhoton()->GetProcessManager();

 if (OpManager) {
     G4int MAXofPostStepLoops =
              OpManager->GetPostStepProcessVector()->entries();
     G4ProcessVector* fPostStepDoItVector =
              OpManager->GetPostStepProcessVector(typeDoIt);

     for ( G4int i=0; i<MAXofPostStepLoops; i++) {
         G4VProcess* fCurrentProcess = (*fPostStepDoItVector)[i];
         fOpProcess = dynamic_cast<G4OpBoundaryProcess*>(fCurrentProcess);
         if (fOpProcess) { theStatus = fOpProcess->GetStatus(); break;}
     }
  }

  
  
  
  //***************************************************************
  //                    SCEPCal ECAL scoring
  //***************************************************************
  
  
  if (PreStepVolume->GetName().contains("crystalECAL"))
  {          
    int nBarrelEtaSeg = 180;
    int nBarrelPhiSeg = 1130;
    
    int nEndcapRings  = 162;
    int nEndcapMaxSeg = 216;
    
    int chOffset = 1000000; //to make sure endcap numbering does not overlap with barrel
    int crystal_0IDR[162] = {0, 1121, 2232, 3333, 4425, 5507, 6580, 7643, 8697, 9742, 10778, 11805, 12823, 13832, 14832, 15823, 16806, 17780, 18745, 19702, 20650, 21590, 22522, 23446, 24361, 25268, 26167, 27058, 27941, 28816, 29683, 30543, 31395, 32239, 33075, 33904, 34725, 35539, 36345, 37144, 37935, 38719, 39496, 40266, 41028, 41783, 42531, 43272, 44006, 44733, 45453, 46166, 46872, 47571, 48264, 48950, 49629, 50301, 50967, 51626, 52278, 52924, 53563, 54196, 54822, 55442, 56056, 56663, 57264, 57859, 58447, 59029, 59605, 60175, 60738, 61295, 61846, 62391, 62930, 63463, 63990, 64511, 65026, 65535, 66038, 66535, 67026, 67511, 67991, 68465, 68933, 69395, 69851, 70302, 70747, 71186, 71620, 72048, 72470, 72887, 73298, 73704, 74104, 74498, 74887, 75270, 75648, 76020, 76387, 76749, 77105, 77456, 77801, 78141, 78476, 78805, 79129, 79447, 79760, 80068, 80371, 80668, 80960, 81247, 81528, 81804, 82075, 82341, 82602, 82857, 83107, 83352, 83592, 83827, 84057, 84281, 84500, 84714, 84923, 85127, 85326, 85520, 85709, 85893, 86072, 86246, 86415, 86578, 86736, 86889, 87037, 87180, 87318, 87451, 87579, 87702, 87820, 87933, 88041, 88144, 88242, 88335};

    G4double chId = 0.;        
    G4String logVolumeName = PreStepVolume->GetLogicalVolume()->GetName();
    
    //in ECAL barrel
    if (PreStepVolume->GetName().contains("Barrel"))
    {
        //number of crystal within phi slice
        G4double crystal_copy_ID = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(0);         
        //number of phi slice within barrel 
        G4double phi_slice_ID = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(1);                 
        //define crystal absolute ID number                        
        G4int sign    = crystal_copy_ID/fabs(crystal_copy_ID);
        chId = sign*(fabs(crystal_copy_ID) + nBarrelEtaSeg*phi_slice_ID );     
    }
    
    //in ECAL endcap
    else if (PreStepVolume->GetName().contains("Endcap"))
    {        
        //number of crystal within eta ring
        G4double crystal_copy_ID = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(0);         
        //number of phi slice within barrel 
        G4double eta_ring_ID = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(1);                 
        //define crystal absolute ID number                        
        G4int sign    = eta_ring_ID/fabs(eta_ring_ID);
        chId = sign*(crystal_copy_ID + crystal_0IDR[int(abs(eta_ring_ID))-1] + chOffset);             
    }

    //count energy deposited through ionization
    if (energydeposited>0)
    {
        if (PreStepVolume->GetName().contains("Front")) 
        {
            fEventAction->AddScepEneF(energydeposited);
            fEventAction->AddScepHit(chId, energydeposited, "FrontEne");
        }
        if (PreStepVolume->GetName().contains("Rear") ) 
        {
            fEventAction->AddScepEneR(energydeposited); 
            fEventAction->AddScepHit(chId, energydeposited, "RearEne");
        }    
    }
        
    // count cherenkov photons
    G4Track* theTrack = step->GetTrack () ;
    G4int nStep = theTrack -> GetCurrentStepNumber();    
    
    if( (nStep == 1) && (theTrack->GetDefinition()->GetParticleName() == "opticalphoton") )
    {
        G4String processName = theTrack->GetCreatorProcess()->GetProcessName();
        
        if (processName == "Cerenkov")
        {
            //kill very long or short wavelengths
            float photWL = MyMaterials::fromEvToNm(theTrack->GetTotalEnergy()/eV);                
            if (photWL> 1000 ||  photWL< 300)  theTrack->SetTrackStatus(fKillTrackAndSecondaries); 
            
            else 
            {
                if (PreStepVolume->GetName().contains("Front")) 
                {
                    fEventAction->AddScepCherF(); 
                    fEventAction->AddScepHit(chId, 1, "FrontCher");
                }
                if (PreStepVolume->GetName().contains("Rear") ) 
                {
                    fEventAction->AddScepCherR(); 
                    fEventAction->AddScepHit(chId, 1, "RearCher");
                }
                //do not propagate the photon
                theTrack->SetTrackStatus(fKillTrackAndSecondaries); 
            }
        }
    }
  }
  
  
  
  //***************************************************************
  //                    SCEPCal Timing scoring
  //***************************************************************
  
  if (PreStepVolume->GetName().contains("crystalTiming"))
  {          
      
    int nBarrelTiming_Z = 29;   //number of modules along z inside a half barrel phi slice (tray)        
    int nBars = 20;
    int chOffset = 1000000; //to make sure endcap numbering does not overlap with barrel
    
    G4double chId = 0.;        
    G4String logVolumeName = PreStepVolume->GetLogicalVolume()->GetName();    
    // bar ID in the array, from to 0 to 19 (bars)
    G4double crystal_copy_ID = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(0);          
    // front layer is 1 , rear layer is 2
    G4double layer_ID = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(1);          
    // module number
    G4double module_ID = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(2);              
    // mother envelope volume name (e.g. was module inside endcap or inside barrel?)
    G4String motherEnv_volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume(2)->GetName();
//     std::cout << "motherEnv_volume = "  << motherEnv_volume << std::endl;
    
    //in Timing barrel
    if (motherEnv_volume.contains("Barrel"))
    {
        //number of phi slice within barrel               
        G4double phi_slice_ID = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(3);                         
        //define crystal absolute ID number                        
        G4int barrel_sign    = module_ID/fabs(module_ID);
        chId = barrel_sign*(crystal_copy_ID + nBars*fabs(module_ID) + nBars*nBarrelTiming_Z*phi_slice_ID );     
//         std::cout << " hit in timing barrel at chId = " << chId << std::endl;
        
    }
    
    //in Timing endcap
    if (motherEnv_volume.contains("Endcap"))
    {        
        //endcap sign
        G4double endcap_sign = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(3);                         
        //define crystal absolute ID number                                
        chId = endcap_sign*(crystal_copy_ID + nBars*fabs(module_ID) + chOffset);             
//         std::cout << " hit in timing endcap at chId = " << chId << std::endl;
    }

    //count energy deposited through ionization
    if (energydeposited>0)
    {            
        G4double time_hit = step->GetPreStepPoint()->GetGlobalTime()/picosecond;
        fEventAction->AddScepTimingEneTime(energydeposited, time_hit, layer_ID);
        fEventAction->AddScepTimingHit(chId, energydeposited, time_hit, layer_ID);        
//         std::cout << "layer: " << layer_ID << " :: ene = " << energydeposited << " :: time " << time_hit << std::endl;
    }
  }
  
//   if ( (step->GetTrack()->GetDefinition()->GetParticleName() == "opticalphoton")) std::cout << "volume = " << PreStepVolume->GetName() << " :: particle = " << step->GetTrack()->GetDefinition()->GetParticleName() << std::endl;
  
  
  
  //IDEA optical photons      
  
  
  std::string SiPMC = "SiPMC";
  std::string SiPMS = "SiPMS";
  std::string SiPMdetection;
  

  //If the particle is an optical photon...
  if(step->GetTrack()->GetDefinition()->GetParticleName() == "opticalphoton"){

     switch (theStatus){

        case TotalInternalReflection: 
           Fiber = PreStepVolume->GetName();

            if(strstr(Fiber.c_str(), C_fiber.c_str()))
            { //it's a Cherenkov fibre
                G4double copynumbertower = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(2); 
                G4double copynumberslice = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(3); 
                G4int Cfibercopynumber = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(1);
                
                std::string LengthFibr =  step->GetPreStepPoint()->GetTouchableHandle()->GetVolume(1)->GetName(); 
                        
                G4double C_fiber_ID = 0;
                int c_signal = cher_distribution(generator);
                if (copynumbertower>0){ //i'm in barrel right or endcap right
                    fEventAction->AddVectorCherPER(c_signal, copynumbertower, copynumberslice);
                    //I want unique Fiber ID: 168750000 is the max of Cfibercopynumber
                    C_fiber_ID = Cfibercopynumber+(168750000*copynumberslice);
                }
                if (copynumbertower<0){ //i'm in barrel left ot endcap left 
                    fEventAction->AddVectorCherPEL(c_signal, copynumbertower, copynumberslice);
                    //I want unique Fiber ID: 168750000 is the max of Cfibercopynumber
                    C_fiber_ID = Cfibercopynumber-(168750000*copynumberslice);
                }
                fEventAction->AddCherenkov(); // add one photoelectron from Cherenkov process in Cherenkov fibers                  
				
                // Fibers routine: fill the C fibres info 
                G4VPhysicalVolume* physVol=step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
		
                G4int k=0;
                G4TouchableHandle theTouchable = step->GetPreStepPoint()->GetTouchableHandle();
                G4ThreeVector origin(0.,0.,0.);
                G4ThreeVector zdir(0.,0.,1.);
                G4ThreeVector vectPos = theTouchable->GetHistory()->
                GetTopTransform().Inverse().TransformPoint(origin);
                G4ThreeVector direction = theTouchable->GetHistory()->
                GetTopTransform().Inverse().TransformAxis(zdir);
                G4double lengthfiber = atof(LengthFibr.c_str());
                G4ThreeVector Halffibervect = direction*lengthfiber/2;
                // Fibre tip position
                G4ThreeVector vectPostip = vectPos-Halffibervect;
                // SiPM position
                G4ThreeVector SiPMvecPos = vectPos+Halffibervect;
                if (c_signal>0){
                    fEventAction->WriteFiber_Info(C_fiber_ID,c_signal,0,vectPostip,copynumberslice,copynumbertower);// 1 == S 0 == C
                    // Extract info for z time
                    std::ofstream TimeFile;
                    TimeFile.open("Time.txt", std::ios_base::app);
                    TimeFile<<"Cher "<<std::fixed << std::setprecision(3) <<C_fiber_ID<<" "<<vectPostip.getX()<<" "<<vectPostip.getY()<<" "<<vectPostip.getZ()<<" "<<c_signal<<" "<<sqrt((SiPMvecPos[0]-step->GetTrack()->GetPosition().getX())*(SiPMvecPos[0]-step->GetTrack()->GetPosition().getX())+(SiPMvecPos[1]-step->GetTrack()->GetPosition().getY())*(SiPMvecPos[1]-step->GetTrack()->GetPosition().getY())+(SiPMvecPos[2]-step->GetTrack()->GetPosition().getZ())*(SiPMvecPos[2]-step->GetTrack()->GetPosition().getZ()))<<" "<<step->GetTrack()->GetGlobalTime()<<G4endl;
                    //TimeFile.close();
                }
                step->GetTrack()->SetTrackStatus(fStopAndKill); //kill photon
				
            }
            break;
            
            
               /*
               Prestep = step->GetPreStepPoint()->GetPosition();   
               Postsep = step->GetPostStepPoint()->GetPosition();
               Momentum = step->GetTrack()->GetMomentumDirection();
               if(Momentum.z()>0.){ //the photon is going towards SiPms
                costheta = Momentum.z();
                if(costheta>0.99){//0.94*/ //if the photon is under the acceptance angle of fibers
                  /* only if you want exponential light attenuation
                  distance = (1560.9-Prestep.z())/costheta;
                  pSurvive = std::exp(-(distance/8900));
                  pTot=PSurvive*pDetection;*/
/*                  pTot =pDetection;
                  if(pRandom<pTot){  
                    fEventAction->AddCherenkov(); // add one photoelectron from Cherenkov process in Cherenkov fibers                  
                    //fEventAction->AddSignalfibre(copynumber); //only if you want SignalFibre
                    fEventAction->AddVectorCherPE(copynumbermodule,copynumber);
                    step->GetTrack()->SetTrackStatus(fStopAndKill); //I kille the photon just after having counted it or excluded
/*                  }
                }
              }
             }
    break;

  case Detection:
  // if you want no parameterization and complete full simulation uncomment this part
    
   /* SiPMdetection = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName();
    if (strstr(SiPMdetection.c_str(),SiPMC.c_str()))
     {
       fEventAction->AddCherenkov();
     } 
   
    if (strstr(SiPMdetection.c_str(),SiPMS.c_str()))
    {
      fEventAction->AddScintillation();
    }
 
  break;*/

        default: 
        //only for parameterization, comment for full simulation
            step->GetTrack()->SetTrackStatus(fStopAndKill);
            break;
    }
  }
  
  



}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......












