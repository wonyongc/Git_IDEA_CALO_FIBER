#include "SCEPCal_GeometryHelper.hh"

#include <iostream>

SCEPCal_GeometryHelper::SCEPCal_GeometryHelper():
// Dual readout HCAL
  m_NbOfBarrel(40),
  m_NbOfEndCap(35),
  m_NZrot(36),
  //SCEPCal
  m_innerR(1800),
  m_nBarrelEtaSeg(180),
  m_nBarrelPhiSeg(1130),
  m_chOffset (1000000),
  m_nEndcapRings(162),
  m_nBarrelTiming_Z(29),
  m_nBarrelTiming_PhiSeg(186),
  m_BarrelTiming_IR(1775),
  m_BarrelTiming_OR(1795),
  m_barLength(60),
  m_barWidth(3),
  m_nBars(20),  
  m_nEndcapModulePerLine(59),
  m_endOIDR{0, 1121, 2232, 3333, 4425, 5507, 6580, 7643, 8697, 9742, 10778, 11805, 12823, 13832, 14832, 15823, 16806, 17780, 18745, 19702, 20650, 21590, 22522, 23446, 24361, 25268, 26167, 27058, 27941, 28816, 29683, 30543, 31395, 32239, 33075, 33904, 34725, 35539, 36345, 37144, 37935, 38719, 39496, 40266, 41028, 41783, 42531, 43272, 44006, 44733, 45453, 46166, 46872, 47571, 48264, 48950, 49629, 50301, 50967, 51626, 52278, 52924, 53563, 54196, 54822, 55442, 56056, 56663, 57264, 57859, 58447, 59029, 59605, 60175, 60738, 61295, 61846, 62391, 62930, 63463, 63990, 64511, 65026, 65535, 66038, 66535, 67026, 67511, 67991, 68465, 68933, 69395, 69851, 70302, 70747, 71186, 71620, 72048, 72470, 72887, 73298, 73704, 74104, 74498, 74887, 75270, 75648, 76020, 76387, 76749, 77105, 77456, 77801, 78141, 78476, 78805, 79129, 79447, 79760, 80068, 80371, 80668, 80960, 81247, 81528, 81804, 82075, 82341, 82602, 82857, 83107, 83352, 83592, 83827, 84057, 84281, 84500, 84714, 84923, 85127, 85326, 85520, 85709, 85893, 86072, 86246, 86415, 86578, 86736, 86889, 87037, 87180, 87318, 87451, 87579, 87702, 87820, 87933, 88041, 88144, 88242, 88335}
{ 
  // Dual readout HCAL
  m_totTower = m_NbOfBarrel + m_NbOfEndCap;
  m_deltaTheta = 45./(m_NbOfBarrel);                              
  m_phiUnit = 360./m_NZrot;
  //SCEPCal
  m_deltaThetaBarrel    = 45./m_nBarrelEtaSeg;
  m_phiUnitBarrel       = 360./m_nBarrelPhiSeg;
  m_deltaThetaEndcap    = 45./m_nBarrelEtaSeg;
  m_phiUnitTimingBarrel = 360./m_nBarrelTiming_PhiSeg;
}
  
// Dual readout HCAL
TVector3 SCEPCal_GeometryHelper::GetTowerVec(unsigned int index, char side)
{  
  // Create an empty TLorentzVector

  TVector3 tower;
  if (side != 'l' && side != 'r') return tower;
  if (index == 0) return tower;

  --index;

  unsigned int sliceindex = index/m_totTower;
  unsigned int towerindex = index-(sliceindex*m_totTower);
  double theta = towerindex*m_deltaTheta+m_deltaTheta/2.;
  double phi = ((double)sliceindex) * m_phiUnit;

  if (side == 'r') theta = theta + 90.;
  else if (side == 'l') theta = 90. - theta;  
  
  tower.SetMagThetaPhi(1,TMath::DegToRad()*(theta),TMath::DegToRad()*phi);

  return tower;
}


  
// ECAL
TVector3 SCEPCal_GeometryHelper::GetCrystalVec(long int index)
{          
    
  // Create an empty TLorentzVector
  TVector3 crystal;
  
  if (index == 0 || fabs(index) > 2000000) return crystal;
    
  int sign = index/fabs(index);  //get index sign (to distinguish left and right sides)  
  if (fabs(index)<=m_chOffset)       // index in barrel
  {
//       std::cout << "in BARREL!" << std::endl;
      int phi_slice_ID   = floor(fabs(index)/m_nBarrelEtaSeg);
      int theta_id       = fabs(index) - m_nBarrelEtaSeg*(phi_slice_ID);
      
      double theta = 90 + sign*(theta_id-0.5)*m_deltaThetaBarrel;
      if (theta_id == 0) 
      {
          phi_slice_ID   = floor(fabs(index)/m_nBarrelEtaSeg)-1;
          theta_id       = fabs(index) - m_nBarrelEtaSeg*(phi_slice_ID);
          theta = 90 + sign*(theta_id-0.5)*m_deltaThetaBarrel;  
      }
      
      double phi = ((double)phi_slice_ID) * m_phiUnitBarrel;      
      crystal.SetMagThetaPhi(1,TMath::DegToRad()*theta,TMath::DegToRad()*phi);       
  }
  
  else if (fabs(index)>=m_chOffset)  //hit in endcap
  {
      
//       std::cout << "in ENDCAP!" << std::endl;
      int ring_ID  = 0;
      int phi_ID   = 0 ;
      
      double theta = 0.;
      double phi   = 0.;
      double thisRingRadius = 0.;
      double centerRingRadius = 0.;
      
      int end_index = fabs(index) - m_chOffset;
      
      for (int iRing = 0; iRing < m_nEndcapRings; iRing++)
      {
        if (iRing < m_nEndcapRings && end_index >= m_endOIDR[iRing] && end_index< m_endOIDR[iRing+1])
        {
            ring_ID = iRing;            
            phi_ID  = end_index - m_endOIDR[iRing];                        
            theta = 45 + (ring_ID+1.)*m_deltaThetaEndcap;
            thisRingRadius = m_innerR/tan(theta*TMath::DegToRad()); //min radius used to computed nDivision in DetectorConstruction.cc
            centerRingRadius = m_innerR/tan((theta-0.5*m_deltaThetaEndcap)*TMath::DegToRad()); //+SCEP_deltatheta_endcap[iRing]/2);
                        
            double crystalSize = 10.;
            int nDivisions = 2.*M_PI*thisRingRadius/crystalSize;
            double thisPhiUnit = 360./nDivisions;
            phi   = double(phi_ID+0.5)*thisPhiUnit;
//             if (phi == 180)          phi = 360-180.
                
//             std::cout << "ring_ID = " << ring_ID << " :: phi_ID = " << phi_ID << "/" << nDivisions << " (exp div = " << m_endOIDR[iRing+1]-m_endOIDR[iRing] << " ) :: phi " << phi << " deg (" << phi*TMath::DegToRad() << " rad) :: thisRingRadius = " << thisRingRadius << " :: thisPhiUnit = " << thisPhiUnit << " :: index = " << index << std::endl;            
            break;
        }
      }
      
      crystal.SetXYZ(centerRingRadius*cos(phi*TMath::DegToRad()), centerRingRadius*sin(phi*TMath::DegToRad()), -sign*m_innerR);
  }  
  

  return crystal;
}



//*******************************************************//
//                          Timing
//******************************************************//

TVector3 SCEPCal_GeometryHelper::GetCrystalTimingVec(long int index, int layer_ID)
{                    
    
  TVector3 crystal;      
  int sign = index/fabs(index);  //get index sign (to distinguish left and right sides)
  double phi;
  double z_pos, x_pos, y_pos;
  double timing_radius = (m_BarrelTiming_IR+m_BarrelTiming_OR)/2.;
  
  if (fabs(index)<m_chOffset)       // index in barrel
  {
//       std::cout << "in BARREL!" << std::endl;
      int phi_slice_ID = floor(fabs(index)/(m_nBarrelTiming_Z*m_nBars));      
      int module_id    = floor((fabs(index) - m_nBarrelTiming_Z*m_nBars*phi_slice_ID)/m_nBars);      
      int crystal_id   = fabs(index) - m_nBarrelTiming_Z*m_nBars*phi_slice_ID - m_nBars*module_id;
      
      //front layer gives fine granularity along z
      if (layer_ID == 1)
      {          
          //use center of module for phi
          phi = ((double)phi_slice_ID) * m_phiUnitTimingBarrel - 90.;          
          
          //use bar position for z
          if (sign<0)      z_pos = -sign*((module_id)*m_barLength + sign*(crystal_id+0.5)*m_barWidth);          
          else if (sign>0) z_pos = -sign*((module_id-1)*m_barLength + sign*(crystal_id+0.5)*m_barWidth);          
          x_pos = (timing_radius-m_barWidth/2)*cos(phi*TMath::DegToRad());
          y_pos = (timing_radius-m_barWidth/2)*sin(phi*TMath::DegToRad());
//           std::cout << "index = " << index << " :: phi_slice_ID = " << phi_slice_ID << " :: module_id = " << module_id << " :: crystal_id = " << crystal_id << " :: phi = " << phi << " :: x,y,z = (" << x_pos << ", " << y_pos << ", " << z_pos << ") " << std::endl;
      }
      //rear layer gives fine granularity along phi
      else if (layer_ID == 2)
      {
          //use bar position for phi
          phi = ((double)phi_slice_ID+0.5) * m_phiUnitTimingBarrel - 90. - (crystal_id+0.5)*m_phiUnitTimingBarrel/m_nBars;                   
          //use center of module for z
          z_pos = -sign*(module_id-0.5)*m_barLength;          
          x_pos = (timing_radius-m_barWidth/2)*cos(phi*TMath::DegToRad());
          y_pos = (timing_radius-m_barWidth/2)*sin(phi*TMath::DegToRad()); 
      }                    
  }
  
  else if (fabs(index)>=m_chOffset)  //hit in endcap
  {
//         std::cout << " endcap" <<std::endl;      
      int endcap_index = fabs(index)-m_chOffset;      
      int module_id  = floor(endcap_index/m_nBars);      
      int crystal_id = endcap_index - m_nBars*module_id;      
      int module_X   = floor(module_id/m_nEndcapModulePerLine);
      int module_Y   = module_id - m_nEndcapModulePerLine*module_X;
            
      if (layer_ID == 1)
      {
            x_pos = (module_X+0.5-m_nEndcapModulePerLine/2)*m_barLength;
            y_pos = -sign*((module_Y+1-m_nEndcapModulePerLine/2)*m_barLength-(crystal_id+0.5)*m_barWidth);
            z_pos = -sign*(timing_radius-m_barWidth/2);
      }
      else if (layer_ID == 2)
      {
            x_pos = ((module_X+1-m_nEndcapModulePerLine/2)*m_barLength-(crystal_id+0.5)*m_barWidth);
            y_pos = -sign*(module_Y+0.5-m_nEndcapModulePerLine/2)*m_barLength;
            z_pos = -sign*(timing_radius+m_barWidth/2);
      }
      
//       std::cout << "endcap_index = " << endcap_index << " :: module_id = " << module_id << " ::  module_X = " << module_X  << " :: module_Y = " << module_Y  << " :: posX = " << x_pos << " :: posY = " << y_pos << " :: posZ = " << z_pos <<  std::endl;
      
  }  

//    if (sign<0)                  
       crystal.SetXYZ(x_pos, y_pos, z_pos);  
  return crystal;
}


TVector3 SCEPCal_GeometryHelper::GetCrystalTimingBothVec(long int index_1, long int index_2)
{                    
    
  TVector3 crystal;      
  int sign = index_1/fabs(index_1);  //get index sign (to distinguish left and right sides)
  double phi;
  double z_pos, x_pos, y_pos;
  double timing_radius = (m_BarrelTiming_IR+m_BarrelTiming_OR)/2.;
  
  if (fabs(index_1)<m_chOffset && fabs(index_2)<m_chOffset)       // index in barrel
  {
      int phi_slice_ID_1 = floor(fabs(index_1)/(m_nBarrelTiming_Z*m_nBars));      
      int module_id_1    = floor((fabs(index_1) - m_nBarrelTiming_Z*m_nBars*phi_slice_ID_1)/m_nBars);      
      int crystal_id_1   = fabs(index_1) - m_nBarrelTiming_Z*m_nBars*phi_slice_ID_1 - m_nBars*module_id_1;
      
      int phi_slice_ID_2 = floor(fabs(index_2)/(m_nBarrelTiming_Z*m_nBars));      
      int module_id_2    = floor((fabs(index_2) - m_nBarrelTiming_Z*m_nBars*phi_slice_ID_2)/m_nBars);      
      int crystal_id_2   = fabs(index_2) - m_nBarrelTiming_Z*m_nBars*phi_slice_ID_2 - m_nBars*module_id_2;
                                  
      //front layer gives fine granularity along z            
      if (sign<0)      z_pos = -sign*((module_id_1)*m_barLength + sign*(crystal_id_1+0.5)*m_barWidth);          
      else if (sign>0) z_pos = -sign*((module_id_1-1)*m_barLength + sign*(crystal_id_1+0.5)*m_barWidth);          
      //rear layer gives fine granularity along phi                
      phi = ((double)phi_slice_ID_2+0.5) * m_phiUnitTimingBarrel - 90. - (crystal_id_2+0.5)*m_phiUnitTimingBarrel/m_nBars;
            
      x_pos = (timing_radius-m_barWidth/2)*cos(phi*TMath::DegToRad());
      y_pos = (timing_radius-m_barWidth/2)*sin(phi*TMath::DegToRad());      
  }
  
  if (fabs(index_1) >=m_chOffset && fabs(index_2)>=m_chOffset)       // index in endcap
  {
//         std::cout << " endcap" <<std::endl;      
      int endcap_index_1 = fabs(index_1)-m_chOffset;      
      int endcap_index_2 = fabs(index_2)-m_chOffset;      
      
      int module_id_1  = floor(endcap_index_1/m_nBars);      
      int crystal_id_1 = endcap_index_1 - m_nBars*module_id_1;
      int module_X_1   = floor(module_id_1/m_nEndcapModulePerLine);
      int module_Y_1   = module_id_1 - m_nEndcapModulePerLine*module_X_1;
      
      int module_id_2  = floor(endcap_index_2/m_nBars);      
      int crystal_id_2 = endcap_index_2 - m_nBars*module_id_2;
      int module_X_2   = floor(module_id_2/m_nEndcapModulePerLine);
      
      x_pos = ((module_X_2+1-m_nEndcapModulePerLine/2)*m_barLength-(crystal_id_2+0.5)*m_barWidth);      
      y_pos = -sign*((module_Y_1+1-m_nEndcapModulePerLine/2)*m_barLength-(crystal_id_1+0.5)*m_barWidth);            
      z_pos = -sign*(timing_radius);
//       std:: cout << "endcap_index = " << endcap_index << " :: module_id = " << module_id << " ::  module_X = " << module_X  << " :: module_Y = " << module_Y  << " :: posX = " << x_pos << " :: posY = " << y_pos << " :: posZ = " << z_pos <<  std::endl;
      
  }

  crystal.SetXYZ(x_pos, y_pos, z_pos);  
  return crystal;
}



int SCEPCal_GeometryHelper::GetTimingPhiSliceID(long int index)
{                    
  int phi_slice_ID = floor(fabs(index)/(m_nBarrelTiming_Z*m_nBars));
  return phi_slice_ID;
}

int SCEPCal_GeometryHelper::GetTimingModuleID(long int index)
{                    
  int phi_slice_ID = floor(fabs(index)/(m_nBarrelTiming_Z*m_nBars));
  int module_id    = floor((fabs(index) - m_nBarrelTiming_Z*m_nBars*phi_slice_ID)/m_nBars);  
  return module_id;
}

int SCEPCal_GeometryHelper::GetTimingCrystalID(long int index)
{                    
  int phi_slice_ID = floor(fabs(index)/(m_nBarrelTiming_Z*m_nBars));
  int module_id    = floor((fabs(index) - m_nBarrelTiming_Z*m_nBars*phi_slice_ID)/m_nBars);
  int crystal_id   = fabs(index) - m_nBarrelTiming_Z*m_nBars*phi_slice_ID - m_nBars*module_id;  
  
//   std::cout << " crystal_id = " << crystal_id << std::endl;
  return crystal_id;
}


/*
void SCEPCal_GeometryHelper::PrintGeometry()
{
  std::cout << "\n\n\nGeometryHelper\n\n\n" << std::endl;
  std::cout << "Nb Of Barrel crystals = " << m_NbOfBarrel << std::endl;      
  std::cout << "Nb Of EndCap crystals = " << m_NbOfEndCap << std::endl;
  std::cout << "Nb Of Phi sliced = " << m_NZrot << std::endl;                  
  std::cout << "Total number of crystals = " << m_totCrystal << std::endl;
  std::cout << "Crystal separation in theta = " << m_deltaTheta << std::endl;
  std::cout << "Slice separation in phi = " << m_phiUnit << std::endl;
}*/
