#include "recoUtils.hh"
#include "TMath.h"




void EcalSeed::Init (int this_crystal_id, float this_theta, float this_phi, float this_ene) 
{
  crystal_id    = this_crystal_id;
  theta         = this_theta;
  phi           = this_phi;
  ene           = this_ene;  
}

void EcalSeed::SetCrystalId (int this_crystal_id) 
{
  crystal_id = this_crystal_id;
}

void EcalSeed::SetTheta (float this_theta) 
{
  theta = this_theta;
}

void EcalSeed::SetPhi (float this_phi) 
{
  phi = this_phi;
}

void EcalSeed::SetEne (float this_ene)
{
  ene = this_ene;
}

int EcalSeed::GetCrystalId () 
{
  return crystal_id;
}

float EcalSeed::GetTheta() 
{
  return theta;
}

float EcalSeed::GetPhi() 
{
  return phi;
}

float EcalSeed::GetEne()
{
  return ene;
}

std::vector<int> EcalSeed::GetGenMatch()
{
  return gen_matched_pdgId;
}



void EcalSeed::AddGenMatch (int thisPdgId)
{
  gen_matched_pdgId.push_back(thisPdgId);
}



// int main () {
//   Rectangle rect;
//   rect.set_values (3,4);
//   cout << "area: " << rect.area();
//   return 0;
// }
