#include "recoUtils.hh"
#include "TMath.h"




void CalHit::Init (int this_hit_id, float this_theta, float this_phi, float this_ene) 
{
  hit_id = this_hit_id;
  theta  = this_theta;
  phi    = this_phi;
  ene    = this_ene;  
  side   = 0;
}

void  CalHit::SetHitId (int this_hit_id) {hit_id = this_hit_id;}
void  CalHit::SetSide (int this_side) {side = this_side;}
void  CalHit::SetTheta (float this_theta)  { theta = this_theta;}
void  CalHit::SetPhi (float this_phi)  {  phi = this_phi;}
void  CalHit::SetEne (float this_ene)  {  ene = this_ene;}

int   CalHit::GetHitId ()  {  return hit_id;}
int   CalHit::GetSide ()  {  return side;}
float CalHit::GetTheta() {  return theta;}
float CalHit::GetPhi() {  return phi;}
float CalHit::GetEne(){  return ene;}



void CalSeed::Init (int this_hit_id, float this_theta, float this_phi, float this_ene) 
{
  hit_id = this_hit_id;
  theta  = this_theta;
  phi    = this_phi;
  ene    = this_ene;  
  side   = 0;
}

void  CalSeed::SetHitId (int this_hit_id) {hit_id = this_hit_id;}
void  CalSeed::SetSide (int this_side) {side = this_side;}
void  CalSeed::SetTheta (float this_theta)  { theta = this_theta;}
void  CalSeed::SetPhi (float this_phi)  {  phi = this_phi;}
void  CalSeed::SetEne (float this_ene)  {  ene = this_ene;}
void  CalSeed::AddGenMatch (int thisPdgId){  gen_matched_pdgId.push_back(thisPdgId);}

int   CalSeed::GetHitId ()  {  return hit_id;}
int   CalSeed::GetSide ()  {  return side;}
float CalSeed::GetTheta() {  return theta;}
float CalSeed::GetPhi() {  return phi;}
float CalSeed::GetEne(){  return ene;}
std::vector<int> CalSeed::GetGenMatch(){  return gen_matched_pdgId;}


// void CalCluster::Init (int this_hit_id, float this_theta, float this_phi, float this_ene) 
// {
//   hit_id = this_hit_id;
//   theta  = this_theta;
//   phi    = this_phi;
//   ene    = this_ene;  
//   side   = 0;
// }



std::vector<CalSeed> CleanSeeds (std::vector<CalSeed> allSeeds, float deltaR)
{

    std::vector<CalSeed> CleanedSeeds;
    for (long unsigned int iseed = 0; iseed < allSeeds.size(); iseed++)
    {
        CalSeed i_seed = allSeeds.at(iseed);
        float i_theta = i_seed.GetTheta();
        float i_phi   = i_seed.GetPhi();         
        bool maxIsolatedHit = true;
//         std::cout << " seed[ " << iseed << "] with Ene = " << i_seed.GetEne() << std::endl;
          
        for (long unsigned int jseed = 0; jseed < allSeeds.size(); jseed++)
        {
            if (jseed == iseed) continue;
            CalSeed j_seed = allSeeds.at(jseed);
            float j_theta = j_seed.GetTheta();
            float j_phi   = j_seed.GetPhi();
            float dd = sqrt(pow(j_theta-i_theta,2) + pow(j_phi-i_phi,2));
            
            if (dd < deltaR)
            {               
                if (i_seed.GetEne()<j_seed.GetEne()) 
                {
                    maxIsolatedHit = false;                
//                     std::cout << " dd = " << dd << " :: j_theta-i_theta= " << j_theta-i_theta << " :: j_phi-i_phi " << j_phi-i_phi << std::endl;
//                     std::cout << "Found a neighboring seed with higher energy: iEne = " << i_seed.GetEne() << " < jEne = " << j_seed.GetEne() << std::endl;
                }
            }            
        }          
          
        if (maxIsolatedHit)  CleanedSeeds.push_back(i_seed);
    }
    
        
    return CleanedSeeds;
}
