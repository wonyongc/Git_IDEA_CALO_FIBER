#include "recoUtils.hh"
#include "TMath.h"
#include "TH2F.h"
#include "TCanvas.h"

#include <unistd.h>
// unsigned int microseconds;

void CalHit::Init (int this_hit_id, float this_theta, float this_phi, float this_ene) 
{
  hit_id = this_hit_id;
  theta  = this_theta;
  phi    = this_phi;
  ene    = this_ene;  
  side   = 0;
}

CalHit::~CalHit(){};

void  CalHit::SetHitId (int this_hit_id) {hit_id = this_hit_id;}
void  CalHit::SetSide (int this_side) {side = this_side;}
void  CalHit::SetTheta (float this_theta)  { theta = this_theta;}
void  CalHit::SetPhi (float this_phi)  {  phi = this_phi;}
void  CalHit::SetEne (float this_ene)  {  ene = this_ene;}

int   CalHit::GetHitId ()  {  return hit_id;}
int   CalHit::GetSide ()  {  return side;}
float CalHit::GetTheta() {  return theta;}
float CalHit::GetEta() {  return -log(tan(theta/2));}
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

CalSeed::~CalSeed(){};

void  CalSeed::SetHitId (int this_hit_id) {hit_id = this_hit_id;}
void  CalSeed::SetSide (int this_side) {side = this_side;}
void  CalSeed::SetTheta (float this_theta)  { theta = this_theta;}
void  CalSeed::SetPhi (float this_phi)  {  phi = this_phi;}
void  CalSeed::SetEne (float this_ene)  {  ene = this_ene;}
void  CalSeed::AddGenMatch (int thisPdgId){  gen_matched_pdgId.push_back(thisPdgId);}

int   CalSeed::GetHitId ()  {  return hit_id;}
int   CalSeed::GetSide ()  {  return side;}
float CalSeed::GetTheta() {  return theta;}
float CalSeed::GetEta() {  return -log(tan(theta/2));}
float CalSeed::GetPhi() {  return phi;}
float CalSeed::GetEne(){  return ene;}
std::vector<int> CalSeed::GetGenMatch(){  return gen_matched_pdgId;}


void CalCluster::Init (CalSeed this_seed, float deltaR, int imageSize) 
{
  
  seed = this_seed;
  maxDeltaR = deltaR;
  myImageSize = imageSize;
    
  EcalClusterEne = 0;  
  HcalClusterEne = 0;  
  EcalClusterNHits = 0;
  HcalClusterNHits = 0;    
  totEne = 0;
  
//   caloType;
}
CalCluster::~CalCluster(){};

void CalCluster::Clusterize (std::vector<CalHit> ecalHits, std::vector<CalHit> hcalHits, std::vector<CalHit> ecalHitsF, std::vector<CalHit> ecalHitsR)
{
      
    EcalClusterEne = 0;  
    HcalClusterEne = 0;  
    EcalClusterNHits = 0;
    HcalClusterNHits = 0;    
    totEne = 0;
    
    TH2F * hImage_E1_temp  = new TH2F ("hImage_E1_temp", "hImage_E1_temp", myImageSize, -bin_width_theta_EC*(myImageSize/2), bin_width_theta_EC*(myImageSize/2), myImageSize, -bin_width_phi_EC*(myImageSize/2), bin_width_phi_EC*(myImageSize/2));
    TH2F * hImage_E2_temp  = new TH2F ("hImage_E2_temp", "hImage_E2_temp", myImageSize, -bin_width_theta_EC*myImageSize/2, bin_width_theta_EC*myImageSize/2, myImageSize, -bin_width_phi_EC*myImageSize/2, bin_width_phi_EC*myImageSize/2);
    TH2F * hImage_HC_temp  = new TH2F ("hImage_HC_temp", "hImage_HC_temp", myImageSize, -bin_width_theta_HC*myImageSize/2, bin_width_theta_HC*myImageSize/2, myImageSize, -bin_width_phi_HC*myImageSize/2, bin_width_phi_HC*myImageSize/2);
    
    for (unsigned int i = 0; i < ecalHits.size(); i++)
    {
        CalHit this_hit = ecalHits.at(i);
                
        float hit_theta = this_hit.GetTheta();
        float hit_phi   = this_hit.GetPhi();
        float dd = sqrt(pow(hit_theta-seed.GetTheta(),2) + pow(hit_phi-seed.GetPhi(),2));

        CalHit this_hitF = ecalHitsF.at(i);
        hImage_E1_temp ->Fill(hit_theta-seed.GetTheta(), hit_phi-seed.GetPhi(), this_hitF.GetEne());
        CalHit this_hitR = ecalHitsR.at(i);
        hImage_E2_temp ->Fill(hit_theta-seed.GetTheta(), hit_phi-seed.GetPhi(), this_hitR.GetEne());
        
        if (dd < maxDeltaR)
        {
            EcalClusterEne+=this_hit.GetEne();
            EcalClusterNHits++;            
        }        
    }
    
    for (unsigned int i = 0; i < hcalHits.size(); i++)
    {
        CalHit this_hit = hcalHits.at(i);
                
        float hit_theta = this_hit.GetTheta();
        float hit_phi   = this_hit.GetPhi();
        float dd = sqrt(pow(hit_theta-seed.GetTheta(),2) + pow(hit_phi-seed.GetPhi(),2));
        
        hImage_HC_temp ->Fill(hit_theta-seed.GetTheta(), hit_phi-seed.GetPhi(), this_hit.GetEne());
        if (dd < maxDeltaR)
        {
            HcalClusterEne+=this_hit.GetEne();
            HcalClusterNHits++;            
        }
    }
    totEne = EcalClusterEne+HcalClusterEne;
    
        
    for (int iBinX = 0; iBinX<myImageSize; iBinX++)
    {
        for (int iBinY = 0; iBinY<myImageSize; iBinY++)
        {
            int pixel = iBinX+iBinY*myImageSize;
            image_E1[pixel] = hImage_E1_temp->GetBinContent(iBinX+1, iBinY+1);
            image_E2[pixel] = hImage_E2_temp->GetBinContent(iBinX+1, iBinY+1);
            image_HC[pixel] = hImage_HC_temp->GetBinContent(iBinX+1, iBinY+1);
//             std::cout << "pixel: " << pixel << ", content: " << image_E1[pixel] << std::endl;
        }
    }
    
    
    hImage_E1_temp->Delete();
    hImage_E2_temp->Delete();
    hImage_HC_temp->Delete();
}



void  CalCluster::SetSeed(CalSeed this_seed) {seed = this_seed;}

CalSeed   CalCluster::GetSeed()  {return seed;}
float CalCluster::GetTotEne(){return totEne;}
float CalCluster::GetEcalClusterEne() {return EcalClusterEne;}
float CalCluster::GetHcalClusterEne() {return HcalClusterEne;}
float CalCluster::GetEcalClusterNHits() {return EcalClusterNHits;}
float CalCluster::GetHcalClusterNHits() {return HcalClusterNHits;}
float* CalCluster::GetImage(std::string segment_name) 
{
    if (segment_name == "E1") return image_E1;
    if (segment_name == "E2") return image_E2;
    if (segment_name == "HC") return image_HC;
    else return NULL;        
}




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
