#include "recoUtils.hh"
#include "TMath.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TGraph.h"

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
void  CalSeed::SetEneFront (float this_eneF)  {  eneF = this_eneF;}
void  CalSeed::SetEneRear (float this_eneR)  {  eneR = this_eneR;}
void  CalSeed::SetEne3x3(float this_ene3x3)      {  ene3x3= this_ene3x3;}
void  CalSeed::SetEne3x3Front(float this_ene3x3F)      {  ene3x3F= this_ene3x3F;}
void  CalSeed::SetEne3x3Rear(float this_ene3x3R)      {  ene3x3R= this_ene3x3R;}
void  CalSeed::SetWeighedPhi(float this_phi)     {  weighed_phi= this_phi;}
void  CalSeed::SetWeighedTheta(float this_theta)  {  weighed_theta= this_theta;}
void  CalSeed::AddGenMatch (int thisPdgId){  gen_matched_pdgId.push_back(thisPdgId);}
void  CalSeed::AddGenEne   (float thisGenEne){  gen_matched_ene.push_back(thisGenEne);}

int   CalSeed::GetHitId ()  {  return hit_id;}
int   CalSeed::GetSide ()   {  return side;}
float CalSeed::GetTheta()   {  return theta;}
float CalSeed::GetEta()     {  return -log(tan(theta/2));}
float CalSeed::GetPhi()     {  return phi;}
float CalSeed::GetEne()     {  return ene;}
float CalSeed::GetEneFront()     {  return eneF;}
float CalSeed::GetEneRear()     {  return eneR;}


float CalSeed::GetEne3x3()      {return ene3x3;}
float CalSeed::GetEne3x3Front()      {return ene3x3F;}
float CalSeed::GetEne3x3Rear()      {return ene3x3R;}
float CalSeed::GetWeighedTheta(){return weighed_theta;}
float CalSeed::GetWeighedPhi()  {return weighed_phi;}

std::vector<int>   CalSeed::GetGenMatch(){return gen_matched_pdgId;}
std::vector<float> CalSeed::GetGenEne()  {return gen_matched_ene;}


void CalCluster::Init (CalSeed this_seed, float deltaR_ECAL, float deltaR_HCAL, int imageSize, std::string type) 
{
  
  seed = this_seed;
  maxDeltaR_ECAL = deltaR_ECAL;
  maxDeltaR_HCAL = deltaR_HCAL;
  myImageSize = imageSize;
  myType = type;
    
  EcalClusterEne = 0;  
  HcalClusterEne = 0;  
  HcalClusterEneNarrow = 0;  
  EcalClusterNHits = 0;
  HcalClusterNHits = 0;    
  totEne = 0;
  
//   caloType;
}
CalCluster::~CalCluster(){};

void CalCluster::Clusterize (std::vector<CalHit> ecalHits, std::vector<CalHit> hcalHits, std::vector<CalHit> ecalHitsF, std::vector<CalHit> ecalHitsR)
{
      
    EcalClusterEne = 0;
    EcalClusterEneF = 0;  
    EcalClusterEneR = 0;  
    HcalClusterEne = 0;  
    HcalClusterEneNarrow = 0;  
    EcalClusterNHits = 0;
    HcalClusterNHits = 0;    
    totEne = 0;
    
    TH2F * hImage_E1_temp  = new TH2F ("hImage_E1_temp", "hImage_E1_temp", myImageSize, -bin_width_theta_EC*myImageSize/2, bin_width_theta_EC*myImageSize/2, myImageSize, -bin_width_phi_EC*myImageSize/2, bin_width_phi_EC*myImageSize/2);
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
        
        if (dd < maxDeltaR_ECAL)
        {
            EcalClusterEne+=this_hit.GetEne();
            EcalClusterEneF+=this_hitF.GetEne();
            EcalClusterEneR+=this_hitR.GetEne();
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
        if (dd < maxDeltaR_ECAL)
        {
            HcalClusterEneNarrow+=this_hit.GetEne();
        }
        if (dd < maxDeltaR_HCAL)
        {
            HcalClusterEne+=this_hit.GetEne();
            HcalClusterNHits++;            
        }
    }
    totEne = EcalClusterEne+HcalClusterEne;
    
    
    /*
    if   (     EcalClusterEne  > 0.25
            && seed.GetEne()/EcalClusterEne < 0.8
//             && EcalClusterEne/(EcalClusterEne+HcalClusterEneNarrow) >0.5
         )
    {
        myType = "ecal";
    }
    else
    {
        myType = "hcal";
    }*/
    
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

CalSeed CalCluster::GetSeed()  {return seed;}
float CalCluster::GetTotEne(){return totEne;}
float CalCluster::GetTotEneNarrow(){return (HcalClusterEneNarrow+EcalClusterEne);}
float CalCluster::GetEcalClusterEne() {return EcalClusterEne;}
float CalCluster::GetEcalClusterEneFront() {return EcalClusterEneF;}
float CalCluster::GetEcalClusterEneRear() {return EcalClusterEneR;}
float CalCluster::GetHcalClusterEne() {return HcalClusterEne;}
float CalCluster::GetHcalClusterEneNarrow() {return HcalClusterEneNarrow;}
float CalCluster::GetEcalClusterNHits() {return EcalClusterNHits;}
float CalCluster::GetHcalClusterNHits() {return HcalClusterNHits;}
std::string CalCluster::GetType() {return myType;}
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





std::vector<CalSeed> MakeSuperSeeds (std::vector<CalSeed> allSeeds, std::vector<CalHit> allHits, std::vector<CalHit> allHitsF, std::vector<CalHit> allHitsR, float deltaR, float superSeedTh)
{

    std::vector<CalSeed> SuperSeeds;
    std::vector<CalHit> thisHits = allHits;
    
    
    float R3x3 = 0.006;
    
    for (long unsigned int iseed = 0; iseed < allSeeds.size(); iseed++)
    {
        CalSeed i_seed = allSeeds.at(iseed);
        float i_theta = i_seed.GetTheta();
        float i_phi   = i_seed.GetPhi();         
        float ene_super_seed  = 0;
        float ene_super_seedF = 0;
        float ene_super_seedR = 0;
        float phi_weighed = 0;
        float theta_weighed = 0;
        float w_tot = 0;
        
//         std::cout << " seed[ " << iseed << "] with Ene = " << i_seed.GetEne() << std::endl;
          
//         for (auto i_hit : thisHits)
        for (unsigned int i = 0; i< thisHits.size(); i++)
        {
            
            CalHit i_hit  = thisHits.at(i);
            CalHit i_hitF = allHitsF.at(i);
            CalHit i_hitR = allHitsR.at(i);
            
            float j_theta = i_hit.GetTheta();
            float j_phi   = i_hit.GetPhi();
            
            float dd = sqrt(pow(j_theta-i_theta,2) + pow(j_phi-i_phi,2));
            
            if (dd < R3x3)
            {
                float this_ene = i_hit.GetEne();
                
                ene_super_seed += this_ene;
                phi_weighed    += j_phi*this_ene;
                theta_weighed  += j_theta*this_ene;
                w_tot          += this_ene;                
                
                ene_super_seedF += i_hitF.GetEne();
                ene_super_seedR += i_hitR.GetEne();
            }
        }
        if (w_tot>0.)
        {
            phi_weighed/=w_tot;
            theta_weighed/=w_tot;
            i_seed.SetWeighedPhi(phi_weighed);
            i_seed.SetWeighedTheta(theta_weighed);
            
            
            i_seed.SetEne3x3(ene_super_seed);
            i_seed.SetEne3x3Front(ene_super_seedF);
            i_seed.SetEne3x3Rear(ene_super_seedR);
            
//             std::cout << "iseed = " << iseed << " \n" << std::endl;
//             std::cout << "      --> ene = " << i_seed.GetEne()     << " :: super_ene = "   << i_seed.GetEne3x3() << std::endl;
//             std::cout << "      --> phi = " << i_seed.GetPhi()     << " :: weigh_phi = "   << i_seed.GetWeighedPhi() << std::endl;
//             std::cout << "      --> theta = " << i_seed.GetTheta() << " :: weigh_theta = " << i_seed.GetWeighedTheta() << std::endl;
//             std::cout << "***************************************************************************** " << std::endl;
//         
            if (ene_super_seed>superSeedTh)  SuperSeeds.push_back(i_seed);
        }
    }    
        
    return SuperSeeds;
}


TGraph * getEquivalentTrajectory (float B, float px, float py, float pz, float charge, float maxR)
{


    float pSum = sqrt(px*px+py*py+pz*pz);
    float pT   = sqrt(px*px+py*py);
    float h = -charge/abs(charge);
    float R;
    if (B>0.) R = pT/fabs(charge)/(0.3*B)*1000;
    else      R = 10000000000;

//     std::cout << "px = " << px << " :: py = " << py << " :: pz = " << pz <<  std::endl;
//     std::cout << "bending radius for pT = " << pT << " : " << R/1000 <<  " m" << std::endl;
    float y0 = 0;
    float x0 = 0;
    float z0 = 0;
    float phi0 = atan(py/px)-h*M_PI/2;
//     float phi0;
    if (px<0. && py <0.)   {phi0 = phi0 - M_PI;}
    if (px<0. && py >0.)   {phi0 = M_PI + phi0;}
//     if =  acos(px/pT)-M_PI/2;// + M_PI;

    float lambda = acos(pT/pSum);


    TGraph* gTraj   = new TGraph();
    TGraph* gEqTraj = new TGraph();

    for (int i = 0; i <100; i++)
    {
        float to_m = 100;
        float x;
        float y;
        float z;

//         std::cout << "i = " << i << std::endl;

        if (B>0.)
        {
            x = x0 + R*(cos(phi0+h*i*to_m*cos(lambda)/R) - cos(phi0) );
            y = y0 + R*(sin(phi0+h*i*to_m*cos(lambda)/R) - sin(phi0) );
            z = z0 + i*to_m*sin(lambda);
//             std::cout <<" x  = " << x << " :: y = " << y << std::endl;
        }
        else if (B == 0.)
        {
            x = x0 + i*px/pSum*to_m;
            y = y0 + i*py/pSum*to_m;
//             std::cout <<" x  = " << x << " :: y = " << y << std::endl;
        }

        if (sqrt(x*x+y*y)<maxR)
        {
            gTraj->SetPoint(gTraj->GetN(), x, y);
        }
        else break;
    }

    Double_t impact_x, impact_y;
    gTraj->GetPoint(gTraj->GetN()-1, impact_x, impact_y);

    gEqTraj->SetPoint(0, 0, 0);
    gEqTraj->SetPoint(1, impact_x, impact_y);



    return gEqTraj;
}

