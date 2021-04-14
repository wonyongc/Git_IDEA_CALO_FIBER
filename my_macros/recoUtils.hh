#include <fstream>
#include <iostream>
#include <map>
#include <vector>
#include <map>
#include <string>
#include <cstdlib>
#include <stdlib.h>
#include <utility>
#include <algorithm>
#include "TMath.h"



class CalHit
{
    int hit_id;
    int side;
    float theta;
    float phi;
    float ene;    
    
    public:
        void Init(int, float, float, float);
        ~CalHit();
        void SetHitId(int);
        void SetSide(int);
        void SetTheta(float);
        void SetPhi(float);
        void SetEne(float);
        
        int GetHitId();
        int GetSide();
        float GetTheta();
        float GetEta();
        float GetPhi();
        float GetEne();            
};


class CalSeed
{
    int hit_id;
    int side;
    float theta;
    float phi;
    float ene;
    std::vector<int> gen_matched_pdgId;
    
    public:
        void Init(int, float, float, float);
        ~CalSeed();
        void SetHitId(int);
        void SetSide(int);
        void SetTheta(float);
        void SetPhi(float);
        void SetEne(float);
        
        int GetHitId();
        int GetSide();
        float GetTheta();
        float GetEta();
        float GetPhi();
        float GetEne();
        std::vector<int> GetGenMatch();
        
        void AddGenMatch(int);
    
};


std::vector<CalSeed> CleanSeeds (std::vector<CalSeed> allSeeds, float deltaR);



class CalCluster
{
    
    CalSeed seed;
    float maxDeltaR_ECAL;
    float maxDeltaR_HCAL;
    
    float EcalClusterEne;
    float EcalClusterEneF;
    float EcalClusterEneR;
    float HcalClusterEne;
    float HcalClusterEneNarrow;
    
    float EcalClusterNHits;
    float HcalClusterNHits;
    
    float totEne;
    std::string myType = "none";
    
    
    int NPHI_EC    = 1130;
    int NTHETA_EC  = 180+180;
    int NPHI_DRT   = 36;
    int NTHETA_DRT = 40+40;
    
    double bin_width_theta_EC = M_PI/NTHETA_EC;
    double bin_width_theta_HC = M_PI/NTHETA_DRT;
    
    double bin_width_phi_EC = 2*M_PI/NPHI_EC;
    double bin_width_phi_HC = 2*M_PI/NPHI_DRT;
    
    int myImageSize = 15;
    float image_E1[225];
    float image_E2[225];
    float image_HC[225];
    
    
    public:
        void Init(CalSeed mySeed, float deltaR_ECAL, float deltaR_HCAL, int imageSize, std::string type);
        void Clusterize (std::vector<CalHit> ecalHits, std::vector<CalHit> hcalHits, std::vector<CalHit> ecalHitsF, std::vector<CalHit> ecalHitsR);
        ~CalCluster();
        void SetSeed(CalSeed this_seed);        
                
        CalSeed GetSeed();
        float GetTotEne();
        float GetTotEneNarrow();
        float GetEcalClusterEne();
        float GetEcalClusterEneFront();
        float GetEcalClusterEneRear();
        float GetHcalClusterEne();
        float GetHcalClusterEneNarrow();
        float GetEcalClusterNHits();
        float GetHcalClusterNHits();
        std::string GetType();
        float* GetImage(std::string segment_name);
                    
};



std::vector<CalSeed> CleanSeeds (std::vector<CalSeed> allSeeds, float deltaR);
