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
        float GetPhi();
        float GetEne();
        std::vector<int> GetGenMatch();
        
        void AddGenMatch(int);
    
};


std::vector<CalSeed> CleanSeeds (std::vector<CalSeed> allSeeds, float deltaR);



class CalCluster
{
    
    CalSeed seed;
    float maxDeltaR;
    
    float EcalClusterEne;
    float HcalClusterEne;
    
    float EcalClusterNHits;
    float HcalClusterNHits;
    
    float totEne;
    
    
    public:
        void Init(CalSeed mySeed, float deltaR);
        void Clusterize (std::vector<CalHit> ecalHits, std::vector<CalHit> hcalHits);
        ~CalCluster();
        void SetSeed(CalSeed this_seed);        
                
        CalSeed GetSeed();
        float GetTotEne();
        float GetEcalClusterEne();
        float GetHcalClusterEne();
        float GetEcalClusterNHits();
        float GetHcalClusterNHits();
                    
};



std::vector<CalSeed> CleanSeeds (std::vector<CalSeed> allSeeds, float deltaR);
