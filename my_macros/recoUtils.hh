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






class EcalSeed
{
    int crystal_id;
    float theta;
    float phi;
    float ene;
    std::vector<int> gen_matched_pdgId;
    
    public:
        void Init(int, float, float, float);
        void SetCrystalId(int);
        void SetTheta(float);
        void SetPhi(float);
        void SetEne(float);
        
        int GetCrystalId();
        float GetTheta();
        float GetPhi();
        float GetEne();
        std::vector<int> GetGenMatch();
        
        void AddGenMatch(int);
    
};
    


// class Rectangle 
// {
//     int width, height;
//   public:
//     void set_values (int,int);
//     int area (void);
// } rect;
