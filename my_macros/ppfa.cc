#include "ppfa.hh"
#include "TMath.h"
#include "TH2F.h"
#include "TCanvas.h"

#include <unistd.h>

#include "fastjet/ClusterSequence.hh"
#include "TRandom.h"
#include "TF1.h"


using namespace std;
using namespace fastjet; 

// unsigned int microseconds;

std::vector<PseudoJet> sorted_by_dd (std::vector<PseudoJet> myHits, PseudoJet refJet)
{
    std::vector<PseudoJet> sortedHits;
                
    while(myHits.size()>0)
    {
        float minDD = 99999;
        int hit_id = 0;
        int minID = -999;
        for (auto hit : myHits)
        {
            float thisDD = hit.delta_R(refJet);
            if (thisDD < minDD)
            {
                minDD = thisDD;
                minID = hit_id;
            }
        
            hit_id++;
        }
        sortedHits.push_back(myHits.at(minID));        
        myHits.erase(myHits.begin() + minID);
    }
    
    return sortedHits;
}




std::vector<PseudoJet> RunProtoPFA (std::vector<PseudoJet> chargedTracks, std::vector<PseudoJet> hitsForJet)
{

    float hcal_stoch = 0.30;
    float hcal_const = 0.023;
    TF1 * funcHcalRes = new TF1 ("funcHcalRes", "sqrt(pow([0]/sqrt(x),2) + pow([1],2))", 0, 300);
    funcHcalRes->SetParameters(hcal_stoch, hcal_const);
    
    float eneResponse = 1;
    float maxDeltaR_ECAL = 0.013;
    float maxDeltaR_HCAL = 0.03;
    
      
    int flag_MCT = 0;
    int flag_JHS = 1;
    int flag_JHC = 2;
    int flag_JES = 3;
    int flag_JEC = 4;
    
    
    std::vector<PseudoJet> leftCaloHits = hitsForJet;
    std::vector<PseudoJet> pfaCollection;
    
    std::vector<PseudoJet> sortedTracks = sorted_by_pt(chargedTracks);
    
    for (auto track : sortedTracks)
    {
    
        float totCaloE = 0.;
        float trueEne = track.E();
//         float targetEne = gRandom->Gaus(trueEne, funcHcalRes->Eval(trueEne)*trueEne)*eneResponse;
        float targetEne = trueEne*eneResponse;
        
        std::vector<PseudoJet> sortedHits = sorted_by_dd(leftCaloHits, track);
        
        
        std::vector<PseudoJet> thisTrackHits;
        std::vector<int> hitsForRemoval;
        
        
        leftCaloHits.clear();
        
        for (auto hit : sortedHits)
        {
            if (totCaloE < targetEne &&
                ((hit.user_index()==flag_JES  && hit.delta_R(track) < maxDeltaR_ECAL) ||
                 (hit.user_index()==flag_JHS  && hit.delta_R(track) < maxDeltaR_HCAL) )
            )
            {
                totCaloE += hit.E();                
            }
            else
            {
                leftCaloHits.push_back(hit);
            }
        }
        
        pfaCollection.push_back(track);
    }
    
    //add neutral hits - left over calo hits not matched to any charged track
    for (auto neutral_hit : leftCaloHits)
    {
        pfaCollection.push_back(neutral_hit);
    }
    
    
    return pfaCollection;
}









