#include "ppfa.hh"
#include "TMath.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TCanvas.h"

#include <unistd.h>

#include "fastjet/ClusterSequence.hh"
#include "TRandom.h"
#include "TF1.h"


using namespace std;
using namespace fastjet; 

// unsigned int microseconds;

std::vector<std::pair<PseudoJet, PseudoJet>> sorted_by_dd (std::vector<std::pair<PseudoJet, PseudoJet>> myHits, PseudoJet refJet)
{
    std::vector<std::pair<PseudoJet, PseudoJet>> sortedHits;

//     std::cout << " start sorting hits by distance, " << myHits.size() << " hits in input collection" << std::endl;
    
    if (myHits.size() ==  1) sortedHits.push_back(myHits.at(0));
    
    while(myHits.size() > 1)
    {
        float minDD = 99999999;
        int hit_id = 0;
        int minID = -999;
        for (auto hit : myHits)
        {
            PseudoJet scint_hit = hit.first;
//             PseudoJet cher_hit = hit.second;
            
            float thisDD = scint_hit.delta_R(refJet);
//             std::cout << " size of collection: " << myHits.size() << " :: hit_id = " << hit_id << " ::  thisDD = " << thisDD << " :: minID = " << minID <<  std::endl;
            if (thisDD < minDD)
            {
                minDD = thisDD;
                minID = hit_id;
            }
        
            hit_id++;
        }
        
        if (minID >-999)
        {
            sortedHits.push_back(myHits.at(minID));        
            myHits.erase(myHits.begin() + minID);
        }
    }
    
    return sortedHits;
}




std::vector<PseudoJet> RunProtoPFA (std::vector<PseudoJet> chargedTracks, std::vector<std::pair<PseudoJet, PseudoJet>> hitsForJet, TH1F* h1SwappedTrackFrac, TH1F *h1ResidualCharged, TH1F *h1ResidualTotCharged)
{

    float hcal_stoch = 0.30;
    float hcal_const = 0.023;
    TF1 * funcHcalRes = new TF1 ("funcHcalRes", "sqrt(pow([0]/sqrt(x),2) + pow([1],2))", 0, 300);
    funcHcalRes->SetParameters(hcal_stoch, hcal_const);
    
    TF1 * funcTotHadRawResponse = new TF1 ("funcTotHadRawResponse", "[3] - [0]*pow(x,[1]) - [2]*x", 0, 200);
    funcTotHadRawResponse->SetParameters(0.46923, -0.962149, -0.000170362 ,0.97);

    
    float eneResponse = 0.95;
    float maxDeltaR_ECAL = 0.05;
    float maxDeltaR_HCAL = 1;
    
      
    int flag_MCT = 0;
    int flag_JHS = 1;
    int flag_JHC = 2;
    int flag_JES = 3;
    int flag_JEC = 4;
    
    float trueTotCharged = 0;
    float recoTotCharged = 0;
    
    double x_factor_hcal = 0.445;
    double x_factor_ecal = 0.371;
    
    
    
    std::vector<std::pair<PseudoJet, PseudoJet>> leftCaloHits = hitsForJet;
    std::vector<PseudoJet> pfaCollection;
    
//     std::vector<PseudoJet> sortedTracks = sorted_by_pt(chargedTracks);
    std::vector<PseudoJet> sortedTracks = chargedTracks;
    
//     std::cout << "pfa algorithm initialized" << std::endl;
    
    int i_track = 0;
    int swappedTrack = 0;
    for (auto track : sortedTracks)
    {
    
        float totCaloE = 0.;
        float trueEne = track.E();
        float targetEne = trueEne*eneResponse;
//         float targetEne = trueEne*funcTotHadRawResponse->Eval(trueEne);
        
//         std::cout << "expected response for " << trueEne << " --> " << funcTotHadRawResponse->Eval(trueEne) << std::endl;
//         float smearedEne = gRandom->Gaus(targetEne, funcHcalRes->Eval(targetEne)*trueEne);
//         if (smearedEne>0) targetEne = smearedEne;
        
        std::vector<std::pair<PseudoJet, PseudoJet>> sortedHits = sorted_by_dd(leftCaloHits, track);
        
        
        std::vector<PseudoJet> thisTrackHits;
        std::vector<int> hitsForRemoval;
        
        
        leftCaloHits.clear();
        std::vector<std::pair<PseudoJet,PseudoJet>> matchedCaloHits;
//         std::cout << i_track << " :: before:: size sorted = " << sortedHits.size() << " :: size left =  " << leftCaloHits.size() << std::endl;
        
        for (auto hit : sortedHits)
        {
//             float deltaDD = sqrt(pow(hit.theta()-track.theta(),2) + pow(hit.phi()-track.phi(),2));
            
            PseudoJet scint_hit = hit.first;
            PseudoJet cher_hit  = hit.second;
            float deltaR  = scint_hit.delta_R(track);
//             std::cout << "deltaR = " << hit.delta_R(track) << " :: deltaDD = "  << deltaDD << std:: endl;
            
            
            if (totCaloE < targetEne &&
                ((scint_hit.user_index()==flag_JES  && deltaR < maxDeltaR_ECAL) ||
                 (scint_hit.user_index()==flag_JHS  && deltaR < maxDeltaR_HCAL) )
            )
            {
                double E_DRO;
                if      (scint_hit.user_index()==flag_JES) E_DRO = (scint_hit.E()-x_factor_ecal*cher_hit.E() )/(1-x_factor_ecal);
                else if (scint_hit.user_index()==flag_JHS) E_DRO = (scint_hit.E()-x_factor_hcal*cher_hit.E() )/(1-x_factor_hcal);
                
//                 totCaloE += scint_hit.E();                
                totCaloE += E_DRO;
                matchedCaloHits.push_back(hit);
            }
            else
            {
                leftCaloHits.push_back(hit);
            }
        }
        
        
//         std::cout << "size sorted = " << sortedHits.size() << " :: size left =  " << leftCaloHits.size() << std::endl;
        if (totCaloE>0.005)
        {
            //matching was good enough
//             if (fabs(totCaloE-trueEne)/trueEne <funcHcalRes->Eval(trueEne)) 
            if (fabs(totCaloE-targetEne)/targetEne <funcHcalRes->Eval(targetEne)*0.75) 
            {
                pfaCollection.push_back(track);
                
//                 h2ScatterPFACharged->Fill(trueEne, totCaloE/trueEne );
                h1ResidualCharged->Fill((totCaloE-trueEne)/trueEne);
                
                trueTotCharged += trueEne;
                recoTotCharged += totCaloE;
                
                swappedTrack++;
            }
            
            //matching rejected
            else
            {
                for (auto non_accepted_hit : matchedCaloHits)
                {
                    leftCaloHits.push_back(non_accepted_hit);
                }
            }
        }
        i_track++;
    }
    
    if (sortedTracks.size()>0) 
    {
//         std::cout << "swapped tracks = " << swappedTrack << " / " << sortedTracks.size() << " = " << float (swappedTrack)/sortedTracks.size() << std::endl;
        h1SwappedTrackFrac->Fill(float (swappedTrack)/sortedTracks.size() );
    }
    
    //add neutral hits - left over calo hits not matched to any charged track
    for (auto neutral_hit : leftCaloHits)
    {
        pfaCollection.push_back(neutral_hit.first);
        pfaCollection.push_back(neutral_hit.second);
    }
    
    
    if (trueTotCharged>0.) h1ResidualTotCharged->Fill((recoTotCharged-trueTotCharged)/trueTotCharged);
    
    return pfaCollection;
}










