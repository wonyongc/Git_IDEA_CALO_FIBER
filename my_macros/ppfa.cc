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


/*
std::vector<PseudoJet> photonFinder (std::vector<PseudoJet> chargedTracks, std::vector<std::pair<PseudoJet, PseudoJet>> hitsForJet, 
                                    float x_factor_ecal, float x_factor_hcal,
                                    TH1F* h1SwappedTrackFrac, TH1F *h1ResidualCharged, TH1F *h1ResidualTotCharged)*/


std::vector<PseudoJet> RunProtoPFA (std::vector<PseudoJet> chargedTracks, std::vector<std::pair<PseudoJet, PseudoJet>> hitsForJet, 
                                    float x_factor_ecal, float x_factor_hcal,
                                    TH1F* h1SwappedTrackFrac, TH1F *h1ResidualCharged, TH1F *h1ResidualTotCharged)
{

    float hcal_stoch = 0.30;
    float hcal_const = 0.023;
    TF1 * funcHcalRes = new TF1 ("funcHcalRes", "sqrt(pow([0]/sqrt(x),2) + pow([1],2))", 0, 300);
    funcHcalRes->SetParameters(hcal_stoch, hcal_const);
    
    TF1 * funcTotHadRawResponse = new TF1 ("funcTotHadRawResponse", "[3] - [0]*pow(x,[1]) - [2]*x", 0, 200);
    funcTotHadRawResponse->SetParameters(0.46923, -0.962149, -0.000170362 ,0.97);

    
    float eneResponse = 0.99;
    float maxDeltaR_ECAL = 0.05;
    float maxDeltaR_HCAL = 0.3;
    
      
    int flag_MCT = 0;
    int flag_JHS = 1;
    int flag_JHC = 2;
    int flag_JES = 3;
    int flag_JEC = 4;
    
    float trueTotCharged = 0;
    float recoTotCharged = 0;
    
    
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
            
            double this_E_DRO;
            if      (scint_hit.user_index()==flag_JES) this_E_DRO = (scint_hit.E()-x_factor_ecal*cher_hit.E() )/(1-x_factor_ecal);
            else if (scint_hit.user_index()==flag_JHS) this_E_DRO = (scint_hit.E()-x_factor_hcal*cher_hit.E() )/(1-x_factor_hcal);
            
            if (totCaloE < targetEne &&
                fabs(totCaloE+this_E_DRO-targetEne) < fabs(totCaloE-targetEne) &&
                ((scint_hit.user_index()==flag_JES  && deltaR < maxDeltaR_ECAL) ||
                 (scint_hit.user_index()==flag_JHS  && deltaR < maxDeltaR_HCAL) )
            )
            {
                totCaloE += this_E_DRO;
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




Double_t rms90(TH1F *h) 
{
    TAxis *axis = h->GetXaxis();
    Int_t nbins = axis->GetNbins();
    Int_t imean = axis->FindBin(h->GetMean());
    Double_t entries = 0.9*h->GetEntries();
    Double_t w = h->GetBinContent(imean);
    Double_t x = h->GetBinCenter(imean);
    Double_t sumw = w;
    Double_t sumwx = w*x;
    Double_t sumwx2 = w*x*x;
    for (Int_t i=1;i<nbins;i++) 
    {
        if (i> 0) 
        {
            w = h->GetBinContent(imean-i);
            x = h->GetBinCenter(imean-i);
            sumw += w;
            sumwx += w*x;
            sumwx2 += w*x*x;
        }
        if (i<= nbins) 
        {
            w = h->GetBinContent(imean+i);
            x = h->GetBinCenter(imean+i);
            sumw += w;
            sumwx += w*x;
            sumwx2 += w*x*x;
        }
        if (sumw > entries) break;
    }
    
    x = sumwx/sumw;
    
    Double_t rms2 = TMath::Abs(sumwx2/sumw -x*x);
    Double_t result = TMath::Sqrt(rms2);
    Double_t rms90 = result;
    
//     printf(“RMS of central 90% = %g, RMS total = %g\n”,result,h->GetRMS());
    
    return rms90;
}




