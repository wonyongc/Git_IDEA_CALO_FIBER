#include "ppfa.hh"
#include "TMath.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TCanvas.h"

#include <unistd.h>

#include "fastjet/ClusterSequence.hh"
#include "TRandom.h"
#include "TF1.h"
#include "TGraph.h"
#include "recoUtils.hh"

// #include "SCEPCal_GeometryHelper.hh"


using namespace std;
using namespace fastjet; 

// int flag_MCT = 100;
int flag_JHS = 1;
int flag_JHC = 2;
int flag_JES = 3;
int flag_JEC = 4;

// if particle reached the calorimeter
float ecal_impact_radius = 1900;
float hcal_impact_radius = 1900;

//radial steps for track-hit matching algo
float R_ECAL_cut[5] = {0.005, 0.01, 0.02, 0.035, 0.07};
float R_HCAL_cut[5] = {0.01,  0.05, 0.10, 0.20,  0.5};
    
    

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
            float thisDD = scint_hit.delta_R(refJet);            
            //sort hits
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

// std::pair<std::vector<PseudoJet>,std::vector<std::pair<PseudoJet, PseudoJet>> > RunProtoPFA (std::vector<PseudoJet> chargedTracks, std::vector<std::pair<PseudoJet, PseudoJet>> hitsForJet,
std::pair<std::vector<PseudoJet>,std::vector<PseudoJet> > RunProtoPFA (std::vector<PseudoJet> chargedTracks, std::vector<std::pair<PseudoJet, PseudoJet>> hitsForJet,
                                    float x_factor_ecal, float x_factor_hcal, float Bfield, float matchPFAcut, bool DRO_ON,
                                    TH1F* h1SwappedTrackFrac, TH1F *h1ResidualCharged, TH1F *h1ResidualTotCharged)
{

    float hcal_stoch = 0.30;
    float hcal_const = 0.023;
    TF1 * funcHcalRes = new TF1 ("funcHcalRes", "sqrt(pow([0]/sqrt(x),2) + pow([1],2))", 0, 300);
    funcHcalRes->SetParameters(hcal_stoch, hcal_const);
    
    TF1 * funcTotHadRawResponse = new TF1 ("funcTotHadRawResponse", "[3] - [0]*pow(x,[1]) - [2]*x", 0, 200);
    //old
//     funcTotHadRawResponse->SetParameters(0.46923, -0.962149, -0.000170362 ,0.97);
    
    //raw new
    if (!DRO_ON) funcTotHadRawResponse->SetParameters(-0.35625, 0.193565, 0.00347531 ,0.367476);
    //dro new
    if (DRO_ON)  funcTotHadRawResponse->SetParameters(-0.364088, 0.0465792, 0.000651904, 0.59792);


    
    float eneResponse = 0.99;
    float maxDeltaR_ECAL = 0.05;
    float maxDeltaR_HCAL = 0.3;
//     float Bfield = 2.0;
    
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
    std::vector<PseudoJet> leftOverTracks;
    
//     std::cout << "pfa algorithm initialized" << std::endl;
    
    int i_track = 0;
    int swappedTrack = 0;
    for (auto track : sortedTracks)
    {
    
        float totCaloE = 0.;
        float totCaloS_EC = 0.;
        float totCaloC_EC = 0.;
        float totCaloS_HC = 0.;
        float totCaloC_HC = 0.;
//         double phi   = track.mcs_phi;
//         double eta   = track.mcs_eta->at(i);
//         double px, py;
//         px = pT*cos(phi);
//         py = pT*sin(phi);        
//         double pz = -pT*sinh(eta);
        
        float trueEne = track.E();
//         float targetEne = trueEne*eneResponse;
        float targetEne = trueEne*funcTotHadRawResponse->Eval(trueEne);
        
//         std::cout << "expected response for " << trueEne << " --> " << funcTotHadRawResponse->Eval(trueEne) << std::endl;
//         float smearedEne = gRandom->Gaus(targetEne, funcHcalRes->Eval(targetEne)*trueEne);
//         if (smearedEne>0) targetEne = smearedEne;
        

        //calculate impact point on calorimeter
        
        float charge = track.user_index()/100;
        float pT = track.perp();
    
        if (Bfield>0 && pT/fabs(charge)/(0.3*Bfield)*1000*2<1800)
        {
//             std::cout << "track with pT = " << pT << "  did not reach the calorimeter" << std::endl;
            pfaCollection.push_back(track);
            continue;
        }
        
        // if particle reached the calorimeter
        float ecal_impact_radius = 1900;
        float hcal_impact_radius = 1900;
        
        PseudoJet effectiveTrackEcal;
        PseudoJet effectiveTrackHcal;
        
        // sort calo hits by distance from impact point of the track
        std::vector<std::pair<PseudoJet, PseudoJet>> sortedHits;
        
        if (Bfield>0)
        {
            TGraph* thisTraj = getEquivalentTrajectory (Bfield, track.px(), track.py(), track.pz(), charge, ecal_impact_radius);
            Double_t impact_x, impact_y;
            thisTraj->GetPoint(thisTraj->GetN()-1, impact_x, impact_y);
        
            float impact_phi = atan(impact_y/impact_x);
            if (impact_x<0. && impact_y <0.)   {impact_phi = impact_phi - M_PI;}
            if (impact_x<0. && impact_y >0.)   {impact_phi = M_PI + impact_phi;}        
            double impact_theta = M_PI- 2*atan(exp(-track.eta()));                
            float scale_p = 1./sqrt(impact_x*impact_x + impact_y*impact_y) * sqrt(pow(track.px(),2)+pow(track.py(),2));
            effectiveTrackEcal = PseudoJet(impact_x*scale_p, impact_y*scale_p, track.pz(), trueEne);
        
        
            thisTraj = getEquivalentTrajectory (Bfield, track.px(), track.py(), track.pz(), charge, hcal_impact_radius);        
            thisTraj->GetPoint(thisTraj->GetN()-1, impact_x, impact_y);
            impact_phi = atan(impact_y/impact_x);
            if (impact_x<0. && impact_y <0.)   {impact_phi = impact_phi - M_PI;}
            if (impact_x<0. && impact_y >0.)   {impact_phi = M_PI + impact_phi;}        
            impact_theta = M_PI- 2*atan(exp(-track.eta()));        
            scale_p = 1./sqrt(impact_x*impact_x + impact_y*impact_y) * sqrt(pow(track.px(),2)+pow(track.py(),2));
            effectiveTrackHcal = PseudoJet(impact_x*scale_p, impact_y*scale_p, track.pz(), trueEne);
            
             sortedHits = sorted_by_dd(leftCaloHits, effectiveTrackHcal);
        }
        else 
        {
            sortedHits = sorted_by_dd(leftCaloHits, track);
        }
        
//         std::cout << "n sorted hits found: " << sortedHits.size() << std::endl;
        
        //clear collection of leftover calo hits
        leftCaloHits.clear();
        std::vector<std::pair<PseudoJet,PseudoJet>> matchedCaloHits;
        
        for (auto hit : sortedHits)
        {
            PseudoJet scint_hit = hit.first;
            PseudoJet cher_hit  = hit.second;
//             float deltaR  = scint_hit.delta_R(track);
            float deltaR;
            if      (Bfield > 0 && scint_hit.user_index()==flag_JES) deltaR = scint_hit.delta_R(effectiveTrackEcal);
            else if (Bfield > 0 && scint_hit.user_index()==flag_JHS) deltaR = scint_hit.delta_R(effectiveTrackHcal);
            else if (Bfield == 0)                                    deltaR = scint_hit.delta_R(track);
            
//             float hit_theta = 2*atan(exp(-scint_hit.eta()));
//             float deltaR  = sqrt(pow(scint_hit.phi()-impact_phi,2) + pow(hit_theta-impact_theta,2));;
//             std::cout << "deltaR = " << hit.delta_R(track) << " :: deltaDD = "  << deltaDD << std:: endl;
            
            double this_E_DRO;
            double this_S = scint_hit.E();
            double this_C = cher_hit.E();
            if      (scint_hit.user_index()==flag_JES) this_E_DRO = (scint_hit.E()-x_factor_ecal*cher_hit.E() )/(1-x_factor_ecal);
            else if (scint_hit.user_index()==flag_JHS) this_E_DRO = (scint_hit.E()-x_factor_hcal*cher_hit.E() )/(1-x_factor_hcal);
            
//             this_E_DRO = scint_hit.E();
            
            totCaloE = (totCaloS_EC-x_factor_ecal*totCaloC_EC )/(1-x_factor_ecal) + (totCaloS_HC-x_factor_hcal*totCaloC_HC )/(1-x_factor_hcal);
//             totCaloE = totCaloS_EC+totCaloS_HC;
            if (totCaloE < targetEne &&
            fabs(totCaloE+this_E_DRO-targetEne) < fabs(totCaloE-targetEne) &&
//             if (current_tot_dro < targetEne &&                
//                 fabs(current_tot_dro+this_E_DRO-targetEne) < fabs(current_tot_dro-targetEne) &&
                ((scint_hit.user_index()==flag_JES  && deltaR < maxDeltaR_ECAL) ||
                 (scint_hit.user_index()==flag_JHS  && deltaR < maxDeltaR_HCAL) )
            )
            {
//                 totCaloE += this_E_DRO;
                if (scint_hit.user_index()==flag_JES)
                {
                    totCaloS_EC += this_S;
                    totCaloC_EC += this_C;
                }
                else if (scint_hit.user_index()==flag_JHS)
                {
                    totCaloS_HC += this_S;
                    totCaloC_HC += this_C;
                }
                matchedCaloHits.push_back(hit);
            }
            else
            {
                leftCaloHits.push_back(hit);
            }
        }
        
        totCaloE = (totCaloS_EC-x_factor_ecal*totCaloC_EC )/(1-x_factor_ecal) + (totCaloS_HC-x_factor_hcal*totCaloC_HC )/(1-x_factor_hcal);
//         totCaloE = totCaloS_EC+totCaloS_HC;
//         std::cout << "size sorted = " << sortedHits.size() << " :: size left =  " << leftCaloHits.size() << std::endl;
        if (totCaloE>0.005)
        {
            //matching was good enough
//             if (fabs(totCaloE-trueEne)/trueEne <funcHcalRes->Eval(trueEne)) 
            if (fabs(totCaloE-targetEne)/targetEne <funcHcalRes->Eval(targetEne)*matchPFAcut) 
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

                leftOverTracks.push_back(track);
            }
        }
//         if (sortedHits.size()==0) pfaCollection.push_back(track);
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
    
//     return pfaCollection;
//     return std::make_pair(pfaCollection, leftCaloHits);
    return std::make_pair(pfaCollection, leftOverTracks);
    
}





std::vector<PseudoJet> RunProtoPFA_Iterative (std::vector<PseudoJet> chargedTracks, std::vector<std::pair<PseudoJet, PseudoJet>> hitsForJet, 
                                    float x_factor_ecal, float x_factor_hcal, float Bfield, float matchPFAcut, bool DRO_ON,
                                    TH1F* h1SwappedTrackFrac, TH1F *h1ResidualCharged, TH1F *h1ResidualTotCharged)
{

    std::cout << "ready to run DR-pPFA..." << std::endl;
    
    float hcal_stoch = 0.30;
    float hcal_const = 0.023;
    TF1 * funcHcalRes = new TF1 ("funcHcalRes", "sqrt(pow([0]/sqrt(x),2) + pow([1],2))", 0, 300);
    funcHcalRes->SetParameters(hcal_stoch, hcal_const);
    
    TF1 * funcTotHadRawResponse = new TF1 ("funcTotHadRawResponse", "[3] - [0]*pow(x,[1]) - [2]*x", 0, 200);
    
    //raw new
    if (!DRO_ON) funcTotHadRawResponse->SetParameters(-0.35625, 0.193565, 0.00347531 ,0.367476);
    //dro new
    if (DRO_ON)  funcTotHadRawResponse->SetParameters(-0.364088, 0.0465792, 0.000651904, 0.59792);
    
//     float eneResponse = 0.99;
//     float maxDeltaR_ECAL = 0.05;
//     float maxDeltaR_HCAL = 0.3;
    
    
    
    float trueTotCharged = 0;
    float recoTotCharged = 0;
    
    
    std::vector<std::pair<PseudoJet, PseudoJet>> leftCaloHits = hitsForJet;
    std::vector<PseudoJet> pfaCollection;
    
//     std::vector<PseudoJet> allSortedTracks = sorted_by_E(chargedTracks);
    std::vector<PseudoJet> allSortedTracks = sorted_by_pt(chargedTracks);
    
    
//     std::vector<PseudoJet> inverse_pt;
//     for (unsigned int it = allSortedTracks.size()-1; it>0; it--)
//     {
//          inverse_pt.push_back(allSortedTracks.at(it));
//     }
//     allSortedTracks = inverse_pt;
//     std::vector<PseudoJet> allSortedTracks = chargedTracks;
    
//     std::cout << "pfa algorithm initialized" << std::endl;
    
    // take out tracks not reaching the calorimeter and put them directly in the pfa collection
    std::vector<PseudoJet> sortedTracks;    //tracks reaching the calorimeter
    std::vector<float> totCaloE;
    std::vector<float> totCaloS_EC;
    std::vector<float> totCaloC_EC;
    std::vector<float> totCaloS_HC;
    std::vector<float> totCaloC_HC;
    std::map<int, std::vector<std::pair<PseudoJet,PseudoJet>>> matchedCaloHits;
    
    
    for (auto track : allSortedTracks)
    {
        float charge = track.user_index()/100;
        float pT = track.perp();
    
        if (Bfield>0 && pT/fabs(charge)/(0.3*Bfield)*1000*2<1800)
        {
//             std::cout << "track with pT = " << pT << "  did not reach the calorimeter" << std::endl;
            pfaCollection.push_back(track);
        }
        else
        {
            sortedTracks.push_back(track);
            totCaloE.push_back(0);
            totCaloS_EC.push_back(0);
            totCaloC_EC.push_back(0);
            totCaloS_HC.push_back(0);
            totCaloC_HC.push_back(0);
        }
    }
    
    std::cout << "track and calo hit lists initialized..." << std::endl;
    

    int swappedTrack = 0;
    std::vector<int> swappedTrackList;
    
    for (int iR = 0; iR < 5; iR++)
    {
        
//         float maxDeltaR_ECAL = R_ECAL_cut[iR];
//         float maxDeltaR_HCAL = R_HCAL_cut[iR];
        std::cout << "**************************************************************************" << std::endl;
        std::cout << iR << " :: matching hits with R_ECAL = " << R_ECAL_cut[iR] << " and R_HCAL = " << R_HCAL_cut[iR] << std::endl;
        
        int i_track = 0;
        
        for (auto track : sortedTracks)
        {
            
           if(std::find(swappedTrackList.begin(), swappedTrackList.end(), i_track) != swappedTrackList.end()) 
           {            
                continue;
           }
//         double phi   = track.mcs_phi;
//         double eta   = track.mcs_eta->at(i);
//         double px, py;
//         px = pT*cos(phi);
//         py = pT*sin(phi);        
//         double pz = -pT*sinh(eta);
        
            float trueEne = track.E();
//         float targetEne = trueEne*eneResponse;
            float targetEne = trueEne*funcTotHadRawResponse->Eval(trueEne);
        
//         std::cout << "expected response for " << trueEne << " --> " << funcTotHadRawResponse->Eval(trueEne) << std::endl;
//         float smearedEne = gRandom->Gaus(targetEne, funcHcalRes->Eval(targetEne)*trueEne);
//         if (smearedEne>0) targetEne = smearedEne;
        

        //calculate impact point on calorimeter
        
            float charge = track.user_index()/100;        
    
//         if (Bfield>0 && pT/fabs(charge)/(0.3*Bfield)*1000*2<1800)
//         {
// //             std::cout << "track with pT = " << pT << "  did not reach the calorimeter" << std::endl;
//             pfaCollection.push_back(track);
//             continue;
//         }
//         
        
        
            PseudoJet effectiveTrackEcal;
            PseudoJet effectiveTrackHcal;
            
            // sort calo hits by distance from impact point of the track
            std::vector<std::pair<PseudoJet, PseudoJet>> closeHits;
            std::vector<std::pair<PseudoJet, PseudoJet>> farHits;
            std::vector<std::pair<PseudoJet, PseudoJet>> sortedHits;
            std::cout <<"sorting calo hits by distance from track " << i_track << std::endl;
        
            if (Bfield>0)
            {
                TGraph* thisTraj = getEquivalentTrajectory (Bfield, track.px(), track.py(), track.pz(), charge, ecal_impact_radius);
                Double_t impact_x, impact_y;
                thisTraj->GetPoint(thisTraj->GetN()-1, impact_x, impact_y);
        
                float impact_phi = atan(impact_y/impact_x);
                if (impact_x<0. && impact_y <0.)   {impact_phi = impact_phi - M_PI;}
                if (impact_x<0. && impact_y >0.)   {impact_phi = M_PI + impact_phi;}        
                double impact_theta = M_PI- 2*atan(exp(-track.eta()));                
                float scale_p = 1./sqrt(impact_x*impact_x + impact_y*impact_y) * sqrt(pow(track.px(),2)+pow(track.py(),2));
                effectiveTrackEcal = PseudoJet(impact_x*scale_p, impact_y*scale_p, track.pz(), trueEne);
        
        
                thisTraj = getEquivalentTrajectory (Bfield, track.px(), track.py(), track.pz(), charge, hcal_impact_radius);        
                thisTraj->GetPoint(thisTraj->GetN()-1, impact_x, impact_y);
                impact_phi = atan(impact_y/impact_x);
                if (impact_x<0. && impact_y <0.)   {impact_phi = impact_phi - M_PI;}
                if (impact_x<0. && impact_y >0.)   {impact_phi = M_PI + impact_phi;}        
                impact_theta = M_PI- 2*atan(exp(-track.eta()));        
                scale_p = 1./sqrt(impact_x*impact_x + impact_y*impact_y) * sqrt(pow(track.px(),2)+pow(track.py(),2));
                effectiveTrackHcal = PseudoJet(impact_x*scale_p, impact_y*scale_p, track.pz(), trueEne);
            
                for (auto hit : leftCaloHits)
                {
                    PseudoJet scint_hit = hit.first;
                    float maxDeltaR = 0;
                    if      (scint_hit.user_index()==flag_JES) maxDeltaR = R_ECAL_cut[iR];
                    if      (scint_hit.user_index()==flag_JHS) maxDeltaR = R_HCAL_cut[iR];
        
                    float thisDD;
                    if      (scint_hit.user_index()==flag_JES) thisDD = scint_hit.delta_R(effectiveTrackEcal);
                    if      (scint_hit.user_index()==flag_JHS) thisDD = scint_hit.delta_R(effectiveTrackHcal);
                    
                    if (maxDeltaR<=0 || thisDD<maxDeltaR) 
                    {
                        closeHits.push_back(hit);
                    }
                    else
                    {
                        farHits.push_back(hit);
                    }
                }
                sortedHits = sorted_by_dd(closeHits, effectiveTrackEcal);//, R_ECAL_cut[iR], R_HCAL_cut[iR]);
            }
            else 
            {
                for (auto hit : leftCaloHits)
                {
                    PseudoJet scint_hit = hit.first;
                    float maxDeltaR = 0;
                    if      (scint_hit.user_index()==flag_JES) maxDeltaR = R_ECAL_cut[iR];
                    if      (scint_hit.user_index()==flag_JHS) maxDeltaR = R_HCAL_cut[iR];
        
                    float thisDD;
                    if      (scint_hit.user_index()==flag_JES) thisDD = scint_hit.delta_R(track);
                    if      (scint_hit.user_index()==flag_JHS) thisDD = scint_hit.delta_R(track);
                    
                    if (maxDeltaR<=0 || thisDD<maxDeltaR) 
                    {
                        closeHits.push_back(hit);
                    }
                    else
                    {
                        farHits.push_back(hit);
                    }
                }
                sortedHits = sorted_by_dd(closeHits, track);//, R_ECAL_cut[iR], R_HCAL_cut[iR]);
            }
        
//         std::cout << "n sorted hits found: " << sortedHits.size() << std::endl;
        
            //clear collection of leftover calo hits
            leftCaloHits = farHits;
            farHits.clear();
            
            for (auto hit : sortedHits)
            {
                PseudoJet scint_hit = hit.first;
                PseudoJet cher_hit  = hit.second;
//             float deltaR  = scint_hit.delta_R(track);
//             float deltaR;
//             if      (Bfield > 0 && scint_hit.user_index()==flag_JES) deltaR = scint_hit.delta_R(effectiveTrackEcal);
//             else if (Bfield > 0 && scint_hit.user_index()==flag_JHS) deltaR = scint_hit.delta_R(effectiveTrackHcal);
//             else if (Bfield == 0)                                    deltaR = scint_hit.delta_R(track);
            
//             float hit_theta = 2*atan(exp(-scint_hit.eta()));
//             float deltaR  = sqrt(pow(scint_hit.phi()-impact_phi,2) + pow(hit_theta-impact_theta,2));;
//             std::cout << "deltaR = " << hit.delta_R(track) << " :: deltaDD = "  << deltaDD << std:: endl;
                
                std::cout <<"i_track = " << i_track << " --> cycling over calo hit" << std::endl;
            
                double this_E_DRO;
                double this_S = scint_hit.E();
                double this_C = cher_hit.E();
                if      (scint_hit.user_index()==flag_JES) this_E_DRO = (scint_hit.E()-x_factor_ecal*cher_hit.E() )/(1-x_factor_ecal);
                else if (scint_hit.user_index()==flag_JHS) this_E_DRO = (scint_hit.E()-x_factor_hcal*cher_hit.E() )/(1-x_factor_hcal);
                
//             this_E_DRO = scint_hit.E();
            
//                 totCaloE[i_track] += (totCaloS_EC-x_factor_ecal*totCaloC_EC )/(1-x_factor_ecal) + (totCaloS_HC-x_factor_hcal*totCaloC_HC )/(1-x_factor_hcal);
//             totCaloE = totCaloS_EC+totCaloS_HC;
            
            if (totCaloE[i_track] < targetEne &&
//                 if (totCaloE[i_track] < targetEne+funcHcalRes->Eval(targetEne)*matchPFAcut*targetEne &&
                fabs(totCaloE[i_track]+this_E_DRO-targetEne) < fabs(totCaloE[i_track]-targetEne)   
//             if (current_tot_dro < targetEne &&                
//                 fabs(current_tot_dro+this_E_DRO-targetEne) < fabs(current_tot_dro-targetEne) &&
//                 ((scint_hit.user_index()==flag_JES  && deltaR < maxDeltaR_ECAL) ||
//                  (scint_hit.user_index()==flag_JHS  && deltaR < maxDeltaR_HCAL) )
                )
                {
//                 totCaloE += this_E_DRO;
                    if (scint_hit.user_index()==flag_JES)
                    {
                        totCaloS_EC[i_track] += this_S;
                        totCaloC_EC[i_track] += this_C;
                    }
                    else if (scint_hit.user_index()==flag_JHS)
                    {
                        totCaloS_HC[i_track] += this_S;
                        totCaloC_HC[i_track] += this_C;
                    }
                    matchedCaloHits[i_track].push_back(hit);                    
                }
                else
                {
                    leftCaloHits.push_back(hit);
                }
            }
        
            //update tot calo energy after DRO correction
            totCaloE[i_track] = (totCaloS_EC[i_track]-x_factor_ecal*totCaloC_EC[i_track] )/(1-x_factor_ecal) + (totCaloS_HC[i_track]-x_factor_hcal*totCaloC_HC[i_track] )/(1-x_factor_hcal);
            
            
            //check if track has enough energy 
//             if (totCaloE[i_track] > targetEne)
            {
                //check if full track is good enough to be swapped out --> no further radial iterations
                if (fabs(totCaloE[i_track]-targetEne)/targetEne < funcHcalRes->Eval(targetEne)*matchPFAcut ) 
                {
                    pfaCollection.push_back(track);
                
                    h1ResidualCharged->Fill((totCaloE[i_track]-trueEne)/trueEne);
                
                    trueTotCharged += trueEne;
                    recoTotCharged += totCaloE[i_track];
                    
                    swappedTrack++;
                    swappedTrackList.push_back(i_track);
                }
                //matching rejected, reinitialize track matching
//                 else
//                 {
//                     for (auto non_accepted_hit : matchedCaloHits[i_track])
//                     {
//                         leftCaloHits.push_back(non_accepted_hit);
//                     }
//                     matchedCaloHits[i_track].clear();          
//                     totCaloE[i_track] = 0;
//                     totCaloS_EC[i_track] = 0;
//                     totCaloC_EC[i_track] = 0;
//                     totCaloS_HC[i_track] = 0;
//                     totCaloC_HC[i_track] = 0;
//                 }
            }
            
            
            
/*            
            if (stopClustering && fabs(totCaloE[i_track]-targetEne)/targetEne > funcHcalRes->Eval(targetEne)*matchPFAcut) 
            {
                for (auto non_accepted_hit : matchedCaloHits[i_track])
                {
                    leftCaloHits.push_back(non_accepted_hit);
                }
            }*/
            
//         totCaloE[i_track] += totCaloS_EC+totCaloS_HC;
//         std::cout << "size sorted = " << sortedHits.size() << " :: size left =  " << leftCaloHits.size() << std::endl;
            
            i_track++;
        }
    }//end of cycle over radial integrations
        
    //final residual matching?
    //[...]
      
      
    //check if accept or reject track-hit matching
//     int i_track = 0;
    
/*
    for (auto track : sortedTracks)
    {
        
        float trueEne = track.E();
        float targetEne = trueEne*funcTotHadRawResponse->Eval(trueEne);
        
        if (totCaloE[i_track]>0.005)
        {
                
            //matching was good enough
            //             if (fabs(totCaloE-trueEne)/trueEne <funcHcalRes->Eval(trueEne)) 
            std::cout << "trueEne[" << i_track << "]: " << trueEne << " :: totCaloE["<<i_track<<"] = " << totCaloE[i_track] << std::endl;
            if (fabs(totCaloE[i_track]-targetEne)/targetEne < funcHcalRes->Eval(targetEne)*matchPFAcut) 
            {
                pfaCollection.push_back(track);
                
    //                 h2ScatterPFACharged->Fill(trueEne, totCaloE[i_track]/trueEne );
                h1ResidualCharged->Fill((totCaloE[i_track]-trueEne)/trueEne);
                
                trueTotCharged += trueEne;
                recoTotCharged += totCaloE[i_track];
                    
                swappedTrack++;
            }
                
            //matching rejected
            else
            {
                for (auto non_accepted_hit : matchedCaloHits[i_track])
                {
                    leftCaloHits.push_back(non_accepted_hit);
                }
            }
        }
    //         if (sortedHits.size()==0) pfaCollection.push_back(track);
        i_track++;
    }
    */
    
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







std::pair<std::vector<PseudoJet>,std::vector<std::pair<PseudoJet, PseudoJet>> > RunNeutralHitsCleanUp (std::pair<std::vector<PseudoJet>,std::vector<std::pair<PseudoJet, PseudoJet>> > pfaCollection, std::vector<PseudoJet> chargedTracks)
{
    
    std::vector<PseudoJet> cleanedPfaCollection;
    std::vector<std::pair<PseudoJet, PseudoJet>> cleanedNeutralHits;
    std::vector<std::pair<PseudoJet, PseudoJet>> neutralHCSeedCandidates;
    std::vector<std::pair<PseudoJet, PseudoJet>> neutralECSeedCandidates;
    std::vector<std::pair<PseudoJet, PseudoJet>> myEcSeeds;
    std::vector<std::pair<PseudoJet, PseudoJet>> myHcSeeds;
    
    float ec_seed_th = 0.1; //mip
    float hc_seed_th = 0.1;
    
    float maxDeltaRSeedEcal = 0.01;
    float maxDeltaRSeedHcal = 0.1;
    
    float ec_hc_R_match = 0.1;
    float ec_mc_R_match = 0.005;
    
    float ec_cluster_R = 0.2;
    float hc_cluster_R = 0.5;
    
    std::cout << "total calo hits in leftover collection: " << pfaCollection.second.size() << std::endl;
    
    for (auto pfaObj : pfaCollection.first)
    {
        if (fabs(pfaObj.user_index()) == 0 || fabs(pfaObj.user_index()) > 99 )
        {
            cleanedPfaCollection.push_back(pfaObj);
        }
    }
    
    for (auto hit : pfaCollection.second)
    {
        PseudoJet scint_hit = hit.first;
        PseudoJet cher_hit  = hit.second;
        
//             float deltaR  = scint_hit.delta_R(track);
//             float deltaR;
//             if      (Bfield > 0 && scint_hit.user_index()==flag_JES) deltaR = scint_hit.delta_R(effectiveTrackEcal);
//             else if (Bfield > 0 && scint_hit.user_index()==flag_JHS) deltaR = scint_hit.delta_R(effectiveTrackHcal);
        // find seeds
//         if (scint_hit.E()>seed_th)
        {
            if (scint_hit.user_index() == flag_JES && scint_hit.E()>ec_seed_th ) myEcSeeds.push_back(std::make_pair(scint_hit, cher_hit));
            if (scint_hit.user_index() == flag_JHS && scint_hit.E()>hc_seed_th) myHcSeeds.push_back(std::make_pair(scint_hit, cher_hit));
        }
    }
    
    std::cout << "ec seeds before cleanup: " << myEcSeeds.size() << std::endl;
    std::cout << "hc seeds before cleanup: " << myHcSeeds.size() << std::endl;
    
    
    //cleanup to find isolated seeds
    std::vector<std::pair<PseudoJet, PseudoJet>> myEcSeedsCleaned = pfaFindIsolatedSeeds(myEcSeeds, maxDeltaRSeedEcal);
    std::vector<std::pair<PseudoJet, PseudoJet>> myHcSeedsCleaned = pfaFindIsolatedSeeds(myHcSeeds, maxDeltaRSeedHcal);
//     
    std::cout << "ec seeds after cleanup: " << myEcSeedsCleaned.size() << std::endl;
    std::cout << "hc seeds after cleanup: " << myHcSeedsCleaned.size() << std::endl;
    
    
    //identify neutral seeds
    //HCAL seeds not matched to any ECAL seed 
    for (auto hc_seed : myHcSeedsCleaned)
    {
        bool unMatchedToEcal = true;
        PseudoJet hc_scint_hit = hc_seed.first;
        PseudoJet hc_cher_hit  = hc_seed.second;
        
        PseudoJet ec_scint_hit;
        PseudoJet ec_cher_hit;
        
        std::pair<PseudoJet,PseudoJet> ecal_seed_matched;
//         for (auto ec_seed : myEcSeedsCleaned)
        //I could run only on front crystal seeds to improve this
        for (auto ec_seed : myEcSeeds)
        {
            ec_scint_hit = ec_seed.first;
            ec_cher_hit  = ec_seed.second;
            float dd = ec_scint_hit.delta_R(hc_scint_hit);
            
            if (dd < ec_hc_R_match)
            {
                ecal_seed_matched = ec_seed;
                unMatchedToEcal = false;
            }
            
            break;
        }
        if (unMatchedToEcal)    //certainly a neutral seed?
        {
            neutralHCSeedCandidates.push_back(std::make_pair(hc_scint_hit, hc_cher_hit));
        }
        else
        {
            bool unMatchedToTrack = true;
            //or with corresponding ecal hit not matched to any charged track?
            for (auto track : chargedTracks)
            {
                float dd_track = ecal_seed_matched.first.delta_R(track);
                if (dd_track < ec_mc_R_match) //ecal hit not matched to any charged track then the hcal hit is from a neutral
                {
                     unMatchedToTrack = false;
                     break;
                }
            }
            if (unMatchedToTrack) neutralHCSeedCandidates.push_back(std::make_pair(hc_scint_hit, hc_cher_hit));
        }
    }
            
    
    
    //ecal seeds not matched to any charged track
    for (auto ec_seed : myEcSeedsCleaned)
    {
        PseudoJet ec_scint_hit = ec_seed.first;
        PseudoJet ec_cher_hit  = ec_seed.second;
        bool unMatchedToTrack = true;
        
        for (auto track : chargedTracks)
        {
            float dd_track = ec_seed.first.delta_R(track);
            if (dd_track < ec_mc_R_match) //ecal hit not matched to any charged track then the hcal hit is from a neutral
            {
                unMatchedToTrack = false;
                break;
            }
        }
        if (unMatchedToTrack) neutralECSeedCandidates.push_back(std::make_pair(ec_scint_hit, ec_cher_hit));
    }
    
    
    std::cout << "neutral ec seeds not matched to charged track: " << neutralECSeedCandidates.size() << std::endl;
    std::cout << "neutral hc seeds not matched to charged track: " << neutralHCSeedCandidates.size() << std::endl;
    
    
    //cluster hits around seeds
    int hits_removed  = 0;
    int hits_remained = 0;
    
    for (auto hit : pfaCollection.second)
    {
        PseudoJet scint_hit = hit.first;
        PseudoJet cher_hit  = hit.second;
        
        //assume hit is not matched to neutral
        bool hitMatchedToNeutral = false;
        
        if (scint_hit.user_index() == flag_JES)
        {
            for (auto ec_seed : neutralECSeedCandidates)
            {
                float dd_track = ec_seed.first.delta_R(scint_hit);
                if (dd_track < ec_cluster_R) //ecal hit not matched to any charged track then the hcal hit is from a neutral
                {
                    hitMatchedToNeutral = true;
                    break;
                }
            }
        }
        
        else if (scint_hit.user_index() == flag_JHS)
        {
            for (auto hc_seed : neutralHCSeedCandidates)
            {
                float dd_track = hc_seed.first.delta_R(scint_hit);
                if (dd_track < hc_cluster_R) //hcal hit not matched to any charged track then the hcal hit is from a neutral
                {
                    hitMatchedToNeutral = true;
                    break;
                }
            }
        }
        
        //if hit was not clustered into any seed it was probably a charged component than ignore it, otherwise add hit to the pfa and neutral collection 
        if (hitMatchedToNeutral)
        {

            cleanedNeutralHits.push_back(std::make_pair(scint_hit, cher_hit));                        
            
            cleanedPfaCollection.push_back(scint_hit);
            cleanedPfaCollection.push_back(cher_hit);
            
            hits_remained++;                                    
        }
        else
        {
            //             std::cout << " unclusteredHit, not added to pfa collection " << std::endl;
            hits_removed++;
        }
    }
    
    std::cout << "fraction of hits removed from collection: " << float(hits_removed) << " / " << pfaCollection.second.size() << " = " << float(hits_removed)/pfaCollection.second.size() << std::endl;
    std::cout << "fraction of hits remaining in collection: " << float(hits_remained) << " / " << pfaCollection.second.size() << " = " << float(hits_remained)/pfaCollection.second.size() << std::endl;
    
    return std::make_pair(cleanedPfaCollection, cleanedNeutralHits);
    
}



//***************************************************************************************************************//
//Parse the collection of ECAL hits to find photon clusters and remove corresponding "photon" hits from collection,
//It returns:
    // a collection with residual calorimeter hits after photons have been removed,
    // and a collection with all the photon hits
//***************************************************************************************************************//

std::pair<std::vector<std::pair<PseudoJet, PseudoJet>>,std::vector<std::pair<PseudoJet, PseudoJet>> > RunNeutralHitEcalCleaning (
std::vector<PseudoJet> chargedTracks, std::vector<std::pair<PseudoJet, PseudoJet>> ecalHitCollection, float Bfield, float maxDeltaRSeedEcal,
TH1F *hNECNeutralSeeds, TH2F *hNeutralSeedShowerShapeScint, TH2F *hNeutralSeedShowerShapeCher, TH2F *hNeutralSeedCSratio)
{
    
    std::vector<std::pair<PseudoJet, PseudoJet>> cleanedEcalHitCollection;
    std::vector<std::pair<PseudoJet, PseudoJet>> photonHits;
    std::vector<std::pair<PseudoJet, PseudoJet>> neutralECSeedCandidates;
    std::vector<std::pair<PseudoJet, PseudoJet>> myEcSeeds;
    
    float ec_seed_th = 0.08; //mip
    float maxDeltaRIsolatedSeed = 0.015;
    float ec_mc_R_match = 0.015; //radius to consider a photon not matched to a charged track
    float ec_cluster_R = 0.013; //maxDeltaRSeedEcal
    
    std::cout << "total calo hits in ecal collection: " << ecalHitCollection.size() << std::endl;
    
    //find ECAL seeds
    for (auto hit : ecalHitCollection)
    {
        PseudoJet scint_hit = hit.first;
        PseudoJet cher_hit  = hit.second;
        {
            if (scint_hit.user_index() == flag_JES && scint_hit.E()>ec_seed_th ) myEcSeeds.push_back(std::make_pair(scint_hit, cher_hit));
        }
    }
    std::cout << "ec seeds before cleanup: " << myEcSeeds.size() << std::endl;
    

    //cleanup to find isolated seeds
    std::vector<std::pair<PseudoJet, PseudoJet>> myEcSeedsCleaned = pfaFindIsolatedSeeds(myEcSeeds, maxDeltaRIsolatedSeed);
    std::cout << "ec seeds after cleanup: " << myEcSeedsCleaned.size() << std::endl;
    
    
    //ecal seeds not matched to any charged track
    for (auto ec_seed : myEcSeedsCleaned)
    {
        PseudoJet ec_scint_hit = ec_seed.first;
        PseudoJet ec_cher_hit  = ec_seed.second;
        bool unMatchedToTrack = true;
        
        for (auto track : chargedTracks)
        {
            //impact point of track on ECAL
            PseudoJet effectiveTrackEcal = track;

            if (Bfield>0)
            {
                float charge = track.user_index()/100;
                TGraph* thisTraj = getEquivalentTrajectory (Bfield, track.px(), track.py(), track.pz(), charge, ecal_impact_radius);
                Double_t impact_x, impact_y;
                thisTraj->GetPoint(thisTraj->GetN()-1, impact_x, impact_y);

                float impact_phi = atan(impact_y/impact_x);
                if (impact_x<0. && impact_y <0.)   {impact_phi = impact_phi - M_PI;}
                if (impact_x<0. && impact_y >0.)   {impact_phi = M_PI + impact_phi;}
                double impact_theta = M_PI- 2*atan(exp(-track.eta()));
                float scale_p = 1./sqrt(impact_x*impact_x + impact_y*impact_y) * sqrt(pow(track.px(),2)+pow(track.py(),2));
                effectiveTrackEcal = PseudoJet(impact_x*scale_p, impact_y*scale_p, track.pz(), track.E());
            }


            float dd_track = ec_seed.first.delta_R(effectiveTrackEcal);
            if (dd_track < ec_mc_R_match) //ecal hit not matched to any charged track then the hcal hit is from a neutral
            {
                unMatchedToTrack = false;
                break;
            }
        }
        if (unMatchedToTrack) neutralECSeedCandidates.push_back(std::make_pair(ec_scint_hit, ec_cher_hit));
    }
    
    
    std::cout << "potential neutral seeds (ec seeds not matched to any charged track): " << neutralECSeedCandidates.size() << std::endl;
    hNECNeutralSeeds->Fill(neutralECSeedCandidates.size());

    //calculate shower shape for seeds
    std::vector<std::pair<PseudoJet, PseudoJet>> photonSeedCandidates;
    for (auto ec_seed : neutralECSeedCandidates)
    {

        float seed_cluster_ene_S = 0;
        float seed_cluster_ene_C = 0;

        for (auto hit : ecalHitCollection)
        {
            PseudoJet scint_hit = hit.first;
            PseudoJet cher_hit = hit.second;

            float dd_seed = ec_seed.first.delta_R(scint_hit);
            if (dd_seed < ec_cluster_R) //hit within seed integration cone
            {
                seed_cluster_ene_S+=scint_hit.E();
                seed_cluster_ene_C+=cher_hit.E();
            }
        }

        float seed_shower_shape_S = ec_seed.first.E()/seed_cluster_ene_S;
        float seed_shower_shape_C = ec_seed.second.E()/seed_cluster_ene_C;
        hNeutralSeedShowerShapeScint->Fill(ec_seed.first.E(),seed_shower_shape_S);
        hNeutralSeedShowerShapeCher ->Fill(ec_seed.second.E(),seed_shower_shape_C);
        hNeutralSeedCSratio->Fill(ec_seed.first.E(),ec_seed.first.E()/ec_seed.second.E());

//         std::cout << "S = " << ec_seed.first.E() << " :: S/C = " << ec_seed.first.E()/ec_seed.second.E() << ":: S_tot/C_tot " << seed_cluster_ene_S/seed_cluster_ene_C<< " :: seed_shower_shape_S = " << seed_shower_shape_S << " :: seed_shower_shape_C = " << seed_shower_shape_C << std::endl;

        if ( (seed_shower_shape_S < 0.95 && seed_cluster_ene_S/seed_cluster_ene_C < 10 )
            || ec_seed.first.E() > 0.3 //include all non mips seeds
        ) photonSeedCandidates.push_back(std::make_pair(ec_seed.first, ec_seed.second));
    }

    std::cout << "potential photon seeds (passing photon shower shape cuts): " << photonSeedCandidates.size() << std::endl;

    
    //cluster hits around seeds to find if transverse shower shape is compatible with a photon
    int hits_remained  = 0;
    int photon_hits = 0;
    
    for (auto hit : ecalHitCollection)
    {
        PseudoJet scint_hit = hit.first;
        PseudoJet cher_hit  = hit.second;
        
        //assume hit is not matched to any neutral seed candidate
        bool hitMatchedPhotonSeed = false;
        
        if (scint_hit.user_index() == flag_JES)
        {
//             for (auto ec_seed : neutralECSeedCandidates)
            for (auto ec_seed : photonSeedCandidates)
            {
                float dd_seed = ec_seed.first.delta_R(scint_hit);
                if (dd_seed < ec_cluster_R) //ecal hit not matched to any charged track then the hcal hit is from a neutral
                {
                    hitMatchedPhotonSeed = true;
                    break;
                }
            }
        }
        
        
        //if hit was not clustered into any seed it was probably a charged component than ignore it, otherwise add hit to the pfa and neutral collection 
        if (hitMatchedPhotonSeed)
        {
            photonHits.push_back(std::make_pair(scint_hit, cher_hit));
            photon_hits++;
        }
        else
        {
            cleanedEcalHitCollection.push_back(std::make_pair(scint_hit, cher_hit));
            hits_remained++;
        }
    }
    
    std::cout << "fraction of photon hits removed from collection: " << photon_hits << " / " << ecalHitCollection.size() << " = " << float(photon_hits)/ecalHitCollection.size() << std::endl;
    std::cout << "fraction of hits remained in pfa collection: " << hits_remained << " / " << ecalHitCollection.size() << " = " << float(hits_remained)/ecalHitCollection.size() << std::endl;
    
    return std::make_pair(cleanedEcalHitCollection, photonHits);
    
}





std::vector<std::pair<PseudoJet, PseudoJet>> pfaFindIsolatedSeeds (std::vector<std::pair<PseudoJet, PseudoJet>> allSeeds, float maxDeltaR)
{

    std::vector<std::pair<PseudoJet, PseudoJet>> isolatedSeeds;
    
    for (auto seed_i : allSeeds)
    {
        PseudoJet hit_i = seed_i.first;
        PseudoJet cher_hit_i = seed_i.second;
        bool maxIsolatedHit = true;
        
        for (auto seed_j : allSeeds)
        {
            PseudoJet hit_j = seed_j.first;
            if (hit_j == hit_i) continue;
            
            float dd = hit_i.delta_R(hit_j);
            if (dd < maxDeltaR && hit_i.E()<hit_j.E())
            {
                maxIsolatedHit = false;
            }
        }
        
        if (maxIsolatedHit)  isolatedSeeds.push_back(std::make_pair(hit_i, cher_hit_i));
    }
    
    return isolatedSeeds;
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


/*
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


*/

