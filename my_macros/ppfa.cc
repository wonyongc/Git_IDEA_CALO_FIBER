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
// #include "SCEPCal_GeometryHelper.hh"


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



