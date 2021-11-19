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
#include "fastjet/ClusterSequence.hh"
#include "TRandom.h"
#include "TF1.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TGraph.h"

// #include "SCEPCal_GeometryHelper.hh"

using namespace fastjet; 

// std::vector<PseudoJet> RunProtoPFA
std::pair<std::vector<PseudoJet>,std::vector<std::pair<PseudoJet, PseudoJet>> > RunProtoPFA(std::vector<PseudoJet> chargedTracks, std::vector<std::pair<PseudoJet, PseudoJet>> hitsForJet, 
                                    float my_x_factor_ecal, float my_x_factor_hcal, float Bfield, float matchPFAcut, bool DRO_ON,
                                    TH1F* check1Histo1D, TH1F *check2Histo1D, TH1F *check3Histo1D);

std::vector<PseudoJet> RunProtoPFA_Iterative (std::vector<PseudoJet> chargedTracks, std::vector<std::pair<PseudoJet, PseudoJet>> hitsForJet, 
                                    float my_x_factor_ecal, float my_x_factor_hcal, float Bfield, float matchPFAcut, bool DRO_ON,
                                    TH1F* check1Histo1D, TH1F *check2Histo1D, TH1F *check3Histo1D);

std::vector<std::pair<PseudoJet, PseudoJet>> sorted_by_dd (std::vector<std::pair<PseudoJet,PseudoJet>> myJets, PseudoJet refJet);//, float Rcut_ECAL, float Rcut_HCAL);



std::pair<std::vector<PseudoJet>,std::vector<std::pair<PseudoJet, PseudoJet>> >  RunNeutralHitsCleanUp (std::pair<std::vector<PseudoJet>,std::vector<std::pair<PseudoJet, PseudoJet>> > pfaCollection, std::vector<PseudoJet> chargedTracks);

std::pair<std::vector<std::pair<PseudoJet, PseudoJet>>,std::vector<std::pair<PseudoJet, PseudoJet>> >  RunNeutralHitEcalCleaning
(std::vector<PseudoJet> chargedTracks, std::vector<std::pair<PseudoJet, PseudoJet>> allEcalHits, float Bfield, float maxDeltaRSeedEcal, TH1F *hNECNeutralSeeds);


std::vector<std::pair<PseudoJet, PseudoJet>> pfaFindIsolatedSeeds(std::vector<std::pair<PseudoJet, PseudoJet>>allSeeds, float maxDeltaRSeed);


Double_t rms90(TH1F *h) ;

// TGraph * getEquivalentTrajectory (float B, float px, float py, float pT, float charge, float maxR);

