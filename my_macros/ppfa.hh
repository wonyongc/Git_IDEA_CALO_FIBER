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


using namespace fastjet; 

std::vector<PseudoJet> RunProtoPFA (std::vector<PseudoJet> chargedTracks, std::vector<std::pair<PseudoJet, PseudoJet>> hitsForJet, 
                                    float my_x_factor_ecal, float my_x_factor_hcal,
                                    TH1F* check1Histo1D, TH1F *check2Histo1D, TH1F *check3Histo1D);

std::vector<std::pair<PseudoJet, PseudoJet>> sorted_by_dd (std::vector<std::pair<PseudoJet,PseudoJet>> myJets, PseudoJet refJet);

Double_t rms90(TH1F *h) ;
