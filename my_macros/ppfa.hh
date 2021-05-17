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


using namespace fastjet; 

std::vector<PseudoJet> RunProtoPFA (std::vector<PseudoJet> chargedTracks, std::vector<PseudoJet> hitsForJet);

std::vector<PseudoJet> sorted_by_dd (std::vector<PseudoJet> myJets, PseudoJet refJet);

