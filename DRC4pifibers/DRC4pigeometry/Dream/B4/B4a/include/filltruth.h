#ifndef filltruth_H
#define filltruth_H
//
// This is the declaration of the filltruth class.
//
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include "TBranch.h"

using namespace std;

class filltruth {

public:

  filltruth();
  ~filltruth();


  TTree * tree;
  TFile * fileout;

  //
  //  Ntuple with only stable particles to be used
  //  for object reconstruction
  //  it does not include include info on partons or 
  //  parent-child relationships
  //
  Int_t           mcs_n;
  vector<float>   mcs_E;
  vector<float>   mcs_pt;
  vector<float>   mcs_m;
  vector<float>   mcs_eta;
  vector<float>   mcs_phi;
  vector<int>     mcs_status;
  vector<int>     mcs_barcode;
  vector<int>     mcs_pdgId;
  vector<float>   mcs_charge;
  vector<float>   mcs_vx_x;
  vector<float>   mcs_vx_y;
  vector<float>   mcs_vx_z;
  bool isTree = false;
  //
  //  Methods for handling ntuple
  // 
  void book_tuple(std::string outputfile);
  void fill_tuple(HepMC::GenEvent* evt);
  void write_tuple();
  
  //
  //  Methods for converting HepMC into ntuple format
  //  Partially inspired by Delphes
  //
  void fillstable( std::vector<HepMC::GenParticle*> evt);
  void AnalyseParticles(const HepMC::GenEvent* evt);
  void getVecFromTuple(const HepMC::GenEvent* evt, vector<int>& movec, vector<int>&  davec, int j) const;
  void ReadStats(const HepMC::GenEvent* evt);
  int find_in_map(const std::map<HepMC::GenParticle*,int>& m,HepMC::GenParticle *p) const;
  typedef std::vector<HepMC::GenParticle*>::const_iterator  particles_in_const_iterator;
  //
  //  Intermediate service vectors 
  //
  std::vector<HepMC::GenParticle*> index_to_particle;
  std::map<HepMC::GenParticle*,int> particle_to_index;
  vector<vector<int> > parents;
  vector<vector<int> > children;

};
#endif
