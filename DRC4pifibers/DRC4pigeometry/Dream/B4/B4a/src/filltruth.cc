#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "HepMC/IO_GenEvent.h"
#include "HepMC/GenEvent.h"
#include "HepMC/Units.h"
#include <cmath>
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "filltruth.h"
#include <iomanip>
#include <string>
double GeV=1000.;
double cfac=1.;
//
//  setup for particle charge
//
TDatabasePDG *fPDG;
TParticlePDG *pdgParticle;

using namespace std;

filltruth::filltruth() {
  fPDG = TDatabasePDG::Instance();
  cout << " create filltruth " << endl;
}

//this is the initalization function, put any new branches you want to create here, remember to declare them in the filltruth.h file as well 
void filltruth::book_tuple(string outfile) {
  std::cout << "Preparing Root Tree" << endl;
  fileout = new TFile(outfile.c_str(),"RECREATE");
  tree = new TTree("truth","truth tree");
  isTree = true;
  
  tree->Branch("mcs_n",&mcs_n,"mcs_n/i");
  tree->Branch("mcs_E",&mcs_E);
  tree->Branch("mcs_pt",&mcs_pt);
  tree->Branch("mcs_m",&mcs_m);
  tree->Branch("mcs_eta",&mcs_eta);
  tree->Branch("mcs_phi",&mcs_phi);
  tree ->Branch("mcs_status", &mcs_status);
  tree->Branch("mcs_barcode",&mcs_barcode);
  tree->Branch("mcs_pdgId", &mcs_pdgId);
  tree->Branch("mcs_charge", &mcs_charge);
  tree->Branch("mcs_vx_x", &mcs_vx_x);
  tree->Branch("mcs_vx_y", &mcs_vx_y);
  tree->Branch("mcs_vx_z", &mcs_vx_z); 
  
}
void filltruth::fill_tuple(HepMC::GenEvent* evt) {

   AnalyseParticles(evt);
//  Clear vectors 
   mcs_n=0;
   mcs_E.clear();
   mcs_pt.clear();
   mcs_m.clear();
   mcs_eta.clear();
   mcs_phi.clear();
   mcs_status.clear();
   mcs_barcode.clear();
   mcs_pdgId.clear();
   mcs_charge.clear();
   mcs_vx_x.clear();
   mcs_vx_y.clear();
   mcs_vx_z.clear();  
//
// Conversion factor if original HEPMC was in GEV
//   HepMC::Units::MomentumUnit mun=evt->momentum_unit();
//   if(mun==1)cfac=1000.;   
   fillstable(index_to_particle);
   tree->Fill();
}
//write the root tree to a file (hpp.root by default)
void filltruth::write_tuple() { 
  if (isTree){
  tree->GetCurrentFile();
  tree->Write();
  fileout->Close();
  cout << "A root tree has been written to a file" << endl;}
}
void filltruth::fillstable( std::vector<HepMC::GenParticle*> evt ) {

  int nstab=0;
  int nbpar=0;
//  cout << "****************************" << endl;
  for(int n=0; n<evt.size(); n++) {
    HepMC::GenParticle* ipan=evt.at(n);
    int status=ipan->status();
    int idp=ipan->pdg_id();
    int aidp=abs(ipan->pdg_id());
    if(status==1) {
      nstab++;
      mcs_E.push_back(ipan->momentum().e()*cfac );
      mcs_pt.push_back(ipan->momentum().perp()*cfac );
      mcs_m.push_back(ipan->momentum().m()*cfac );
      mcs_eta.push_back(ipan->momentum().eta() );
      mcs_phi.push_back(ipan->momentum().phi() );
      mcs_status.push_back(ipan->status() );
      int id=ipan->pdg_id();
      pdgParticle = fPDG->GetParticle(id);
      float charge=pdgParticle ? int(pdgParticle->Charge()/3.0) : -999;
      mcs_pdgId.push_back(id);
      mcs_charge.push_back(charge);
      mcs_vx_x.push_back(0.);
      mcs_vx_y.push_back(0.);
      mcs_vx_z.push_back(0.);
      mcs_barcode.push_back(ipan->barcode());
    }
  }
  mcs_n=nstab;
}
//-------------------------------------------------------------------------
int filltruth::find_in_map( const std::map<HepMC::GenParticle*,int>& m, HepMC::GenParticle *p) const
{
  std::map<HepMC::GenParticle*,int>::const_iterator iter = m.find(p);
  return (iter == m.end()) ? 0 : iter->second;
}

//--------------------------------------------------------------------------
void  filltruth::ReadStats(const HepMC::GenEvent* evt) {

 unsigned int particle_counter=0;
 index_to_particle.clear();

 HepMC::GenEvent::vertex_const_iterator v;

 // std::cout << " evt->size() = " << evt.size() <<std::endl;
 for (v = evt->vertices_begin(); v != evt->vertices_end(); ++v ) {
// making a list of incoming particles of the vertices
// so that the mother indices in HEPEVT can be filled properly
     HepMC::GenVertex::particles_out_const_iterator p1;

     for (p1 = (*v)->particles_in_const_begin();p1 != (*v)->particles_in_const_end(); ++p1 ) {
	 ++particle_counter;
         index_to_particle.push_back(*p1);
	 particle_to_index[*p1] = particle_counter-1;

     }   
// daughters are entered only if they aren't a mother of
// another vertex
     HepMC::GenVertex::particles_out_const_iterator p2;

     for (p2 = (*v)->particles_out_const_begin();p2 != (*v)->particles_out_const_end(); ++p2) {
       if (!(*p2)->end_vertex()) {
         ++particle_counter;
         index_to_particle.push_back(*p2);
         particle_to_index[*p2] = particle_counter-1;
       }
     }
   } 
   
}

void filltruth::getVecFromTuple(const HepMC::GenEvent* evt, vector<int>& movec, vector<int>& davec, int j) const {
  davec.clear();
  movec.clear();
  if (!evt) {
    cout <<  "HepMCFileReader: Got no event :-(  Game over already  ?" <<endl; 
  } 
  else {
    HepMC::GenParticle* ipaj=index_to_particle[j];
    if ( ipaj->production_vertex() ) { 
      HepMC::GenVertex::particle_iterator ic;
      for (ic = ipaj->production_vertex()->particles_begin(HepMC::parents);
           ic != ipaj->production_vertex()->particles_end(HepMC::parents); ++ic) {
           int indp=find_in_map( particle_to_index,*ic);
           movec.push_back(indp);
      }
    }
    if (ipaj->end_vertex()) {
      HepMC::GenVertex::particle_iterator ic;
      for (ic = ipaj->end_vertex()->particles_begin(HepMC::children);
           ic != ipaj->end_vertex()->particles_end(HepMC::children); ++ic) {
           int indp=find_in_map( particle_to_index,*ic);
           davec.push_back(indp);
      }
    } 
  }
//  std::cout << " nmoth " << movec.size() << " nch " << davec.size() << std::endl; 
} 
void filltruth::AnalyseParticles(const HepMC::GenEvent* evt)
{
//
//  fills vector of parents and children
//
  ReadStats(evt);
  parents.clear();
  children.clear();
  for(int n=0; n<index_to_particle.size(); n++) {
    vector<int> movec;
    vector<int> davec;
    HepMC::GenParticle* ipan=index_to_particle[n];
    getVecFromTuple(evt, movec, davec, n);
    int nmo=movec.size();    
    int nch=davec.size();    
    vector<int> ch_index;
    vector<int> par_index;
    for (Int_t j=0; j<nch;j++) {
      ch_index.push_back(davec.at(j));
    }
    for (Int_t j=0; j<nmo;j++) {
      par_index.push_back(movec.at(j));
    }
    parents.push_back(par_index);
    children.push_back(ch_index);
  }
}
