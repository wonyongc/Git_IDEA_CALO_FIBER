#ifndef B4MyMaterials_hh
#define B4MyMaterials_hh

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"



class MyMaterials
{
private:
  
public:
  MyMaterials();
  ~MyMaterials();
  
  
  static G4Material* OpticalGrease();
  static G4Material* OpticalGrease155();
  static G4Material* MeltMount168();
  
  static G4Material* LYSO();
  static G4Material* PWO();
  static G4Material* BGO();
  static G4Material* CsI();

  
  static G4double fromNmToEv(G4double wavelength);
  static G4double fromEvToNm(G4double energy);
  static G4double CalculateSellmeier(int size, G4double indexZero, G4double *nVec, G4double *lVec, G4double wavelength);
};

#endif
