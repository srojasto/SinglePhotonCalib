
//#ifndef SPECALCULATION_H
//#define SPECALCULATION_H

struct outValues{
    Double_t fraction = 0;
    Double_t belowTrs = 0;
    Double_t overTrs = 0;
    Double_t totalN = 0;
};

struct SPEValues{
  Double_t occupancy = 0;
  Double_t mean =0;
  Double_t variance =0;
  Double_t stdDev = 0;
  Double_t meanUncertainty = 0;
  Double_t VarianceUncertainty = 0;
};

outValues CalculateFraction(TH1* ,Double_t ,Bool_t );
SPEValues CalculateSPE(TH1*, TH1*, Double_t, Bool_t );
Double_t NormHisto(TH1*, TH1*);
TH1* GetQHisto(const Char_t* fileRoot, const Char_t* channel);
TH1* GetAmpHisto(const Char_t* fileRoot, const Char_t* channel);
