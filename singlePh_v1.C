// -*- C++ -*-
// $Id$
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TStyle.h"

struct outValues{
    Double_t fraction = 0;
    Double_t belowTrs = 0;
    Double_t overTrs = 0;
    Double_t totalN = 0;
};


outValues CalculateFraction(TH1* histo, Double_t threshold, Bool_t print = kTRUE)
{
  // Calculate occupancy using the scaled histograms
  //test second branch modification
  outValues result;

  Int_t TrsBin = histo ->  FindBin(threshold);
  Int_t MinBin = histo ->  GetMinimumBin();
  Int_t MaxBin = histo ->  GetMaximumBin();

  result.belowTrs = histo->Integral(MinBin,TrsBin);
  result.overTrs = histo->Integral(TrsBin+1,MaxBin);
  result.totalN = histo->Integral();
  result.fraction = result.belowTrs/result.totalN;

  if(print){
    cout << "Threshold bin = " << histo ->  FindBin(threshold);
    cout << "\tMinimum bin = " << histo ->  GetMinimumBin();
    cout << "\tMaximum bin = " << histo ->  GetMaximumBin() << endl;
    cout << "belowTrs integral = " << result.belowTrs;
    cout << "\toverTrs integral = " << result.overTrs << endl;
    cout << "belowTrs + overTrs = " << result.belowTrs + result.overTrs << " while the totalN = " <<  result.totalN << endl;
    cout << "fraction = " << result.fraction <<"\n\n";
  }

  return result;
}

void singlePh_v1(const Char_t* fileRoot="results.root",const Char_t* filePedestal="pedestal.root"){

  TString aFile = fileRoot;
  TFile *f = TFile::Open(aFile);
  TTree *TESig = (TTree*)f->Get("RawDataTree");
  // TE->Print();
  // TE->Show(10);

  TString pedFile = filePedestal;
  TFile *fped = TFile::Open(pedFile);
  TTree * TEped = (TTree*)fped->Get("RawDataTree");

  gStyle->SetOptStat(0);

  TESig->Draw("main_FDDref.Charge_GATE>>h1(1200,-300,0)", "", "GOFF");
  TH1 *h1 = TESig->GetHistogram();
  h1->SetTitle("Charge distribution;Charge;Entries");

  Double_t w = 1400;
  Double_t h = 1000;
  Double_t threshold = -250;
  TCanvas * c1 = new TCanvas("c1","Signal charge", w, h);
  c1->SetLogy();
  h1->Draw();


  TEped->Draw("main_FDDref.Charge_GATE>>h1Ped(1200,-300,0)", "", "GOFF");
  TH1 *h1Ped = TEped->GetHistogram();
  h1Ped->SetTitle("Charge distribution;Charge;Entries");

  cout << "before scaling" << endl;

  // Calculation of fraction of the defined threshold
  outValues blankResult = CalculateFraction(h1Ped, threshold);
  outValues signalResult = CalculateFraction(h1, threshold);



  // Getting the maximum bin and factor for normalization
  Int_t binmaxSig = h1->GetMaximumBin();
  Int_t nSig = h1->GetBinContent(binmaxSig);
  cout << "binmax="<< binmaxSig << endl;
  cout << "xbin="<< nSig << endl;

  Int_t binmax = h1Ped->GetMaximumBin();
  Int_t n = h1Ped->GetBinContent(binmax);
  cout << "binmax="<< binmax << endl;
  cout << "xbin="<< n << endl;

  // Normalization of blank data respect to laser data set
  Double_t NormFactor = float(nSig)/float(n);
  cout << "NormFactor="<< NormFactor << endl;
  h1Ped->Scale(NormFactor);

  h1Ped->SetLineColor(kRed);
  h1Ped->SetLineStyle(1);

  // Printing histograms in the same canvas
  TCanvas * cped = new TCanvas("cped","Pedestal", w, h);
  cped->SetLogy();
  h1->Draw();
  h1Ped->Draw("B SAME");

  // Obtaining the mean and stdDev of the distributions
  Double_t pedMean = h1Ped -> GetMean();
  Double_t pedMeanError = h1Ped -> GetMeanError();
  Double_t pedStdDev = h1Ped -> GetStdDev();
  Double_t pedStdDevError = h1Ped -> GetStdDevError();

  Double_t signalMean = h1 -> GetMean();
  Double_t signalMeanError = h1 -> GetMeanError();
  Double_t signalStdDev = h1 -> GetStdDev();
  Double_t signalStdDevError = h1 -> GetStdDevError();


  TLine *pedLine = new TLine(pedMean ,0 ,pedMean,nSig+nSig*1);
  pedLine -> SetLineColor(kRed);
  pedLine -> SetLineWidth(2);
  pedLine -> SetLineStyle(2);
  pedLine -> Draw();

  TLine *signalLine = new TLine(signalMean ,0 ,signalMean,nSig+nSig*1);
  signalLine -> SetLineColor(kBlue);
  signalLine -> SetLineWidth(2);
  signalLine -> SetLineStyle(2);
  signalLine -> Draw();

  TLatex pedLtx;
  pedLtx.SetTextSize(0.025);
  pedLtx.SetTextAlign(13);  //align at top
  pedLtx.SetTextColor(kRed);  //align at top
  pedLtx.DrawLatexNDC(.3,.85,"Background");
  pedLtx.DrawLatexNDC(.3,.8,TString::Format("E[B] = %0.2f #pm %0.2f", pedMean, pedMeanError));
  pedLtx.DrawLatexNDC(.3,.75,TString::Format("V[B] = %0.2f #pm %0.2f", pedStdDev, pedStdDevError));


  TLatex signalLtx;
  signalLtx.SetTextSize(0.025);
  signalLtx.SetTextAlign(13);  //align at top
  signalLtx.SetTextColor(kBlue);  //align at top
  signalLtx.DrawLatexNDC(.5,.85,"Total");
  signalLtx.DrawLatexNDC(.5,.8,TString::Format("E[T] = %0.2f #pm %0.2f", signalMean, signalMeanError));
  signalLtx.DrawLatexNDC(.5,.75,TString::Format("V[T] = %0.2f #pm %0.2f", signalStdDev, signalStdDevError));

  // Calculation of ocupancy
  // nf = number of events falling below the charge threshold cut in blank data

  // Plots to checck the occupancy and fraction behavior at different thresholds
  TGraph *grOccupancy = new TGraph();
  grOccupancy -> SetTitle("Ocupancy plot;Threshold (mV);Ocupancy");

  TGraph *grFraction = new TGraph();
  grFraction -> SetTitle("Fraction plot;Threshold (mV);Fraction");

  TGraph2D *dt = new TGraph2D();


  // Calculate occupancy using the scaled histograms
  cout << "after scaling" << endl;
  Int_t Index = 0;
  for (Double_t i = -290; i < -100; i += 1, Index++){

    outValues afblankResult = CalculateFraction(h1Ped, i, kFALSE);
    outValues afsignalResult = CalculateFraction(h1, i, kFALSE);
    Double_t occupancy;
    afsignalResult.belowTrs != 0? occupancy = -TMath::Log(afsignalResult.belowTrs/(afblankResult.fraction*afsignalResult.totalN)) : 0;

    grOccupancy -> SetPoint(Index, i, occupancy);
    grFraction-> SetPoint(Index, i, afblankResult.fraction);
    dt -> SetPoint(Index, i, occupancy, afblankResult.fraction);

    // cout << Index;
    // cout << ") thrs = " << i << "\toccupancy = " << occupancy;
    // cout << "\tf = " << afblankResult.fraction << endl;
  }

  TCanvas * cThreshold = new TCanvas("cThreshold"," Threshold", w, h);
  cThreshold -> Divide(1,2);
  cThreshold -> cd(1);
  grOccupancy -> Draw("AL*");
  cThreshold -> cd(2);
  grFraction -> Draw("AL*");

  TCanvas * cg2 = new TCanvas("cg2"," Threshold, occupancy and fraction", w, h);
  gStyle->SetPalette(1);
  dt->Draw("surf1");

	//TE->StartViewer();

}
