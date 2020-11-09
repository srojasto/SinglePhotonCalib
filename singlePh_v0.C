// -*- C++ -*-
// $Id$
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TStyle.h"

double CalculateFraction(TH1 * histo,Double_t threshold)
{
  // Calculate occupancy using the scaled histograms
  Int_t blankTrsBin = histo ->  FindBin(threshold);
  Int_t blankMinBin = histo ->  GetMinimumBin();
  Int_t blankMaxBin = histo ->  GetMaximumBin();
  cout << "Threshold bin = " << histo ->  FindBin(threshold);
  cout << "\tMinimum bin = " << histo ->  GetMinimumBin();
  cout << "\tMaximum bin = " << histo ->  GetMaximumBin() << endl;

  Double_t nfIntegral = histo->Integral(blankMinBin,blankTrsBin);
  Double_t ovfIntegral = histo->Integral(blankTrsBin+1,blankMaxBin);
  cout << "nf integral = " << nfIntegral;
  cout << "\tovf integral = " << ovfIntegral << endl;

  Double_t Integral = histo->Integral();
  cout << "integral = " << Integral << endl;
  cout << "nf + ovf = " << nfIntegral + ovfIntegral << " while the integral = " << Integral;
  cout << endl << endl;
  return nfIntegral/Integral;
}

void singlePh_v0(const Char_t* fileRoot="results.root",const Char_t* filePedestal="pedestal.root"){

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
  TCanvas * c1 = new TCanvas("c1","Signal charge", w, h);
  c1->SetLogy();
  h1->Draw();


  TEped->Draw("main_FDDref.Charge_GATE>>h1Ped(1200,-300,0)", "", "GOFF");
  TH1 *h1Ped = TEped->GetHistogram();
  h1Ped->SetTitle("Charge distribution;Charge;Entries");


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

  TCanvas * cped = new TCanvas("cped","Pedestal", w, h);
  cped->SetLogy();
  h1->Draw();
  h1Ped->Draw("B SAME");

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

  Double_t threshold = -250;
  TString trsCondition = TString::Format("main_FDDref.Charge_GATE < %f",threshold);

  Int_t nf = TEped->Draw("main_FDDref.Charge_GATE>>hAux(1200,-300,0)", trsCondition, "GOFF");
  cout << "nf = " << nf << endl;


  // At = number of events falling below the charge threshold cut in laser induces PE
  Int_t At = TESig->Draw("main_FDDref.Charge_GATE>>hAux(1200,-300,0)", trsCondition, "GOFF");
  cout << "At = " << At << endl;


  // Checking total number of entries of pedestal and signal run
  Int_t nBlank = TEped->Draw("main_FDDref.Charge_GATE>>hAux(1200,-300,0)","", "GOFF");
  cout << "nBlank = " << nBlank << endl;

  Double_t Nzero = double(At)*nBlank/double(nf);
  cout << "Nzero = " << Nzero << endl;

  Int_t nSignal = TESig->Draw("main_FDDref.Charge_GATE>>hAux(1200,-300,0)", "", "GOFF");
  cout << "nSignal = " << nSignal << endl;

  // Calculation of occupancy
  Double_t occupancy = -TMath::Log((Nzero/double(nSignal)));
  cout << "occupancy = " << occupancy << endl << endl;

  // Calculate occupancy using the scaled histograms
  Int_t blankTrsBin = h1Ped ->  FindBin(threshold);
  Int_t blankMinBin = h1Ped ->  GetMinimumBin();
  Int_t blankMaxBin = h1Ped ->  GetMaximumBin();
  cout << "Values for blank data:" << endl;
  cout << "Blank Threshold bin = " << h1Ped ->  FindBin(threshold);
  cout << "\tBlank Minimum bin = " << h1Ped ->  GetMinimumBin();
  cout << "\tBlank Maximum bin = " << h1Ped ->  GetMaximumBin() << endl;

  Double_t nfIntegral = h1Ped->Integral(blankMinBin,blankTrsBin);
  Double_t ovfIntegral = h1Ped->Integral(blankTrsBin+1,blankMaxBin);
  cout << "nf integral = " << nfIntegral;
  cout << "\tovf integral = " << ovfIntegral << endl;

  Double_t Integral = h1Ped->Integral();
  cout << "integral = " << Integral << endl;
  cout << "sum of the parts = " << nfIntegral + ovfIntegral << " while the integral = " << Integral;
  cout << endl << endl;


  Int_t sigTrsBin = h1 ->  FindBin(threshold);
  Int_t sigMinBin = h1 ->  GetMinimumBin();
  Int_t sigMaxBin = h1 ->  GetMaximumBin();
  cout << "Values for signal data:" << endl;
  cout << "Blank Threshold bin = " << h1 ->  FindBin(threshold);
  cout << "\tBlank Minimum bin = " << h1 ->  GetMinimumBin();
  cout << "\tBlank Maximum bin = " << h1 ->  GetMaximumBin() << endl;

  Double_t ATIntegral = h1->Integral(sigMinBin,sigTrsBin);
  Double_t ATovIntegral = h1->Integral(sigTrsBin+1,sigMaxBin);
  cout << "AT integral = " << ATIntegral;
  cout << "\tATov integral = " << ATovIntegral << endl;

  Double_t sigIntegral = h1->Integral();
  cout << "signal integral = " << sigIntegral << endl;
  cout << "sum of the parts = " << ATIntegral + ATovIntegral << " while the integral = " << sigIntegral << endl;

  cout << endl << "Calculation of scaled ocupancy:" << endl;
  Double_t scaledNzero = ATIntegral*Integral/nfIntegral;
  cout << "scaledNzero = " << scaledNzero << endl;

  Double_t scaledOccupancy = -TMath::Log((scaledNzero/double(nSignal)));
  cout << "scaledOccupancy = " << scaledOccupancy << endl << endl;

  cout << "fraction function test: " << CalculateFraction(h1, threshold) << endl;

	//TE->StartViewer();

}
