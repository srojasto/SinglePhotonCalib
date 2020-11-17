// -*- C++ -*-
// $Id$
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLine.h"
#include "TMath.h"
#include <iostream>
#include "TLatex.h"
#include "TGraph.h"

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
  Double_t stdDevUncertainty = 0;
};

outValues CalculateFraction(TH1* ,Double_t ,Bool_t );
SPEValues CalculateSPE(TH1* h1Signal, TH1* h1Blank, Double_t threshold, Bool_t print = kTRUE);


void singlePh_v1(const Char_t* fileRoot="results.root",const Char_t* filePedestal="pedestal.root"){

  TString aFile = fileRoot;
  TFile *f = TFile::Open(aFile);
  TTree *TELaser = (TTree*)f->Get("RawDataTree");
  // TE->Print();
  // TE->Show(10);

  TString pedFile = filePedestal;
  TFile *fped = TFile::Open(pedFile);
  TTree * TEped = (TTree*)fped->Get("RawDataTree");

  gStyle->SetOptStat(0);

  TELaser->Draw("main_FDDref.Amplitude>>h1LaserAmp(1200,0,1200)", "", "GOFF");
  TH1 *h1LaserAmp = TELaser->GetHistogram();
  h1LaserAmp->SetTitle("Amplitude distribution Laser Induced;Amplitude (mV);Entries");

  TELaser->Draw("main_FDDref.Charge_GATE:main_FDDref.Amplitude>>h2AmpCharge(1200,0,1200, 1200,-300,0)", "", "GOFF");
  TH1 *h2AmpCharge = TELaser->GetHistogram();
  h2AmpCharge->SetTitle("Amplitude distribution Laser Induced;Amplitude (mV);Charge (fC);Entries");

  TEped->Draw("main_FDDref.Charge_GATE:main_FDDref.Amplitude>>h2BlankAmpCharge(1200,0,1200, 1200,-300,0)", "", "GOFF");
  TH1 *h2BlankAmpCharge = TEped->GetHistogram();
  h2BlankAmpCharge->SetTitle("Amplitude distribution Blank;Amplitude (mV);Charge (fC);Entries");

  TELaser->Draw("main_FDDref.Charge_GATE>>h1(1200,-300,0)", "", "GOFF");
  //TELaser->Draw("main_FDDref.Charge_GATE>>h1(1200,-300,0)", " main_FDDref.Amplitude < 150 && main_FDDref.Amplitude > 50 ", "GOFF");
  TH1 *h1 = TELaser->GetHistogram();
  h1->SetTitle("Charge distribution;Charge (fC);Entries");

  Double_t w = 1400;
  Double_t h = 1000;
  Double_t threshold = -250;
  TCanvas * c1 = new TCanvas("c1","Signal charge", w, h);
  c1 -> Divide(2,2);
  c1 -> cd(1)->SetLogy();
  h1->Draw();
  c1 -> cd(2)->SetLogy();
  h1LaserAmp -> Draw();
  c1 -> cd(3) -> SetLogz();
  h2AmpCharge -> Draw("COLZ");
  h2BlankAmpCharge -> SetLineColor(kRed);
  h2BlankAmpCharge -> Draw("SAME Cont2");
  c1 -> cd(4) -> SetLogz();


  TEped->Draw("main_FDDref.Charge_GATE>>h1Ped(1200,-300,0)", "", "GOFF");
  TH1 *h1Ped = TEped->GetHistogram();
  h1Ped->SetTitle("Charge distribution;Charge (fC);Entries");

  cout << "before scaling" << endl;

  // Calculation of fraction of the defined threshold
  outValues blankResult = CalculateFraction(h1Ped, threshold, kTRUE);
  outValues signalResult = CalculateFraction(h1, threshold, kTRUE);



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
  TCanvas * cped = new TCanvas("cped","Pedestal", w, h*0.6);
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
  grOccupancy -> SetTitle("Ocupancy plot;Threshold (pC);Ocupancy (#lambda) [PE/trigger]");

  TGraph *grFraction = new TGraph();
  grFraction -> SetTitle("Fraction plot;Threshold (pC);Fraction (f)");

  TGraph *grFracOccupancy = new TGraph();
  grFracOccupancy -> SetTitle("Fraction vs Occupancy;Threshold fraction (f);Ocupancy (#lambda) [PE/trigger]");


  // Calculate occupancy using the scaled histograms
  cout << "after scaling" << endl;
  Int_t Index = 0;
  Double_t occupancy;

  for (Double_t i = -290; i < -200; i += 1, Index++){

    outValues afblankResult = CalculateFraction(h1Ped, i, kFALSE);
    outValues afsignalResult = CalculateFraction(h1, i, kFALSE);
    afsignalResult.belowTrs != 0? occupancy = -TMath::Log(double(afsignalResult.belowTrs)/(double(afblankResult.fraction)*double(afsignalResult.totalN))) : 0;

    grOccupancy -> SetPoint(Index, i, occupancy);
    grFraction -> SetPoint(Index, i, double(afblankResult.fraction));
    grFracOccupancy -> SetPoint(Index, double(afblankResult.fraction), occupancy);

    // cout << Index;
    // cout << ") thrs = " << i << "\toccupancy = " << occupancy;
    // cout << "\tf = " << afblankResult.fraction << endl;
  }

  TCanvas * cThreshold = new TCanvas("cThreshold"," Threshold", w, h);
  cThreshold -> Divide(1,3);
  cThreshold -> cd(1);
  grOccupancy -> Draw("AL*");
  cThreshold -> cd(2);
  grFraction -> Draw("AL*");
  cThreshold -> cd(3);
  grFracOccupancy -> Draw("AL*");


  CalculateSPE(h1, h1Ped, -280, kTRUE);

	//TE->StartViewer();

}

outValues CalculateFraction(TH1* histo, Double_t threshold, Bool_t print = kTRUE){
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

SPEValues CalculateSPE(TH1* h1Signal, TH1* h1Blank, Double_t threshold, Bool_t print = kTRUE){

  // Obtaining the mean and stdDev of the distributions
  Double_t BlankMean = h1Blank -> GetMean();
  Double_t BlankMeanError = h1Blank -> GetMeanError();
  Double_t BlankStdDev = h1Blank -> GetStdDev();
  Double_t BlankStdDevError = h1Blank -> GetStdDevError();

  Double_t signalMean = h1Signal -> GetMean();
  Double_t signalMeanError = h1Signal -> GetMeanError();
  Double_t signalStdDev = h1Signal -> GetStdDev();
  Double_t signalStdDevError = h1Signal -> GetStdDevError();

  // Calculation of the SPE properties: occupancy, mean, standard deviation, etc..
  SPEValues result;

  // Single photo-electron mean calculation E[\psi]= (E[T]-E[B])/E[L]
  outValues BlankResult = CalculateFraction(h1Blank, -280, kFALSE);
  outValues SignalResult = CalculateFraction(h1Signal, -280, kFALSE);
  result.occupancy = -TMath::Log(double(SignalResult.belowTrs)/(double(BlankResult.fraction)*double(SignalResult.totalN)));
  result.mean = (h1Signal->GetMean() - h1Blank->GetMean()) / result.occupancy;

  // Single photo-electron Variance calculation V[\psi]= ((V[T]-V[B])/E[L]) - E[psi]
  result.variance = ((signalStdDev - BlankStdDev) / result.occupancy) - result.mean*result.mean;
  result.stdDev = TMath::Sqrt( TMath::Abs(result.variance));

  //Statitical uncertainties Eq. 16 from paper
  result.meanUncertainty = (result.occupancy*(result.mean*result.mean + result.variance) + 2*BlankStdDev)/(double(SignalResult.totalN)*result.occupancy*result.occupancy) + (result.mean*result.mean * (TMath::Exp(result.occupancy) + 1 - (double(BlankResult.fraction)) ))/ (double(BlankResult.fraction)*double(SignalResult.totalN)*result.occupancy*result.occupancy);

  result.stdDevUncertainty = 0;

  if(print){
    cout << "\n> Occupancy = " << result.occupancy;
    cout << "\n> Mean: E[psi] = " << result.mean;
    cout << "\n> Variance: V[psi] = " << result.variance;
    cout << "\n> Sandard Dev.: STDev[psi] = " << result.stdDev;
    cout << "\n> Mean Uncertainties: V[E[psi]] = " << result.meanUncertainty;
    cout << "-> " << result.meanUncertainty*100 << "%\n";
    cout << endl;
  }

  return result;
}
