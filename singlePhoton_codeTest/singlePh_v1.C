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
  Double_t VarianceUncertainty = 0;
};

outValues CalculateFraction(TH1* ,Double_t ,Bool_t );
SPEValues CalculateSPE(TH1*, TH1*, Double_t, Bool_t );
Double_t NormHisto(TH1*, TH1*);

void singlePh_v1(const Char_t* SChannel = "main_FDDref", Double_t AmpWindowMin = 50, Double_t AmpWindowMax = 1010, const Char_t* fileRoot="results.root",const Char_t* filePedestal="pedestal.root"){

  TString aFile = fileRoot;
  TFile *f = TFile::Open(aFile);
  TTree *TELaser = (TTree*)f->Get("RawDataTree");
  // TE->Print();
  // TE->Show(10);

  TString pedFile = filePedestal;
  TFile *fped = TFile::Open(pedFile);
  TTree * TEped = (TTree*)fped->Get("RawDataTree");

  gStyle->SetOptStat(0);

  //------------------------------------------------------
  // Parameters secction
  //------------------------------------------------------
  //Char_t* SChannel = (char*)"main_FDDref";
  //Double_t AmpWindowMin = 50;
  //Double_t AmpWindowMax = 1000;
  Double_t AmpWindowSaturation = 1080;
  Double_t threshold = -283;

  //------------------------------------------------------
  //------------------------------------------------------

  TELaser->Draw(TString::Format("%s.Amplitude>>h1LaserAmp(1200,0,1200)", SChannel), "", "GOFF");
  TH1 *h1LaserAmp = TELaser->GetHistogram();
  h1LaserAmp->SetTitle("Amplitude distribution Laser Induced;Amplitude (mV);Entries");

  TELaser->Draw(TString::Format("%s.Charge_GATE:%s.Amplitude>>h2AmpCharge(1200,0,1200, 1200,-300,0)", SChannel, SChannel), "", "GOFF");
  TH1 *h2AmpCharge = TELaser->GetHistogram();
  h2AmpCharge->SetTitle("Amplitude distribution Laser Induced;Amplitude (mV);Charge (fC);Entries");

  TEped->Draw(TString::Format("%s.Charge_GATE:%s.Amplitude>>h2BlankAmpCharge(1200,0,1200, 1200,-300,0)",SChannel, SChannel), "", "GOFF");
  TH1 *h2BlankAmpCharge = TEped->GetHistogram();
  h2BlankAmpCharge->SetTitle("Amplitude distribution Blank;Amplitude (mV);Charge (fC);Entries");

  TELaser->Draw(TString::Format("%s.Charge_GATE>>h1(1200,-300,0)", SChannel), "", "GOFF");
  //TELaser->Draw("main_FDDref.Charge_GATE>>h1(1200,-300,0)", " main_FDDref.Amplitude < 150 && main_FDDref.Amplitude > 50 ", "GOFF");
  TH1 *h1 = TELaser->GetHistogram();
  h1->SetTitle("Charge distribution;Charge (fC);Entries");

  TELaser->Draw(TString::Format("%s.Charge_GATE>>h1ChargeCut(1200,-300,0)", SChannel), TString::Format(" %s.Amplitude > %f && %s.Amplitude < %f ", SChannel, AmpWindowMin, SChannel, AmpWindowMax), "GOFF");
  TH1 *h1ChargeCut = TELaser->GetHistogram();
  h1ChargeCut->SetTitle("Charge distribution;Charge (fC);Entries");
  h1ChargeCut -> SetLineColor(kMagenta);

  TELaser->Draw(TString::Format("%s.Charge_GATE>>h1ChargeCutOver(1200,-300,0)", SChannel), TString::Format(" %s.Amplitude > %f && %s.Amplitude < %f ", SChannel, AmpWindowMax, SChannel, AmpWindowSaturation), "GOFF");
  TH1 *h1ChargeCutOver = TELaser->GetHistogram();
  h1ChargeCutOver->SetTitle("Charge distribution;Charge (fC);Entries");
  h1ChargeCutOver -> SetLineColor(kMagenta+2);

  TELaser->Draw(TString::Format("%s.Charge_GATE>>h1ChargeCutUnder(1200,-300,0)", SChannel), TString::Format(" %s.Amplitude < %f  ",SChannel ,AmpWindowMin), "GOFF");
  TH1 *h1ChargeCutUnder = TELaser->GetHistogram();
  h1ChargeCutUnder->SetTitle("Charge distribution;Charge (fC);Entries");
  h1ChargeCutUnder -> SetLineColor(kMagenta+4);

  TEped->Draw(TString::Format("%s.Charge_GATE>>h1ChargeCutPedestal(1200,-300,0)", SChannel), TString::Format(" %s.Amplitude > %f && %s.Amplitude < %f ", SChannel, AmpWindowMin, SChannel, AmpWindowMax), "GOFF");
  TH1 *h1ChargeCutPedestal = TEped->GetHistogram();
  h1ChargeCutPedestal->SetTitle("Charge distribution;Charge (fC);Entries");
  h1ChargeCutPedestal -> SetLineColor(kRed);

  h1ChargeCutPedestal->Scale(NormHisto(h1ChargeCut,h1ChargeCutPedestal));

  Double_t w = 1400;
  Double_t h = 1000;
  TCanvas * c1 = new TCanvas("c1","Signal charge", w, h);
  c1 -> Divide(2,2);
  c1 -> cd(1)->SetLogy();
  h1->Draw();
  h1ChargeCut -> Draw("SAME");
  h1ChargeCutOver -> Draw("SAME");
  h1ChargeCutUnder -> Draw("SAME");
  c1 -> cd(2) -> SetLogz();
  h2AmpCharge -> Draw("COLZ");
  h2BlankAmpCharge -> SetLineColor(kRed);
  h2BlankAmpCharge -> Draw("SAME Cont2");
  c1 -> cd(3) -> SetLogz();
  // Draw charge after cutting
  h1ChargeCut -> Draw();
  h1ChargeCutPedestal -> Draw("SAME");
  c1 -> cd(4)->SetLogy();
  h1LaserAmp -> Draw();


  TEped->Draw(TString::Format("%s.Charge_GATE>>h1Ped(1200,-300,0)", SChannel), "", "GOFF");
  TH1 *h1Ped = TEped->GetHistogram();
  h1Ped->SetTitle("Charge distribution;Charge (fC);Entries");

  cout << "before scaling" << endl;

  // Calculation of fraction of the defined threshold
  outValues blankResult = CalculateFraction(h1Ped, threshold, kTRUE);
  outValues signalResult = CalculateFraction(h1, threshold, kTRUE);

  // Getting the maximum bin and factor for normalization
  Int_t binmaxSig = h1->GetMaximumBin();
  Int_t nSig = h1->GetBinContent(binmaxSig);
  Int_t binmax = h1Ped->GetMaximumBin();
  Int_t n = h1Ped->GetBinContent(binmax);

  h1Ped->Scale(NormHisto(h1,h1Ped));

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
  printf("\nRMS = %f StdDev = %f\n",h1 -> GetRMS(), h1 -> GetStdDev());

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

  TGraph *grMean = new TGraph();
  grMean -> SetTitle("Mean ;Threshold (pC); Mean Chargue (pC)");

  TGraph *grVariance = new TGraph();
  grVariance -> SetTitle("Standrd Deviation;Threshold (pC);Variance (pC)");

  TGraph *grMeanUncert = new TGraph();
  grMeanUncert -> SetTitle("Mean Uncertainty;Threshold (pC);Uncertainty");

  TGraph *grVarianceUncert = new TGraph();
  grVarianceUncert -> SetTitle("Standrd Deviation Uncertainty;Threshold (pC);Uncertainty");

  // Calculate occupancy using the scaled histograms
  cout << "after scaling" << endl;
  Int_t Index = 0;
  Double_t occupancy;
  SPEValues SPEResult;

  for (Double_t i = threshold; i < -200; i += 1, Index++){

    outValues afblankResult = CalculateFraction(h1Ped, i, kFALSE);
    outValues afsignalResult = CalculateFraction(h1, i, kFALSE);

    SPEResult = CalculateSPE(h1, h1Ped, i, kFALSE);
    grOccupancy -> SetPoint(Index, i, SPEResult.occupancy);
    grFraction -> SetPoint(Index, i, double(afblankResult.fraction));

    grMean -> SetPoint(Index, i, SPEResult.mean);
    grVariance -> SetPoint(Index, i, SPEResult.variance);

    grMeanUncert -> SetPoint(Index, i, SPEResult.meanUncertainty);
    grVarianceUncert -> SetPoint(Index, i, SPEResult.VarianceUncertainty);

     // cout << Index;
     // cout << ") thrs = " << i << "\toccupancy = " << occupancy << "\t  SPEOcupancy = " << SPEResult.occupancy;
     // cout << "denominador = " << double(afblankResult.fraction)*double(afsignalResult.totalN);
     // cout << "\tf = " << afblankResult.fraction << endl;
  }

  cout << "Calculation of SPE (Full distribution)";
  SPEResult = CalculateSPE(h1, h1Ped, -279, kTRUE);
  cout << "Calculation of SPE (Disstribution with charge cut)";
  SPEResult = CalculateSPE(h1ChargeCut, h1ChargeCutPedestal, -279, kTRUE);

  TCanvas * cThreshold = new TCanvas("cThreshold"," Threshold", w, h);
  cThreshold -> Divide(1,2);
  cThreshold -> cd(1);
  grOccupancy -> Draw("AL*");
  cThreshold -> cd(2);
  grFraction -> Draw("AL*");

  TCanvas * cMeanVar = new TCanvas("cMeanVar","Mean and Variance", w, h);
  cMeanVar -> Divide(1,2);
  cMeanVar -> cd(1);
  grMean -> Draw("AL*");
  cMeanVar -> cd(2);
  grVariance -> Draw("AL*");


  TCanvas * cUncert = new TCanvas("cUncert"," Uncertainty", w, h);
  cUncert -> Divide(1,2);
  cUncert -> cd(1) -> SetLogy();
  grMeanUncert -> Draw("AL*");
  cUncert -> cd(2) -> SetLogy();
  grVarianceUncert -> Draw("AL*");

	//TE->StartViewer();

}


Double_t NormHisto(TH1* h1, TH1* h1Ped){
  Int_t binmaxSig = h1->GetMaximumBin();
  Int_t nSig = h1->GetBinContent(binmaxSig);
  cout << "binmax="<< binmaxSig << endl;
  cout << "xbin="<< nSig << endl;

  Int_t binmax = h1Ped->GetMaximumBin();
  Int_t n = h1Ped->GetBinContent(binmax);
  cout << "binmax="<< binmax << endl;
  cout << "xbin="<< n << endl;

  Double_t Factor = float(nSig)/float(n);
  cout << "NormFactor="<< Factor << endl;
  return Factor;
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
    cout << "\nThreshold bin = " << histo ->  FindBin(threshold);
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
  Double_t BlankVariance = h1Blank -> GetStdDev()*h1Blank -> GetStdDev(); //Variance = StdDev^2
  Double_t BlankVarianceError = h1Blank -> GetStdDevError();

  Double_t signalMean = h1Signal -> GetMean();
  Double_t signalMeanError = h1Signal -> GetMeanError();
  Double_t signalVariance = h1Signal -> GetStdDev()*h1Signal -> GetStdDev();//Variance = StdDev^2
  Double_t signalVarianceError = h1Signal -> GetStdDevError();

  // Calculation of the SPE properties: occupancy, mean, standard deviation, etc..
  SPEValues result;

  // Single photo-electron mean calculation E[\psi]= (E[T]-E[B])/E[L]
  outValues BlankResult = CalculateFraction(h1Blank, threshold, print);
  outValues SignalResult = CalculateFraction(h1Signal, threshold, print);
  result.occupancy = -TMath::Log(double(SignalResult.belowTrs)/(double(BlankResult.fraction)*double(SignalResult.totalN)));

  if (result.occupancy != 0 && SignalResult.totalN != 0){
    result.mean = (h1Signal->GetMean() - h1Blank->GetMean()) / result.occupancy;

    // Single photo-electron Variance calculation V[\psi]= ((V[T]-V[B])/E[L]) - E[psi]^2
    result.variance = ((signalVariance - BlankVariance) / result.occupancy) - result.mean*result.mean;
    result.stdDev = TMath::Sqrt( TMath::Abs(result.variance));

    //Statitical uncertainties Eq. 16 from paper
    result.meanUncertainty = (result.occupancy*(result.mean*result.mean + result.variance) + 2*BlankVariance)/(double(SignalResult.totalN)*result.occupancy*result.occupancy) + (result.mean*result.mean * (TMath::Exp(result.occupancy) + 1 - 2*double(BlankResult.fraction)))/ (double(BlankResult.fraction)*double(SignalResult.totalN)*result.occupancy*result.occupancy);

    //Statitical uncertainties Eq. 16 from paper
    result.VarianceUncertainty = ((result.mean*result.mean + result.variance)*(result.mean*result.mean + result.variance)*(TMath::Exp(result.occupancy) + 1 - 2*double(BlankResult.fraction)) )/(double(BlankResult.fraction)*double(SignalResult.totalN)*result.occupancy*result.occupancy);
  }
  else{
    result.mean  = 0;
    result.variance = 0;
    result.stdDev = 0;
    result.meanUncertainty = 0;
    result.VarianceUncertainty = 0;
  }

  if(print){
    cout << "\n> Occupancy = " << result.occupancy;
    cout << "\n> Mean: E[psi] = " << result.mean;
    cout << "\n> Variance: V[psi] = " << result.variance;
    cout << "\n> Standard Dev.: STDev[psi] = " << result.stdDev;
    cout << "\n> Mean Uncertainties: V[E[psi]] = " << result.meanUncertainty;
    cout << "-> " << result.meanUncertainty*100/result.mean << "%";
    cout << "\n> StdDev Uncertainties: V[V[psi]] = " << result.VarianceUncertainty;
    cout << "-> " << result.VarianceUncertainty*100/result.variance << "%\n";
    cout << endl;
  }

  return result;
}
