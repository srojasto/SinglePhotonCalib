#include "SPECalculation.h"
#include "/home/solangel/sw/digdataprocessing/DigVariables.cpp"

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
    cout << "> Occupancy = " << result.occupancy;
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

TH1* GetQHisto(const Char_t* fileRoot="root/11otacek_2400v_4PMTs_2nd.root", const Char_t* channel = "main_FDDA"){
  TString aFile = fileRoot;
  TFile *f = TFile::Open(aFile);
  TTree *TE = (TTree*)f->Get("RawDataTree");

  MAIN_signal PMT1signal; // structure from "digdataprocessing/DigVariables.cpp"
  TE->SetBranchAddress(channel, &PMT1signal);
  Int_t nEntry=TE->GetEntries();

  Double_t binMax, binMin, binSize; //Variables to store the maximum and minimum range and number of bins for the histogram
  TE->GetEntry(0);
  binMax = PMT1signal.Charge_GATE;
  binMin = PMT1signal.Charge_GATE;

  Double_t binMaxAux = PMT1signal.Charge_GATE;
  Double_t binMinAux = PMT1signal.Charge_GATE;

  for(Int_t i=1; i<nEntry; i++){
    TE->GetEntry(i);
    if(PMT1signal.Charge_GATE == -1) continue;

    // Store current value in the auxiliar variables
    binMaxAux = PMT1signal.Charge_GATE;
    binMinAux = PMT1signal.Charge_GATE;

    // if the current value is grater (lower) than the previous, then stores the value
    if(binMaxAux>binMax) binMax = binMaxAux;
    if(binMaxAux<binMin) binMin = binMinAux;

  }
  cout << "Minimum value: " << binMin;
  cout << "\tMaximum value: " << binMax;
  binSize=4096;
  cout << "\tNumber of bins: " << binSize;
  cout << endl;

  TH1* h1 = new TH1F("h1", "title", binSize+20, binMin-10, binMax+10);
  for(Int_t i=0; i<nEntry; i++){
    TE->GetEntry(i);
    if(PMT1signal.Charge_GATE == -1) continue;
    h1 -> Fill(PMT1signal.Charge_GATE);
  }
  return h1;
}


TH1* GetAmpHisto(const Char_t* fileRoot="root/11otacek_2400v_4PMTs_2nd.root", const Char_t* channel = "main_FDDA"){
  TString aFile = fileRoot;
  TFile *f = TFile::Open(aFile);
  TTree *TE = (TTree*)f->Get("RawDataTree");

  MAIN_signal PMT1signal; // structure from "digdataprocessing/DigVariables.cpp"
  TE->SetBranchAddress(channel, &PMT1signal);
  Int_t nEntry=TE->GetEntries();

  Double_t binMax, binMin, binSize; //Variables to store the maximum and minimum range and number of bins for the histogram
  TE->GetEntry(0);
  binMax = PMT1signal.Amplitude;
  binMin = PMT1signal.Amplitude;

  Double_t binMaxAux = PMT1signal.Amplitude;
  Double_t binMinAux = PMT1signal.Amplitude;

  for(Int_t i=1; i<nEntry; i++){
    TE->GetEntry(i);
    if(PMT1signal.Amplitude == -1) continue;

    // Store current value in the auxiliar variables
    binMaxAux = PMT1signal.Amplitude;
    binMinAux = PMT1signal.Amplitude;

    // if the current value is grater (lower) than the previous, then stores the value
    if(binMaxAux>binMax) binMax = binMaxAux;
    if(binMaxAux<binMin) binMin = binMinAux;

  }
  cout << "Minimum value: " << binMin;
  cout << "\tMaximum value: " << binMax;
  binSize=4096;
  cout << "\tNumber of bins: " << binSize;
  cout << endl;

  TH1* h1 = new TH1F("h1", "title", binSize+20, binMin-10, binMax+10);
  for(Int_t i=0; i<nEntry; i++){
    TE->GetEntry(i);
    if(PMT1signal.Amplitude == -1) continue;
    h1 -> Fill(PMT1signal.Amplitude);
  }
  return h1;
}
