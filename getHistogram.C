//#include "/home/solangel/sw/digdataprocessing/DigVariables.cpp"
#include "SPECalculation.C"

void getHistogram(){
  printf("Program \n");

  const Char_t* fileRoot="root/11otacek_2300v_4PMTs_3rd.root";

  TString aFile = fileRoot;
  TFile *f = TFile::Open(aFile);
  TTree *TE = (TTree*)f->Get("RawDataTree");

  MAIN_signal PMT1signal, PMT2signal;
  TE->SetBranchAddress("main_FDDA", &PMT1signal);

  Int_t nEntry=TE->GetEntries();

  TH1* h1 = new TH1F("h1", "title", 400,-10,400);

  Double_t mean = 0;
  Double_t variance = 0;
  Double_t sum = 0;
  Int_t ValidEntries = 0;


  //Calculate Charge Mean
  for(Int_t i=0; i<nEntry; i++){
    TE->GetEntry(i);
    if(PMT1signal.Charge_GATE == -1) continue;
    sum += PMT1signal.Charge_GATE;
    ValidEntries++;
    h1 -> Fill(PMT1signal.Charge_GATE);
  }
  mean = sum/ValidEntries;
  printf("\nMean(Q) = Sum/validEntires = %f/%i = %f\n",sum, ValidEntries, mean);

  //Calculate Charge  variance
  sum = 0;
  ValidEntries = 0;
  for(Int_t i=0; i<nEntry; i++){
    TE->GetEntry(i);
    if(PMT1signal.Charge_GATE == -1) continue;
    sum += (PMT1signal.Charge_GATE - mean)*(PMT1signal.Charge_GATE - mean) ;
    ValidEntries++;
  }
  variance = sum/(double)ValidEntries;
  printf("Variance(Q) = Sum/validEntires = %f/%i = %f\n",sum, ValidEntries, variance);

  printf("From histogram\nMean(Q) = %f\tVarience(Q) = %f\n", h1->GetMean(), h1->GetStdDev()*h1->GetStdDev());


  TCanvas * c1 = new TCanvas("c1","Signal charge", 800,800);
  c1 -> SetLogy();
  h1 -> Draw();


  //Calculate Amplitude Mean
  mean = 0;
  variance = 0;
  sum = 0;
  ValidEntries = 0;

  for(Int_t i=0; i<nEntry; i++){
    TE->GetEntry(i);
    if(PMT1signal.Amplitude == -1) continue;
    sum += PMT1signal.Amplitude;
    ValidEntries++;
    h1 -> Fill(PMT1signal.Amplitude);
  }
  mean = sum/ValidEntries;
  printf("\nMean(Amp) = Sum/validEntires = %f/%i = %f\n",sum, ValidEntries, mean);

  //Calculate Amplitude variance
  sum = 0;
  ValidEntries = 0;
  for(Int_t i=0; i<nEntry; i++){
    TE->GetEntry(i);
    if(PMT1signal.Amplitude == -1) continue;
    sum += (PMT1signal.Amplitude - mean)*(PMT1signal.Amplitude - mean) ;
    ValidEntries++;
  }
  variance = sum/(double)ValidEntries;
  printf("Variance(Amp) = Sum/validEntires = %f/%i = %f\n",sum, ValidEntries, variance);

  TH1* hq = GetQHisto(fileRoot, "main_FDDA");
  hq->SetLineColor(kRed);
  hq->Draw("SAME");

  TCanvas * c2 = new TCanvas("c2","Signal Amplitude", 800,800);
  c2 -> SetLogy();

  TH1* hAmp= GetAmpHisto(fileRoot, "main_FDDA");
  hAmp->SetLineColor(kRed);
  hAmp->Draw("SAME");
  printf("From histogram\nMean = %f\tVarience = %f\n", hAmp->GetMean(), hAmp->GetStdDev()*hAmp->GetStdDev());
}
