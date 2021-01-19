#include "/home/solangel/sw/digdataprocessing/DigVariables.cpp"

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

  for(Int_t i=0; i<nEntry; i++){
    TE->GetEntry(i);
    if(PMT1signal.Charge_GATE == -1) continue;
    sum += PMT1signal.Charge_GATE;
    ValidEntries++;
    h1 -> Fill(PMT1signal.Charge_GATE);
    //printf("Amplitude %e\r", Amplitude);
  }
  mean = sum/ValidEntries;
  printf("\nMean = Sum/validEntires = %f/%i = %f\n",sum, ValidEntries, mean);

  sum = 0;
  ValidEntries = 0;
  for(Int_t i=0; i<nEntry; i++){
    TE->GetEntry(i);
    if(PMT1signal.Charge_GATE == -1) continue;
    sum += (PMT1signal.Charge_GATE - mean)*(PMT1signal.Charge_GATE - mean) ;
    ValidEntries++;
  }
  variance = sum/(double)ValidEntries;
  printf("Variance = Sum/validEntires = %f/%i = %f\n",sum, ValidEntries, variance);

  printf("From histogram\nMean = %f\tVarience = %f\n", h1->GetMean(), h1->GetStdDev()*h1->GetStdDev());

  TCanvas * c1 = new TCanvas("c1","Signal charge", 800,800);
  c1 -> SetLogy();
  h1 -> Draw();


}
