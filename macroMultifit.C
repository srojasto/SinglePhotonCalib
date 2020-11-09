// -*- C++ -*-
// $Id$
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sys/stat.h>

bool IsPathExist(const std::string &s)
{
  struct stat buffer;
  return (stat (s.c_str(), &buffer) == 0);
}

void macroMultifit(const Char_t* inputfile="inputSN.txt"){

  string SNumber;
  Double_t w = 1400;
  Double_t h = 1000;
  ifstream myfile (inputfile);
	gStyle->SetOptStat("nerm");
  gStyle->SetOptFit(1111);

  if (myfile.is_open())
  {

    while ( getline (myfile,SNumber) )
    {
      TString path2("results/"+SNumber);
      TString path1("results/"+SNumber);

      if (!IsPathExist("results/"+SNumber)){
        cout << " SN does not exist, creating result folder" << endl;
      std::string res = std::string("results/") + SNumber;
      mkdir(res.c_str(), 0777);
    }


      ofstream outfile ("results/"+SNumber+"/"+SNumber+".txt");
      outfile << "HV/D:FDDAmp/D:FDDAmpError/D:ADrefAmp/D:ADrefAmpError/D:FDDQ/D:FDDQError/D:ADrefQ/D:ADrefQError/D";
      outfile << ":FDDAmp2/D:FDDAmpError2/D:ADrefAmp2/D:ADrefAmpError2/D:FDDQ2/D:FDDQError2/D:ADrefQ2/D:ADrefQError2/D"<< endl;
      if(outfile.is_open()){

        TCanvas * c1 = new TCanvas("c1","ADref PMT Analysis", w, h);
        c1->Divide(2,3);
        TCanvas * c2 = new TCanvas("c2","FDDref PMT Analysis", w, h);
        c2->Divide(2,3);
        TCanvas * c3 = new TCanvas("c3","ADref PMT Analysis - pedestal", w, h);
        c3->Divide(2,3);
        TCanvas * c4 = new TCanvas("c4","FDDref PMT Analysis - pedestal", w, h);
        c4->Divide(2,3);



      for (int i = 180; i < 241; i++){
        std::string voltage = std::to_string(i);


        if(!IsPathExist(SNumber+"_"+voltage+"0v"+"/result.root")){
          cout << "voltage does not exist" << endl;

        }
        else{

        string fileRoot=SNumber+"_"+voltage+"0v"+"/result.root";
      	TString aFile = fileRoot;

      	TFile *f = TFile::Open(aFile);
      	TTree *TE = (TTree*)f->Get("RawDataTree");
      	TE->Print();
      	TE->Show(10);



          TE->Draw("main_FDDref.Amplitude>>hAmp(600,0,600)", "", "GOFF");
        TH1 *hAmp = TE->GetHistogram();
        TF1 *fFDDAmp = new TF1("fFDDAmp", "landau", hAmp->GetMean()-hAmp->GetStdDev(), hAmp->GetMean()+hAmp->GetStdDev());
        hAmp->SetTitle("Amplitude;Volts;Entries");


        	TE->Draw("(main_FDDref.time_end-main_FDDref.time_begin) >> hdT(100,0,50)", "", "GOFF");
      	TH1 *hdT = TE->GetHistogram();
      	hdT -> SetTitle("Signal width at 10% of peak;time (ns);Entries");

          TE->Draw("main_FDDref.Charge_GATE >> hQ(800,-100,100)", "", "GOFF");
        TH1 *hQ = TE->GetHistogram();
         Double_t par[9];
    TF1 *g1    = new TF1("g1","gaus",85,95);
    TF1 *g2    = new TF1("g2","gaus",98,108);
    TF1 *g3    = new TF1("g3","gaus",110,121);
    TF1 *fFDDQ = new TF1("fFDDQ","gaus(0)+gaus(3)+gaus(6)",85,125);
    fFDDQ->SetLineColor(2);
    hQ->Fit(g1,"R");
    hQ->Fit(g2,"R+");
    hQ->Fit(g3,"R+");
    g1->GetParameters(&par[0]);
    g2->GetParameters(&par[3]);
    g3->GetParameters(&par[6]);
    fFDDQ->SetParameters(par);
    hQ->Fit(fFDDQ,"R+");
        hQ->SetTitle("Signal FDD Charge;Charge;Entries");

        	TE->Draw("(main_FDDref.Charge_GATE):main_FDDref.Amplitude >> h2AmpQ(600,0,600, 200,-100,100)", "", "GOFF");
      	TH1 *h2AmpQ = TE->GetHistogram();
        h2AmpQ->SetTitle("Charge and Amplitude;Amplitude (mV);Charge");

        	TE->Draw("(main_FDDref.time_end-main_FDDref.time_begin):main_FDDref.Amplitude >> h2AmpdT(600,0,600, 100,0,50)", "", "GOFF");
      	TH1 *h2AmpdT = TE->GetHistogram();
        h2AmpdT->SetTitle("Width (10%) and Amplitude;Amplitude (volts);#Deltat (ns)");

        	TE->Draw("(main_FDDref.time_end-main_FDDref.time_begin):main_FDDref.Charge_GATE >> h2QdT(200,-100,100, 100,0,50)", "", "GOFF");
      	TH1 *h2QdT = TE->GetHistogram();
        h2QdT->SetTitle("Width (10%) and Charge;Charge;#Deltat (ns)");


      	c2->cd(1);
      	hAmp->Draw();
        hAmp->Fit(fFDDAmp,"R");
        fFDDAmp->SetRange(fFDDAmp->GetParameter(1)-fFDDAmp->GetParameter(2),fFDDAmp->GetParameter(1)+fFDDAmp->GetParameter(2));
        hAmp->Fit(fFDDAmp,"R");
        fFDDAmp->Draw("L SAME");
      	c2->cd(2);
      	hdT->Draw();
      	c2->cd(3);
      	c2->cd(3)->SetLogz();
      	hQ->Draw("COLZ");
        hQ->Fit(fFDDQ,"R");
        fFDDQ->SetRange(fFDDQ->GetParameter(1)-fFDDQ->GetParameter(2),fFDDQ->GetParameter(1)+fFDDQ->GetParameter(2));
        hQ->Fit(fFDDQ,"R");
        fFDDQ->Draw("L SAME");
      	c2->cd(4);
      	h2AmpQ->Draw("COLZ");
      	c2->cd(5);
      	c2->cd(5)->SetLogz();
      	h2AmpdT->Draw("COLZ");
      	c2->cd(6);
      	h2QdT->Draw("COLZ");

  c2->SaveAs(path2+"/"+SNumber+"_FDDref.pdf","Title:FDDref");


/////////////////////////////


        TE->Draw("main_ADref.Amplitude>>hRefAmp(600,0,600)", "", "GOFF");
      TH1 *hRefAmp = TE->GetHistogram();
      TF1 *fRefAmp = new TF1("fRefAmp", "landau", hRefAmp->GetMean()-hRefAmp->GetStdDev(), hRefAmp->GetMean()+hRefAmp->GetStdDev());
      hRefAmp->SetTitle("Amplitude;Volts;Entries");

        TE->Draw("(main_ADref.time_end-main_ADref.time_begin) >> hRefdT(100,0,50)", "", "GOFF");
      TH1 *hRefdT = TE->GetHistogram();
      hRefdT -> SetTitle("Signal width at 10% of peak;time (ns);Entries");

        TE->Draw("main_ADref.Charge_GATE >> hRefQ(800,-100,100)", "", "GOFF");
      TH1 *hRefQ = TE->GetHistogram();
      TF1 *fRefQ = new TF1("fRefQ", "landau", hRefQ->GetMean()-hRefQ->GetStdDev(), hRefQ->GetMean()+hRefQ->GetStdDev());
      hRefQ->SetTitle("Signal Ref Charge;Charge;Entries");

        TE->Draw("(main_ADref.Charge_GATE):main_ADref.Amplitude >> h2RefAmpQ(600,0,600, 200,-100,100)", "", "GOFF");
      TH1 *h2RefAmpQ = TE->GetHistogram();
      h2RefAmpQ->SetTitle("Charge and Amplitude;Amplitude (mV);Charge");

        TE->Draw("(main_ADref.time_end-main_ADref.time_begin):main_ADref.Amplitude >> h2RefAmpdT(600,0,600, 100,0,50)", "", "GOFF");
      TH1 *h2RefAmpdT = TE->GetHistogram();
      h2RefAmpdT->SetTitle("Width (10%) and Amplitude;Amplitude (volts);#Deltat (ns)");

        TE->Draw("(main_ADref.time_end-main_ADref.time_begin):main_ADref.Charge_GATE >> h2RefQdT(200,-100,100, 100,0,50)", "", "GOFF");
      TH1 *h2RefQdT = TE->GetHistogram();
      h2RefQdT->SetTitle("Width (10%) and Charge;Charge;#Deltat (ns)");


      c1->cd(1);
      hRefAmp->Draw();
      hRefAmp->Fit(fRefAmp,"R");
      fRefAmp->SetRange(fRefAmp->GetParameter(1)-fRefAmp->GetParameter(2),fRefAmp->GetParameter(1)+fRefAmp->GetParameter(2));
      hRefAmp->Fit(fRefAmp,"R");
      fRefAmp->Draw("L SAME");
      c1->cd(2);
      hRefdT->Draw();
      c1->cd(3);
      c1->cd(3)->SetLogz();
      hRefQ->Draw("COLZ");
      hRefQ->Fit(fRefQ,"R");
      fRefQ->SetRange(fRefQ->GetParameter(1)-fRefQ->GetParameter(2),fRefQ->GetParameter(1)+fRefQ->GetParameter(2));
      hRefQ->Fit(fRefQ,"R");
      fRefQ->Draw("L SAME");
      c1->cd(4);
      h2RefAmpQ->Draw("COLZ");
      c1->cd(5);
      c1->cd(5)->SetLogz();
      h2RefAmpdT->Draw("COLZ");
      c1->cd(6);
      h2RefQdT->Draw("COLZ");
  c1->SaveAs(path1+"/"+SNumber+"_ADref.pdf","Title:ADref");

//**************************************************


string fileRoot2=SNumber+"_"+voltage+"0v"+"/pedestal.root";
TString aFile2 = fileRoot2;

TFile *f2 = TFile::Open(aFile2);
TTree *TE2 = (TTree*)f2->Get("RawDataTree");
TE2->Print();
TE2->Show(10);

  TE2->Draw("main_FDDref.Amplitude>>hAmp(600,0,600)", "", "GOFF");
TH1 *hAmp2 = TE2->GetHistogram();
TF1 *fFDDAmp2 = new TF1("fFDDAmp2", "gaus", hAmp2->GetMean()-hAmp2->GetStdDev(), hAmp2->GetMean()+hAmp2->GetStdDev());
hAmp2->SetTitle("Amplitude;Volts;Entries");

  TE2->Draw("(main_FDDref.time_end-main_FDDref.time_begin) >> hdT(100,0,50)", "", "GOFF");
TH1 *hdT2 = TE2->GetHistogram();
hdT2 -> SetTitle("Signal width at 10% of peak;time (ns);Entries");

  TE2->Draw("main_FDDref.Charge_GATE >> hQ(800,-100,100)", "", "GOFF");
TH1 *hQ2 = TE2->GetHistogram();
TF1 *fFDDQ2 = new TF1("fFDDQ2", "gaus", hQ2->GetMean()-hQ2->GetStdDev(), hQ2->GetMean()+hQ2->GetStdDev());
hQ2->SetTitle("Signal FDD Charge;Charge;Entries");

  TE2->Draw("(main_FDDref.Charge_GATE):main_FDDref.Amplitude >> h2AmpQ(600,0,600, 200,-100,100)", "", "GOFF");
TH1 *h2AmpQ2 = TE2->GetHistogram();
h2AmpQ2->SetTitle("Charge and Amplitude;Amplitude (mV);Charge");

  TE2->Draw("(main_FDDref.time_end-main_FDDref.time_begin):main_FDDref.Amplitude >> h2AmpdT(600,0,600, 100,0,50)", "", "GOFF");
TH1 *h2AmpdT2 = TE2->GetHistogram();
h2AmpdT2->SetTitle("Width (10%) and Amplitude;Amplitude (volts);#Deltat (ns)");

  TE2->Draw("(main_FDDref.time_end-main_FDDref.time_begin):main_FDDref.Charge_GATE >> h2QdT(200,-100,100, 100,0,50)", "", "GOFF");
TH1 *h2QdT2 = TE2->GetHistogram();
h2QdT2->SetTitle("Width (10%) and Charge;Charge;#Deltat (ns)");


c4->cd(1);
hAmp2->Draw();
hAmp2->Fit(fFDDAmp2,"R");
fFDDAmp2->SetRange(fFDDAmp2->GetParameter(1)-fFDDAmp2->GetParameter(2),fFDDAmp2->GetParameter(1)+fFDDAmp2->GetParameter(2));
hAmp2->Fit(fFDDAmp2,"R");
fFDDAmp2->Draw("L SAME");
c4->cd(2);
hdT2->Draw();
c4->cd(3);
c4->cd(3)->SetLogz();
hQ2->Draw("COLZ");
hQ2->Fit(fFDDQ2,"R");
fFDDQ2->SetRange(fFDDQ2->GetParameter(1)-fFDDQ2->GetParameter(2),fFDDQ2->GetParameter(1)+fFDDQ2->GetParameter(2));
hQ2->Fit(fFDDQ2,"R");
fFDDQ2->Draw("L SAME");
c4->cd(4);
h2AmpQ2->Draw("COLZ");
c4->cd(5);
c4->cd(5)->SetLogz();
h2AmpdT2->Draw("COLZ");
c4->cd(6);
h2QdT2->Draw("COLZ");

c4->SaveAs(path2+"/"+SNumber+"_FDDref_ped.pdf","Title:FDDref");


/////////////////////////////


TE2->Draw("main_ADref.Amplitude>>hRefAmp(600,0,600)", "", "GOFF");
TH1 *hRefAmp2 = TE2->GetHistogram();
TF1 *fRefAmp2 = new TF1("fRefAmp2", "gaus", hRefAmp2->GetMean()-hRefAmp2->GetStdDev(), hRefAmp2->GetMean()+hRefAmp2->GetStdDev());
hRefAmp2->SetTitle("Amplitude;Volts;Entries");

TE2->Draw("(main_ADref.time_end-main_ADref.time_begin) >> hRefdT(100,0,50)", "", "GOFF");
TH1 *hRefdT2 = TE2->GetHistogram();
hRefdT2 -> SetTitle("Signal width at 10% of peak;time (ns);Entries");

TE2->Draw("main_ADref.Charge_GATE >> hRefQ(800,-100,100)", "", "GOFF");
TH1 *hRefQ2 = TE2->GetHistogram();
TF1 *fRefQ2 = new TF1("fRefQ2", "gaus", hRefQ2->GetMean()-hRefQ2->GetStdDev(), hRefQ2->GetMean()+hRefQ2->GetStdDev());
hRefQ2->SetTitle("Signal Ref Charge;Charge;Entries");

TE2->Draw("(main_ADref.Charge_GATE):main_ADref.Amplitude >> h2RefAmpQ(600,0,600, 200,-100,100)", "", "GOFF");
TH1 *h2RefAmpQ2 = TE2->GetHistogram();
h2RefAmpQ2->SetTitle("Charge and Amplitude;Amplitude (mV);Charge");

TE2->Draw("(main_ADref.time_end-main_ADref.time_begin):main_ADref.Amplitude >> h2RefAmpdT(600,0,600, 100,0,50)", "", "GOFF");
TH1 *h2RefAmpdT2 = TE2->GetHistogram();
h2RefAmpdT2->SetTitle("Width (10%) and Amplitude;Amplitude (volts);#Deltat (ns)");

TE2->Draw("(main_ADref.time_end-main_ADref.time_begin):main_ADref.Charge_GATE >> h2RefQdT(200,-100,100, 100,0,50)", "", "GOFF");
TH1 *h2RefQdT2 = TE2->GetHistogram();
h2RefQdT2->SetTitle("Width (10%) and Charge;Charge;#Deltat (ns)");


c3->cd(1);
hRefAmp2->Draw();
hRefAmp2->Fit(fRefAmp2,"R");
fRefAmp2->SetRange(fRefAmp2->GetParameter(1)-fRefAmp2->GetParameter(2),fRefAmp2->GetParameter(1)+fRefAmp2->GetParameter(2));
hRefAmp2->Fit(fRefAmp2,"R");
fRefAmp2->Draw("L SAME");
c3->cd(2);
hRefdT2->Draw();
c3->cd(3);
c3->cd(3)->SetLogz();
hRefQ2->Draw("COLZ");
hRefQ2->Fit(fRefQ2,"R");
fRefQ2->SetRange(fRefQ2->GetParameter(1)-fRefQ2->GetParameter(2),fRefQ2->GetParameter(1)+fRefQ2->GetParameter(2));
hRefQ2->Fit(fRefQ2,"R");
fRefQ2->Draw("L SAME");
c3->cd(4);
h2RefAmpQ2->Draw("COLZ");
c3->cd(5);
c3->cd(5)->SetLogz();
h2RefAmpdT2->Draw("COLZ");
c3->cd(6);
h2RefQdT2->Draw("COLZ");
c3->SaveAs(path1+"/"+SNumber+"_ADref_ped.pdf","Title:ADref");


        outfile << voltage << "0" << "\t" << fFDDAmp->GetParameter(1) << "\t" << fFDDAmp->GetParError(1);
        outfile << "\t" << fRefAmp->GetParameter(1) << "\t" << fRefAmp->GetParError(1);
        outfile << "\t" << fFDDQ->GetParameter(1) << "\t" << fFDDQ->GetParError(1);
        outfile << "\t" << fRefQ->GetParameter(1) << "\t" << fRefQ->GetParError(1);
        outfile << "\t" << fFDDAmp2->GetParameter(1) << "\t" << fFDDAmp2->GetParError(1);
        outfile << "\t" << fRefAmp2->GetParameter(1) << "\t" << fRefAmp2->GetParError(1);
        outfile << "\t" << fFDDQ2->GetParameter(1) << "\t" << fFDDQ2->GetParError(1);
        outfile << "\t" << fRefQ2->GetParameter(1) << "\t" << fRefQ2->GetParError(1) << std::endl;




          }
        }

        outfile.close();
      }

    }
    myfile.close();

  }
  else cout << "Unable to open file";





}
