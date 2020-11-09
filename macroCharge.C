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

void macroCharge(const Char_t* inputfile="inputSN.txt"){

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


        TCanvas * c2 = new TCanvas("c2","FDDref PMT Analysis", w, h);





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


        	TE->Draw("(main_FDDref.time_end-main_FDDref.time_begin) >> hdT(100,0,50)", "", "GOFF");
      	TH1 *hdT = TE->GetHistogram();
      	hdT -> SetTitle("Signal width at 10% of peak;time (ns);Entries");

          TE->Draw("main_FDDref.Charge_GATE >> hQ(1500,0,270)", "", "GOFF");
        TH1 *hQ = TE->GetHistogram();
        hQ->SetTitle("Signal FDD Charge;Charge;Entries");

         Double_t par[11];

          TF1 *g1 = new TF1("g1","gaus(0)+expo(3)",40,55);
          g1->SetParameters(314,44,0.8,5,-0.01);
          g1->SetParLimits(4,-11,-0.001);
          g1->SetLineColor(kGreen);


          TF1 *g2 = new TF1("g2","gaus",73,118);
            g2->SetLineColor(kBlack);
          g2->SetParameters(600,100,70);

          TF1 *g3 = new TF1("g3","gaus",120,240);
            g3->SetLineColor(kBlue);
          g3->SetParameters(700,140,90);


          TF1 *fFDDQ = new TF1("fFDDQ","g1(0)+gaus(5)+gaus(8)",40,170);
          fFDDQ->SetLineColor(kMagenta);




      	c2->SetLogz();
      	hQ->Draw("COLZ");
        hQ->Fit(g1,"R");
        hQ->Fit(g2,"R+");
        hQ->Fit(g3,"R+");
        g1->GetParameters(&par[0]);
        g2->GetParameters(&par[5]);
        g3->GetParameters(&par[8]);
        fFDDQ->SetParameters(par);
        hQ->Fit(fFDDQ,"R+");
        fFDDQ->Draw("L SAME");


          }
        }

        outfile.close();
      }

    }
    myfile.close();

  }
  else cout << "Unable to open file";





}
