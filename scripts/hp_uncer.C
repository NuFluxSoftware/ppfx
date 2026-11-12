#include <vector>

const int Nhel = 4;
const char* hel[Nhel]     = {"numu","numubar","nue","nuebar"};
const char* fhel[Nhel]    = {"#nu_{#mu}","#bar{#nu}_{#mu}","#nu_{e}","#bar{#nu}_{e}"};

const int Nuniv = 100;

const int Ncorr = 11;  //Initially it was 10

//const char* corr[Ncorr]  = {"pcpi","pck","ncpi","pcnu","pipinc","mes","nua","att","abs","oth","tot"};
const char* fcorr[Ncorr] = {"ThinTargetpCPion","ThinTargetpCKaon","ThinTargetnCPion","ThinTargetpCNucleon", "ThinTargetpipCpip"
                            "ThinTargetMesonIncident","ThinTargetnucleonA","TargetAttenuation","TotalAbsorption","Other","total"};
//const int corr_col[Ncorr] = {kBlue,kBlue,kOrange,kOrange,2,44,44,3,3,9,1};
const int corr_col[Ncorr] = {44,3,9,1,3,4,4,kOrange, kOrange,44,2};

const int corr_sty[Ncorr] = {1,1,1,1,2,1,2,1,2,2,1};


std::vector<TH1D*> vhist[Ncorr];

std::vector<TH1D*> getvhist(TFile* fIn, int icorr, const char* chel);
TH1D* get_cv_with_error_band(std::vector<TH1D*> vhistos);
TH1D* get_fract_error(TH1D* hin);
void plotting(int ihel);
void drawleg(TLegend* leg, TH1D* href, int id);

void hp_uncer(const char* fnameIn){

  std::cout<<"starting() "<<fnameIn<<std::endl;
  TFile* fileIn  = new TFile(fnameIn,"read");
  for(int i=0;i<Nhel;i++){
    std::cout<<"plotting "<<hel[i]<<std::endl;

    for(int j=0;j<Ncorr;j++){
      vhist[j] = getvhist(fileIn,j,hel[i]);
    }
    plotting(i);

    for(int j=0;j<Ncorr;j++)vhist[j].clear();
  }

}
std::vector<TH1D*> getvhist(TFile* fIn, int icorr, const char* chel){

  std::vector<TH1D*> vtmp;

      const char* folder_name_to_exclude = "nom";
  TDirectory* udir = fIn->GetDirectory(Form("%s_%s",chel,corr[icorr]));
  if (strcmp(udir->GetName(), folder_name_to_exclude) == 0) {
        std::cout << "Skipping folder: " << Form("%s_%s", chel, corr[icorr]) << std::endl;
        return vtmp;
    }

  std::cout<<"=> Reading: "<<Form("%s_%s",chel,corr[icorr])<<std::endl;
  TIter next(udir->GetListOfKeys());
  TKey* key;
  while((key= (TKey*)next())){
    TClass* cl = gROOT->GetClass(key->GetClassName());
    if(!cl->InheritsFrom("TH1D"))continue;
    vtmp.push_back((TH1D*)key->ReadObj());
  }
  return vtmp;
}

void plotting(int ihel){

  TGaxis::SetMaxDigits(5);

  TH1D* herr[Ncorr];
  for(int kk=0;kk<Ncorr;kk++)herr[kk] = get_cv_with_error_band(vhist[kk]);

  //fractional error:
  TH1D* hfes[Ncorr];
  std::vector<TH1D*> vec_fe;
  for(int kk=0;kk<Ncorr;kk++){
    hfes[kk] = get_fract_error(herr[kk]);
    vec_fe.push_back(hfes[kk]);
  }

  gROOT->SetStyle("Plain");
  gStyle->SetLegendBorderSize(0);

double x1leg = 0.16;
double deltax = 0.18;
double sizex = 0.30;
double y1leg = 0.70;
double y1aleg = 0.70;
double y2leg = 0.80;

TCanvas* cfe = new TCanvas();
TLatex* lat = new TLatex();
lat->SetNDC();

// Create Legends
TLegend* leg1 = new TLegend(x1leg, y1leg, x1leg + sizex, y2leg);

x1leg += deltax;
TLegend* leg2 = new TLegend(x1leg, y1leg, x1leg + sizex, y2leg);

x1leg += deltax;
TLegend* leg3 = new TLegend(x1leg, y1leg, x1leg + sizex, y2leg);

x1leg += deltax;
TLegend* aleg = new TLegend(x1leg, y1aleg, x1leg + sizex, y2leg);

// Set fill color (if needed)
leg1->SetFillColor(0);
leg2->SetFillColor(0);
leg3->SetFillColor(0);
aleg->SetFillColor(0);

// Set line properties for histograms
for (int kk = 0; kk < Ncorr; kk++) {
    hfes[kk]->SetLineColor(corr_col[kk]);
    hfes[kk]->SetLineStyle(corr_sty[kk]);
    hfes[kk]->SetLineWidth(2);
}

// draw histograms
for (int kk = 0; kk < Ncorr; kk++) {
    if (kk == 0) {
        //hfes[kk]->getxaxis()->setrangeuser(0, 12 - 0.0001);
        hfes[kk]->GetXaxis()->SetRangeUser(0, 6 - 0.0001);
        // set the minimum y-axis range to a small positive value
        //const double minvisiblevalue = 0.01; // adjust this value as needed
        //const double maxvisiblevalue = 0.25; // adjust this value as needed

        //hfes[kk]->GetYaxis()->SetRangeUser(minvisiblevalue, maxvisiblevalue);
        hfes[kk]->GetYaxis()->SetRangeUser(0, 0.45 - 0.0001);

        hfes[kk]->GetYaxis()->SetRangeUser(0.001, 0.20);
        hfes[kk]->GetXaxis()->SetTitle("Neutrino energy (GeV)");
        hfes[kk]->GetYaxis()->SetTitle("Fractional Uncertainties");
        hfes[kk]->Draw("hist");
    } else {
        hfes[kk]->Draw("histsame");
    }
}

// Draw legends
drawleg(leg1, hfes[4], 4); drawleg(leg1, hfes[0], 0); drawleg(leg1, hfes[1], 1);
drawleg(leg2, hfes[6], 6); drawleg(leg2, hfes[2], 2); drawleg(leg2, hfes[3], 3);
drawleg(leg3, hfes[7], 7); drawleg(leg3, hfes[5], 5); drawleg(leg3, hfes[8], 8);
drawleg(aleg, hfes[9], 9); drawleg(aleg, hfes[10], 10);



  double size_leg = 0.020;
  leg1->SetTextSize(size_leg);
  leg2->SetTextSize(size_leg);

  leg3->SetTextSize(size_leg);
  aleg->SetTextSize(size_leg);

  leg1->Draw();
  leg2->Draw();
  leg3->Draw();
  aleg->Draw();

  lat->DrawLatex(0.2,0.86,"Hadron Production Uncertainties");
  lat->DrawLatex(0.78,0.85,Form("#scale[1.2]{%s}",fhel[ihel]));

  std::cout<<"=> "<<Form("tot_fe_%s.png",hel[ihel])<<std::endl;
  cfe->SaveAs(Form("tot_fe_%s.png",hel[ihel]));

}

TH1D* get_cv_with_error_band(std::vector<TH1D*> vhistos){
  int Nuniv = vhistos.size();
  if(Nuniv==0)exit(1);
  TH1D* hband = (TH1D*)vhistos[0]->Clone();

    int Nbins   = hband->GetXaxis()->GetNbins();

    for(int jj=1;jj<=Nbins;jj++){

      //mean:
      double mean = 0;
      for(int ii=0;ii<Nuniv;ii++){
        mean += vhistos[ii]->GetBinContent(jj);
      }
      mean /= double(Nuniv);

      //error:
      double err = 0.0;
      for(int ii=0;ii<Nuniv;ii++){
        err += pow(vhistos[ii]->GetBinContent(jj) - mean,2);
      }
      err /= double(Nuniv);
      err = sqrt(err);
      hband->SetBinContent(jj,mean);
      hband->SetBinError(jj,err);
    }

    return hband;

  }

TH1D* get_fract_error(TH1D* hin){

    TH1D* hfe = (TH1D*)hin->Clone();

    int nbins   = hfe->GetXaxis()->GetNbins();
    for(int jj=1;jj<=nbins;jj++){

      double cont = hin->GetBinContent(jj);
      double err = hin->GetBinError(jj);
      if(cont>0)hfe->SetBinContent(jj,err/cont);

    }

    return hfe;
  }

void drawleg(TLegend* leg, TH1D* href, int id){
  leg->AddEntry(href,Form("%s",cleg[id]),"l");
}

