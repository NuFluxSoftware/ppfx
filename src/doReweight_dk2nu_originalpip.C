
#include "dk2nu/tree/dkmeta.h"
#include "dk2nu/tree/dk2nu.h"
#include "MakeReweight.h"
#include "NuWeight.h"

#include "TFile.h"
#include "TH2D.h"
#include "TChain.h"
#include <fstream>
#include <string.h>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <iomanip>

#include "TH1D.h"

//Some constants::

const int NbinsE   = 1200;        //Change this as desired,IMPORTANT  
const double emin  =   0.;
const double emax  = 120.;
const int Nnuhel   = 4;
const char* nuhel[Nnuhel] = {"numu","numubar","nue","nuebar"};
const char* nulabel[Nnuhel] = {"#nu_{#mu}","#bar#nu_{#mu}","#nu_{e}","#bar#nu_{e}"};
const double pi    = 3.14159265358979323846;

using namespace NeutrinoFluxReweight;

class nu_g4numi;
class Numi2Pdg;

int idx_hel(int pdgdcode);
std::string change_to_xrootd_path(std::string temp);

/*!
 * Run the reweighting for a single file (inputFile) for
 * particular MIPP covariance matrix given in input.xml file. 
 */
void doReweight_dk2nu(const char* inputFile, const char* outputFile, const char* optionsFile, const char* cxxdet, const char* cyydet, const char* czzdet){ 
  
  //const char* thisDir = getenv("PPFX_DIR"); // unused
  int idet = -1;    //this is identity of the detector
  bool doing_precalculated_pos = false;  //Initially we are assuming that there is no precalculated position of the detector
  
  NeutrinoFluxAuxiliar::NuWeight* nuweight;		//namespace::class* pointer, the poineter points to an object of the type in class	  
  if( std::string(cyydet) == "none" && std::string(czzdet) == "none" ){
    idet = atoi(cxxdet);
    doing_precalculated_pos = true;
  }
  else{
    std::vector<double> vdet;
    double xxdet = atof(cxxdet); vdet.push_back(xxdet);
    double yydet = atof(cyydet); vdet.push_back(yydet);
    double zzdet = atof(czzdet); vdet.push_back(zzdet);
    nuweight = new NeutrinoFluxAuxiliar::NuWeight(vdet);
  }
  
  std::cout<< "Instance of MakeReweight()" <<std::endl;
  MakeReweight* makerew = MakeReweight::getInstance();
  makerew->SetOptions(optionsFile); 
  std::cout<< "Theoptions file setup are" <<makerew<<std::endl;   

  std::cout<<"Making an output file to store histograms"<<std::endl;
  TFile* fOut = new TFile(outputFile,"recreate");
  std::cout<<"File name: "<<fOut->GetName()<<std::endl;
  
  int Nuniverses = makerew->GetNumberOfUniversesUsed();
  std::cout<<"The number of universe we are using are:"<< Nuniverses<<std::endl; 
  //printf("\n Number of universes: %d \n", Nuniverses);

//Adding some more from dune script, to get all the values separtely
  TH1D* hnom[Nnuhel];                                                           //1
  TH1D* hcv[Nnuhel];                                                            //2
  TH1D* hthin[Nnuhel][Nuniverses];                                              //3
  TH1D* hpipCpipinc[Nnuhel][Nuniverses];                                       //3
  TH1D* hmesinc[Nnuhel][Nuniverses];                                            //4
  //TH1D* hmipp[Nnuhel][Nuniverses];                                              
  TH1D* hatt[Nnuhel][Nuniverses];                                               //5 
  TH1D* hothers[Nnuhel][Nuniverses];                                            //6
  TH1D* htotal[Nnuhel][Nuniverses];                                             //7
  TH1D* htotalabs[Nnuhel][Nuniverses];                                          //8
  TH1D* httpcpi[Nnuhel][Nuniverses];                                            //9
  TH1D* httpck[Nnuhel][Nuniverses];                                             //10
  TH1D* httncpi[Nnuhel][Nuniverses];                                            //11 
  TH1D* httpcnu[Nnuhel][Nuniverses];                                            //12 
  TH1D* httnua[Nnuhel][Nuniverses];                                            //13 

  for(int ii=0;ii<Nnuhel;ii++){
   const TString xtitle = "E_{#nu} [GeV]";
   const TString ytitle = "##nu_{unoscillated} [m^{-2}]";
    hnom[ii] = new TH1D(Form("hnom_%s",nuhel[ii]), Form("Uncorrected %s flux",nulabel[ii]), NbinsE,emin,emax);
    hnom[ii]->SetXTitle(xtitle);
    hnom[ii]->SetYTitle(ytitle);

    hcv[ii]  = new TH1D(Form("hcv_%s",nuhel[ii]), Form("Fully PPFX-corrected %s flux (central value)", nulabel[ii]),NbinsE,emin,emax);
    hcv[ii]->SetXTitle(xtitle);
    hcv[ii]->SetYTitle(ytitle);


   for(int jj=0;jj< Nuniverses; jj++){
      hthin[ii][jj]   = new TH1D(Form("hthin_%s_%d",  nuhel[ii],jj),Form("%s flux corrected based on thin target data, univ. #%i", nulabel[ii], jj), NbinsE,emin,emax);
      hthin[ii][jj]->SetXTitle(xtitle);
      hthin[ii][jj]->SetYTitle(ytitle);

    
    hpipCpipinc[ii][jj]   = new TH1D(Form("hpipCpipninc_%s_%d",  nuhel[ii],jj),Form("%s flux corrected based on NA61 pip data, univ. #%i", nulabel[ii], jj), NbinsE,emin,emax);
      hpipCpipinc[ii][jj]->SetXTitle(xtitle);
      hpipCpipinc[ii][jj]->SetYTitle(ytitle);


      hmesinc[ii][jj]   = new TH1D(Form("hmesinc_%s_%d",  nuhel[ii],jj),Form("%s flux corrected based on thin target incident meson data, univ. #%i", nulabel[ii], jj), NbinsE,emin,emax);
      hmesinc[ii][jj]->SetXTitle(xtitle);
      hmesinc[ii][jj]->SetYTitle(ytitle);


     /* hmipp[ii][jj]   = new TH1D(Form("hmipp_%s_%d",  nuhel[ii],jj),Form("%s flux corrected based on MIPP NuMI target data, univ. #%i", nulabel[ii], jj),NbinsE,emin,emax);
      hmipp[ii][jj]->SetXTitle(xtitle);
      hmipp[ii][jj]->SetYTitle(ytitle);  */  

      hatt[ii][jj]    = new TH1D(Form("hatt_%s_%d",   nuhel[ii],jj),Form("%s flux corrected for attenuation in target, univ. #%i", nulabel[ii], jj),NbinsE,emin,emax);
      hatt[ii][jj]->SetXTitle(xtitle);
      hatt[ii][jj]->SetYTitle(ytitle);

      hothers[ii][jj] = new TH1D(Form("hothers_%s_%d",nuhel[ii],jj),Form("%s flux corrected for other effects, univ. #%i", nulabel[ii], jj),NbinsE,emin,emax);
      hothers[ii][jj]->SetXTitle(xtitle);
      hothers[ii][jj]->SetYTitle(ytitle);

      htotal[ii][jj]  = new TH1D(Form("htotal_%s_%d", nuhel[ii],jj),Form("%s flux corrected for all effects, univ. #%i", nulabel[ii], jj),NbinsE,emin,emax);
      htotal[ii][jj]->SetXTitle(xtitle);
      htotal[ii][jj]->SetYTitle(ytitle);

      htotalabs[ii][jj]  = new TH1D(Form("htotalabs_%s_%d", nuhel[ii],jj),Form("%s flux corrected for total absorption, univ. #%i", nulabel[ii], jj),NbinsE,emin,emax);
      htotalabs[ii][jj]->SetXTitle(xtitle);
      htotalabs[ii][jj]->SetYTitle(ytitle);

      httpcpi[ii][jj]  = new TH1D(Form("httpcpi_%s_%d", nuhel[ii],jj),Form("%s flux corrected for pcpi, univ. #%i", nulabel[ii], jj),NbinsE,emin,emax);
      httpcpi[ii][jj]->SetXTitle(xtitle);
      httpcpi[ii][jj]->SetYTitle(ytitle);
   
      httpck[ii][jj]  = new TH1D(Form("httpck_%s_%d", nuhel[ii],jj),Form("%s flux corrected for pck, univ. #%i", nulabel[ii], jj),NbinsE,emin,emax);
      httpck[ii][jj]->SetXTitle(xtitle);
      httpck[ii][jj]->SetYTitle(ytitle);

      httncpi[ii][jj]  = new TH1D(Form("httncpi_%s_%d", nuhel[ii],jj),Form("%s flux corrected for ncpi, univ. #%i", nulabel[ii], jj),NbinsE,emin,emax);
      httncpi[ii][jj]->SetXTitle(xtitle);
      httncpi[ii][jj]->SetYTitle(ytitle);

      httpcnu[ii][jj]  = new TH1D(Form("httpcnu_%s_%d", nuhel[ii],jj),Form("%s flux corrected for pcnu, univ. #%i", nulabel[ii], jj),NbinsE,emin,emax);
      httpcnu[ii][jj]->SetXTitle(xtitle);
      httpcnu[ii][jj]->SetYTitle(ytitle);

      httnua[ii][jj]  = new TH1D(Form("httnua_%s_%d", nuhel[ii],jj),Form("%s flux corrected for nua, univ. #%i", nulabel[ii], jj),NbinsE,emin,emax);
      httnua[ii][jj]->SetXTitle(xtitle);
      httnua[ii][jj]->SetYTitle(ytitle);
    }
  }


 for(int ii=0;ii<Nnuhel;ii++){
    for(int jj=0;jj<Nuniverses;jj++){
      hthin[ii][jj]->Sumw2();
     // hmipp[ii][jj]->Sumw2();
      hatt[ii][jj]->Sumw2();
      hothers[ii][jj]->Sumw2();
      htotal[ii][jj]->Sumw2();
      htotalabs[ii][jj]->Sumw2();
    
      httpcpi[ii][jj]->Sumw2();
      httpck[ii][jj]->Sumw2();
      httncpi[ii][jj]->Sumw2();
      httpcnu[ii][jj]->Sumw2();
      hpipCpipinc[ii][jj]->Sumw2();
      hmesinc[ii][jj]->Sumw2();
      httnua[ii][jj]->Sumw2();

    }
    }

  
  //Loading ntuples:
  TChain* chain_evts   = new TChain("dk2nuTree");  
  TChain* chain_meta   = new TChain("dkmetaTree");  
  bsim::Dk2Nu*  dk2nu  = new bsim::Dk2Nu;  
  bsim::DkMeta* dkmeta = new bsim::DkMeta;
  
  std::cout<<" Adding ntuple at: "<<inputFile<<std::endl;  

  std::ifstream ifs(inputFile);
  //  ifs.open(inputFile);
  std::string line;
  int counter = 0;
  while (ifs.good()) {
    getline(ifs,line);
    if(line.find(".root")>10000)continue;
    chain_evts->Add(line.c_str());   
    if(counter==0)chain_meta->Add(line.c_str());
    std::cout<<"Entering: "<<line<<std::endl;
    counter++;
  }
  ifs.close();

  chain_evts->Add(inputFile);
  chain_evts->SetBranchAddress("dk2nu",&dk2nu);
  int nentries  = chain_evts->GetEntries();

  chain_meta->Add(inputFile);
  chain_meta->SetBranchAddress("dkmeta",&dkmeta);
  chain_meta->GetEntry(0); //all entries are the same     
  
  std::string detname = "UserPosition";
  if(doing_precalculated_pos){
    detname = (dkmeta->location)[idet].name;
  }
  std::cout<<"=> Doing the analysis for: "<< detname <<std::endl; 

 // std::vector<double> vwgt_mipp_pi;
 // std::vector<double> vwgt_mipp_K;
  std::vector<double> vwgt_abs;
  std::vector<double> vwgt_att;
  std::vector<double> vwgt_ttpCpi;
  std::vector<double> vwgt_ttpCk;
  std::vector<double> vwgt_ttnCpi;
  std::vector<double> vwgt_ttpCnu;
  std::vector<double> vwgt_ttnua;
  std::vector<double> vwgt_ttpipinc; 
  std::vector<double> vwgt_ttmesinc;
  std::vector<double> vwgt_oth;
 
  std::cout<<"N of entries: "<<nentries<<std::endl;
  double fluxWGT  = 0;
  double nuenergy = 0;

  //for(int ii=0;ii<nentries;ii++){  
    for(int ii=0;ii<nentries;ii++){
    if(ii%1000==0)std::cout<<ii/1000<<" k evts"<<std::endl;
   // std::cout<<ii<<"evts"<<std::endl;

   // vwgt_mipp_pi.clear();  
   // vwgt_mipp_K.clear();  
    vwgt_abs.clear();  
    vwgt_att.clear();  
    vwgt_ttpCpi.clear();  
    vwgt_ttpCk.clear();  
    vwgt_ttnCpi.clear();  
    vwgt_ttpCnu.clear();  
    vwgt_ttpipinc.clear();   
    vwgt_ttmesinc.clear();
    vwgt_ttnua.clear();  
    vwgt_oth.clear(); 
     
    chain_evts->GetEntry(ii);     
    makerew->calculateWeights(dk2nu,dkmeta);
    if(doing_precalculated_pos){
      fluxWGT = ( (dk2nu->nuray)[idet].wgt )*(dk2nu->decay.nimpwt)/pi;   
      nuenergy = (dk2nu->nuray)[idet].E; 
    }
    else{
      nuweight->calculate_weight(dk2nu);
      fluxWGT  = (nuweight->wgt)*(dk2nu->decay.nimpwt)/pi;
      nuenergy = nuweight->enu;
    }

    int nuidx = idx_hel(dk2nu->decay.ntype);
        
    if(nuidx<0){
      std::cout<<"=> Wrong neutrino file"<<std::endl;
    }
    hnom[nuidx]->Fill(nuenergy,fluxWGT);  
    hcv[nuidx]->Fill(nuenergy,fluxWGT*makerew->GetCVWeight());  
    //Get the vector  of universes per reweighter: 12 reweighters
   // vwgt_mipp_pi = makerew->GetWeights("MIPPNumiPionYields");
   // vwgt_mipp_K  = makerew->GetWeights("MIPPNumiKaonYields"); 
    vwgt_abs     = makerew->GetWeights("TotalAbsorption");
    vwgt_att     = makerew->GetWeights("TargetAttenuation");
    vwgt_ttpCpi  = makerew->GetWeights("ThinTargetpCPion");
    vwgt_ttpCk   = makerew->GetWeights("ThinTargetpCKaon");
    vwgt_ttnCpi  = makerew->GetWeights("ThinTargetnCPion");
    vwgt_ttpCnu  = makerew->GetWeights("ThinTargetpCNucleon");
    vwgt_ttpipinc= makerew->GetWeights("ThinTargetpipCpip");   
    vwgt_ttmesinc= makerew->GetWeights("ThinTargetMesonIncident");
    vwgt_ttnua   = makerew->GetWeights("ThinTargetnucleonA");
    vwgt_oth     = makerew->GetWeights("Other"); 

    for(int jj=0;jj<Nuniverses;jj++){
      double wgt_thin = vwgt_ttpCpi[jj]*vwgt_ttpCk[jj]*vwgt_ttnCpi[jj]*vwgt_ttpCnu[jj]*vwgt_ttnua[jj];
     // double wgt_mipp = vwgt_mipp_pi[jj]*vwgt_mipp_K[jj];
      //double wgt_na61 = vwgt_ttpipinc[jj]; //unused
      double wgt_att = vwgt_att[jj]*vwgt_abs[jj];
      //double wgt_mes = vwgt_ttmesinc[jj];   // unused
   
      hthin[nuidx][jj]->Fill(nuenergy,fluxWGT*wgt_thin);
      hmesinc[nuidx][jj]->Fill(nuenergy,fluxWGT*vwgt_ttmesinc[jj]);
      hpipCpipinc[nuidx][jj]->Fill(nuenergy,fluxWGT*vwgt_ttpipinc[jj]); 
      //hmipp[nuidx][jj]->Fill(nuenergy,fluxWGT*wgt_mipp);
      hatt[nuidx][jj]->Fill(nuenergy,fluxWGT*vwgt_att[jj]);
      hothers[nuidx][jj]->Fill(nuenergy,fluxWGT*vwgt_oth[jj]);
      htotal[nuidx][jj]->Fill(nuenergy,fluxWGT*wgt_thin*wgt_att*vwgt_ttpipinc[jj]*vwgt_oth[jj]*vwgt_ttmesinc[jj]); 
      htotalabs[nuidx][jj]->Fill(nuenergy,fluxWGT*vwgt_abs[jj]);
      httpcpi[nuidx][jj]->Fill(nuenergy,fluxWGT*vwgt_ttpCpi[jj]);
      httpck[nuidx][jj]->Fill(nuenergy,fluxWGT*vwgt_ttpCk[jj]);
      httncpi[nuidx][jj]->Fill(nuenergy,fluxWGT*vwgt_ttnCpi[jj]);
      httpcnu[nuidx][jj]->Fill(nuenergy,fluxWGT*vwgt_ttpCnu[jj]);
      httnua[nuidx][jj]->Fill(nuenergy,fluxWGT*vwgt_ttnua[jj]);


    }
  }

  std::cout<<"storing general histos"<<std::endl;
  fOut->cd();
  
  fOut->mkdir("nom"); 
  for(int ii=0;ii<Nnuhel;ii++){
    fOut->mkdir(Form("%s_na61pionInc",nuhel[ii]));
    fOut->mkdir(Form("%s_targetmesonInc",nuhel[ii]));
    fOut->mkdir(Form("%s_attenuation",nuhel[ii]));
    fOut->mkdir(Form("%s_others",nuhel[ii]));
    fOut->mkdir(Form("%s_total",nuhel[ii]));
    fOut->mkdir(Form("%s_absorption",nuhel[ii]));
    fOut->mkdir(Form("%s_targetpcpi",nuhel[ii]));
    fOut->mkdir(Form("%s_targetpck",nuhel[ii]));
    fOut->mkdir(Form("%s_targetncpi",nuhel[ii]));
    fOut->mkdir(Form("%s_targetpcnu",nuhel[ii]));
    fOut->mkdir(Form("%s_targetnuA",nuhel[ii]));
} 
  std::cout<<"Done making the output files"<<std::endl;

  for(int ii=0;ii<Nnuhel;ii++){
    fOut->cd("nom");   
    hnom[ii]->Write();
    hcv[ii]->Write();
    for(int jj=0;jj< Nuniverses; jj++){
      fOut->cd(Form("%s_targetmesonInc" ,nuhel[ii]));  hmesinc[ii][jj]->Write();
      fOut->cd(Form("%s_na61pionInc" ,nuhel[ii]));  hpipCpipinc[ii][jj]->Write();
      fOut->cd(Form("%s_attenuation"  ,nuhel[ii]));  hatt[ii][jj]->Write();
      fOut->cd(Form("%s_others"       ,nuhel[ii]));  hothers[ii][jj]->Write();
      fOut->cd(Form("%s_total"        ,nuhel[ii]));  htotal[ii][jj]->Write();
      fOut->cd(Form("%s_absorption",nuhel[ii]));  htotalabs[ii][jj]->Write();
      fOut->cd(Form("%s_targetpcpi",nuhel[ii]));  httpcpi[ii][jj]->Write();
      fOut->cd(Form("%s_targetpck" ,nuhel[ii]));  httpck[ii][jj]->Write();
      fOut->cd(Form("%s_targetncpi",nuhel[ii]));  httncpi[ii][jj]->Write();
      fOut->cd(Form("%s_targetpcnu",nuhel[ii]));  httpcnu[ii][jj]->Write();
      fOut->cd(Form("%s_targetnuA",nuhel[ii]));  httnua[ii][jj]->Write();
      
    }
  }
  
  //Releasing memory:
  makerew->resetInstance();
  fOut->Close();
  
  std::cout<<"End of run()"<<std::endl;

}
  
int idx_hel(int pdgcode){
  int idx = -1;
  if(pdgcode ==  14)idx = 0;
  if(pdgcode == -14)idx = 1;
  if(pdgcode ==  12)idx = 2;
  if(pdgcode == -12)idx = 3;
  return idx;
}
void usage(){
  std::cout<<"This script calculates the flux at one position in the NuMI beamline: "<<std::endl;
  std::cout<<"Using a precalculated detector positions:"<<std::endl;
  std::cout<<"  bin/doReweight_dk2nu [inputFile] [outputFile] [optionsFile] [idet]"<<std::endl;
  std::cout<<"Using a user input position:"<<std::endl;
  std::cout<<"  bin/doReweight_dk2nu [inputFile] [outputFile] [optionsFile] [xpos] [ypos] [zpos]"<<std::endl;
  std::cout<<"  "<<std::endl;
  std::cout<<"Inputs  "<<std::endl;
  std::cout<<"[inputFile] : g4numi ntuple in dk2nu/dkmeta format (v6 minerva branch is recommended)"<<std::endl;
  std::cout<<"[outputFile] : user definied output file name."<<std::endl;
  std::cout<<"[optionsFile] :xml file with the ppfx input parameters (look at ${PPFX_DIR}/script/input_default.xml)"<<std::endl;
  std::cout<<"[idet] : index of the precalculated detector (look at the location.name in the dkmeta tree of the g4numi ntuple)"<<std::endl;
  std::cout<<"[xpos], [ypos], [zpos] : position (cm) respect to the MC NuMI coordinate system to calculate the flux"<<std::endl;
  std::cout<<"  "<<std::endl;
}
//////////////

#ifndef __CINT__
int main(int argc, const char* argv[]){

  if(argc==5){   
    doReweight_dk2nu(argv[1],argv[2],argv[3],argv[4],"none","none");    
  }
  else if(argc==7){
    doReweight_dk2nu(argv[1],argv[2],argv[3],argv[4],argv[5],argv[6]);
  }
  else{
    usage();
    exit (1);
    }
  
  return 0;
}
#endif
