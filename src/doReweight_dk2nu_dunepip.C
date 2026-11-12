
#include "dk2nu/tree/dkmeta.h"
#include "dk2nu/tree/dk2nu.h"
#include "MakeReweight.h"
#include "NuWeight.h"
#include <string>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <iomanip>

#include "TH1D.h"

//Some constants::

const int NbinsE   = 120;
const double emin  =   0.;
const double emax  = 120.;
const int Nnuhel   = 4;
const char* nuhel[Nnuhel] = {"numu","numubar","nue","nuebar"};
const double pi    = 3.14159265358979323846;

using namespace NeutrinoFluxReweight;

class nu_g4numi;
class Numi2Pdg;

int idx_hel(int pdgdcode);
//! returns string with "proper" xrootd enabled path if reading from PNFS (dCache)
std::string change_to_xrootd_path(std::string temp);


/*!
 * Run the reweighting for a single file (inputFile) for
 * particular MIPP covariance matrix given in input.xml file. 
 */
void doReweight_dk2nu(const char* inputFile, const char* outputFile, const char* optionsFile){ 

std::vector<Double_t>var_bins_numu;
std::vector<Double_t>var_bins_nue;
/*
for(int i =0;i<6/0.5;i++)var_bins_numu.push_back(i*0.5);
for(int i = 0;i<(8-6)/1;i++)var_bins_numu.push_back(6.0+i*1);
for(int i = 0;i<(20-8)/4;i++)var_bins_numu.push_back(8.0+i*4);
for(int i = 0;i<(40-20)/20;i++)var_bins_numu.push_back(20+i*20);
for(int i = 0;i<(100-40)/60;i++)var_bins_numu.push_back(40+i*40);
var_bins_numu.push_back(100.0);

//now for nue
for(int i = 0;i<10/2;i++)var_bins_nue.push_back(i*2);
for(int i = 0;i<(20-10)/10;i++)var_bins_nue.push_back(10+i*10);
for(int i = 0;i<(100-20)/80;i++)var_bins_nue.push_back(20+i*80);
var_bins_nue.push_back(100.0);
  
*/
  std::vector<Double_t>var_bins;
  //1 GeV bin width upto 30 GeV
 // for(int i = 0;i<5/1;i++)var_bins.push_back(i*1);
  //2 GeV bin width upto 120 GeV
 // for (int i = 0;i<(120-5)/2;i++)var_bins.push_back(5+i*2);
for(int i = 0;i<120/0.5;i++)var_bins.push_back(i*0.5);
//for(int i = 0;i<120/1;i++)var_bins_nue.push_back(i*1);    
  const char* thisDir = getenv("PPFX_DIR");
  const char* OutputDir=thisDir;
 

  // pnfs - dCache access via XROOTD logic
  std::string infile = change_to_xrootd_path( std::string(inputFile) ) ;		//inputFile);

  std::string outfile =std::string(outputFile);
  //Create 2 outputs for both near detector and far detector
  std::string near = "DUNEND";  // or LBNFND
  std::string far  = "DUNEFD";  // or LBNFFD
  int detindex = 1; // 1 = Near Detector, 2 = Far Detector
  std::size_t foundnear = outfile.find(near);
  std::size_t foundfar = outfile.find(far);
  if(foundnear != std::string::npos)detindex=1;
  if(foundfar != std::string::npos)detindex=2;
  
  std::cout<<" The detector index is "<<detindex<<std::endl;
  std::cout<< "Instance of MakeReweight()" <<std::endl;
  MakeReweight* makerew = MakeReweight::getInstance();
  makerew->SetOptions(optionsFile); 
  
  std::cout<<"Making an output file to store histograms"<<std::endl;
  TFile* fOut = new TFile(outputFile,"recreate");
  std::cout<<"File name: "<<fOut->GetName()<<std::endl;
  
  fOut->mkdir("nom");
  for(int ii=0;ii<Nnuhel;ii++){
    fOut->mkdir(Form("%s_na61pionInc",nuhel[ii]));
    fOut->mkdir(Form("%s_thintarget",nuhel[ii]));
  //  fOut->mkdir(Form("%s_mippnumi",nuhel[ii]));
    fOut->mkdir(Form("%s_attenuation",nuhel[ii]));
    fOut->mkdir(Form("%s_others",nuhel[ii]));
    fOut->mkdir(Form("%s_total",nuhel[ii]));
    fOut->mkdir(Form("%s_totabsorption",nuhel[ii]));
    fOut->mkdir(Form("%s_tgtabsorption",nuhel[ii]));
    fOut->mkdir(Form("%s_targetpcpi",nuhel[ii]));
    fOut->mkdir(Form("%s_targetpck",nuhel[ii]));      
    fOut->mkdir(Form("%s_targetncpi",nuhel[ii]));    
    fOut->mkdir(Form("%s_targetpcnu",nuhel[ii]));      
    fOut->mkdir(Form("%s_targetmesonInc",nuhel[ii]));      
    fOut->mkdir(Form("%s_targetnuA",nuhel[ii]));     
    
    
  }
  std::cout<<"Done making the output file"<<std::endl;
  
  int Nuniverses = makerew->GetNumberOfUniversesUsed();
  
  TH1D* hnom[Nnuhel];
  TH1D* hcv[Nnuhel];
  TH1D* hthin[Nnuhel][Nuniverses];
  TH1D* hpipCpipinc[Nnuhel][Nuniverses];
//  TH1D* hmipp[Nnuhel][Nuniverses];
  TH1D* hatt[Nnuhel][Nuniverses];
  TH1D* hothers[Nnuhel][Nuniverses];
  TH1D* htotal[Nnuhel][Nuniverses];
  TH1D* htotab[Nnuhel][Nuniverses];
  TH1D* htgtab[Nnuhel][Nuniverses];
  TH1D* httpcpi[Nnuhel][Nuniverses];
  TH1D* httpck[Nnuhel][Nuniverses]; 
  TH1D* httncpi[Nnuhel][Nuniverses]; 
  TH1D* httpcnu[Nnuhel][Nuniverses];  
  TH1D* httmesinc[Nnuhel][Nuniverses];  
  TH1D* httnua[Nnuhel][Nuniverses]; 
  
  for(int ii=0;ii<Nnuhel;ii++){
 // std::vector<Double_t>var_bins;
  //clone the binning vector
 // if((ii==0)||(ii==1))var_bins = var_bins_numu;
 // if((ii==2)||(ii==3))var_bins = var_bins_nue;
  std::cout<<" ii "<<var_bins.at(ii);
    hnom[ii] = new TH1D(Form("hnom_%s",nuhel[ii]),"",var_bins.size()-1,&var_bins[0]);
    hcv[ii]  = new TH1D(Form("hcv_%s",nuhel[ii]),"",var_bins.size()-1,&var_bins[0]);
    for(int jj=0;jj< Nuniverses; jj++){
      hthin[ii][jj]   = new TH1D(Form("hthin_%s_%d",  nuhel[ii],jj),"",var_bins.size()-1,&var_bins[0]);
   hpipCpipinc[ii][jj]= new TH1D(Form("hpipCpipninc_%s_%d",  nuhel[ii],jj),"",var_bins.size()-1,&var_bins[0]);
  //    hmipp[ii][jj]   = new TH1D(Form("hmipp_%s_%d",  nuhel[ii],jj),"",var_bins.size()-1,&var_bins[0]);
      hatt[ii][jj]    = new TH1D(Form("hatt_%s_%d",   nuhel[ii],jj),"",var_bins.size()-1,&var_bins[0]);
      hothers[ii][jj] = new TH1D(Form("hothers_%s_%d",nuhel[ii],jj),"",var_bins.size()-1,&var_bins[0]);
      htotal[ii][jj]  = new TH1D(Form("htotal_%s_%d", nuhel[ii],jj),"",var_bins.size()-1,&var_bins[0]);
      htotab[ii][jj]  = new TH1D(Form("habsorption_%s_%d",nuhel[ii],jj),"",var_bins.size()-1,&var_bins[0]);
      htgtab[ii][jj]  = new TH1D(Form("htargetabsorption_%s_%d",nuhel[ii],jj),"",var_bins.size()-1,&var_bins[0]);
      httpcpi[ii][jj]  = new TH1D(Form("htargetpcpi_%s_%d",nuhel[ii],jj),"",var_bins.size()-1,&var_bins[0]);
      httpck[ii][jj]  = new TH1D(Form("htargetpcK_%s_%d",nuhel[ii],jj),"",var_bins.size()-1,&var_bins[0]);  
      httncpi[ii][jj]  = new TH1D(Form("htargetncpi_%s_%d",nuhel[ii],jj),"",var_bins.size()-1,&var_bins[0]);
      httpcnu[ii][jj]  = new TH1D(Form("htargetpcnu_%s_%d",nuhel[ii],jj),"",var_bins.size()-1,&var_bins[0]);      
      httmesinc[ii][jj]  = new TH1D(Form("htargetmesonInc_%s_%d",nuhel[ii],jj),"",var_bins.size()-1,&var_bins[0]);  
      httnua[ii][jj]  = new TH1D(Form("htargetnucleonA_%s_%d",nuhel[ii],jj),"",var_bins.size()-1,&var_bins[0]);    
    
    
    
    }
  }
  
    
    for(int ii=0;ii<Nnuhel;ii++){
    for(int jj=0;jj<Nuniverses;jj++){
      hthin[ii][jj]->Sumw2();
      hpipCpipinc[ii][jj]->Sumw2();
    //  hmipp[ii][jj]->Sumw2();
      hatt[ii][jj]->Sumw2();
      hothers[ii][jj]->Sumw2();
      htotal[ii][jj]->Sumw2();
      htotab[ii][jj]->Sumw2();
      htgtab[ii][jj]->Sumw2();
      httpcpi[ii][jj]->Sumw2();
      httpck[ii][jj]->Sumw2(); 
      httncpi[ii][jj]->Sumw2();
      httpcnu[ii][jj]->Sumw2();
      httmesinc[ii][jj]->Sumw2();
      httnua[ii][jj]->Sumw2();    
    
    }
    }
  
  //Loading ntuples:
  TChain* chain_evts   = new TChain("dk2nuTree");  
  TChain* chain_meta   = new TChain("dkmetaTree");  
  bsim::Dk2Nu*  dk2nu  = new bsim::Dk2Nu;  
  bsim::DkMeta* dkmeta = new bsim::DkMeta;
  
  std::cout<<" Adding ntuple at: "<< infile <<std::endl;	//inputFile<<std::endl;

  chain_evts->Add(infile.c_str());
  chain_evts->SetBranchAddress("dk2nu",&dk2nu);
  int nentries  = chain_evts->GetEntries();

  chain_meta->Add(infile.c_str());
  chain_meta->SetBranchAddress("dkmeta",&dkmeta);
  chain_meta->GetEntry(0); //all entries are the same     
  
 // std::vector<double> vwgt_mipp_pi;
 // std::vector<double> vwgt_mipp_K;
  std::vector<double> vwgt_abs;
  std::vector<double> vwgt_att;
  std::vector<double> vwgt_ttpCpi;
  std::vector<double> vwgt_ttpCk;
  std::vector<double> vwgt_ttnCpi;
  std::vector<double> vwgt_ttpCnu;
  std::vector<double> vwgt_ttnua;
  std::vector<double> vwgt_ttpipinc; //bhu
  std::vector<double> vwgt_ttmesinc;
  std::vector<double> vwgt_oth;
 
  std::cout<<"N of entries: "<<nentries<<std::endl;
  
  for(int ii=0;ii<nentries;ii++){  
    if(ii%1000==0)std::cout<<ii/1000<<" k evts"<<std::endl;
   // vwgt_mipp_pi.clear();  
   // vwgt_mipp_K.clear();  
    vwgt_abs.clear();  
    vwgt_att.clear();  
    vwgt_ttpCpi.clear();  
    vwgt_ttpCk.clear();  
    vwgt_ttnCpi.clear();  
    vwgt_ttpCnu.clear();  
     vwgt_ttpipinc.clear();   //bhu
    vwgt_ttmesinc.clear();
    vwgt_ttnua.clear();  
    vwgt_oth.clear(); 
     
    chain_evts->GetEntry(ii);     
 //   std::cout<<"I got the chain events"<<std::endl; 
   makerew->calculateWeights(dk2nu,dkmeta);
  //  std::cout<<"I have calculated the weights"<<std::endl;
    double fluxWGT = ( (dk2nu->nuray)[detindex].wgt )*(dk2nu->decay.nimpwt)/pi;	// 3.1416;
    int nuidx = idx_hel(dk2nu->decay.ntype);
    double nuenergy = (dk2nu->nuray)[detindex].E;
    double vz = (dk2nu->decay.vz);
    double size = (dk2nu->ancestor).size();
    double startz = (dk2nu->ancestor)[size-2].startz; 
    
    if(nuidx<0){
      std::cout<<"=> Wrong neutrino file"<<std::endl;
    }
    hnom[nuidx]->Fill(nuenergy,fluxWGT);
    hcv[nuidx]->Fill(nuenergy,fluxWGT*makerew->GetCVWeight());

  //  vwgt_mipp_pi = makerew->GetWeights("MIPPNumiPionYields");
  //  vwgt_mipp_K  = makerew->GetWeights("MIPPNumiKaonYields"); 
    vwgt_abs     = makerew->GetWeights("TotalAbsorption");
    vwgt_att     = makerew->GetWeights("TargetAttenuation");
    vwgt_ttpCpi  = makerew->GetWeights("ThinTargetpCPion");
    vwgt_ttpCk   = makerew->GetWeights("ThinTargetpCKaon");
    vwgt_ttnCpi  = makerew->GetWeights("ThinTargetnCPion");
    vwgt_ttpCnu  = makerew->GetWeights("ThinTargetpCNucleon");
    vwgt_ttpipinc= makerew->GetWeights("ThinTargetpipCpip");   //BHU
    vwgt_ttmesinc= makerew->GetWeights("ThinTargetMesonIncident");
    vwgt_ttnua   = makerew->GetWeights("ThinTargetnucleonA");
    vwgt_oth     = makerew->GetWeights("Other"); 
    
    
    for(int jj=0;jj<Nuniverses;jj++){
      double wgt_thin = vwgt_ttpCpi[jj]*vwgt_ttpCk[jj]*vwgt_ttnCpi[jj]*vwgt_ttpCnu[jj]*vwgt_ttmesinc[jj]*vwgt_ttnua[jj];
    //  double wgt_mipp = vwgt_mipp_pi[jj]*vwgt_mipp_K[jj];
      double wgt_att = vwgt_att[jj]*vwgt_abs[jj];
      double wgt_na61 = vwgt_ttpipinc[jj]; //BHU
      hpipCpipinc[nuidx][jj]->Fill(nuenergy,fluxWGT*vwgt_ttpipinc[jj]);
      hthin[nuidx][jj]->Fill(nuenergy,fluxWGT*wgt_thin);
  //    hmipp[nuidx][jj]->Fill(nuenergy,fluxWGT*wgt_mipp);
      hatt[nuidx][jj]->Fill(nuenergy,fluxWGT*wgt_att);
      hothers[nuidx][jj]->Fill(nuenergy,fluxWGT*vwgt_oth[jj]);
      htotal[nuidx][jj]->Fill(nuenergy,fluxWGT*wgt_thin*wgt_att*vwgt_oth[jj]*vwgt_ttpipinc[jj]);
      htotab[nuidx][jj]->Fill(nuenergy,fluxWGT*vwgt_abs[jj]);
      htgtab[nuidx][jj]->Fill(nuenergy,fluxWGT*vwgt_att[jj]);
      httpcpi[nuidx][jj]->Fill(nuenergy,fluxWGT*vwgt_ttpCpi[jj]);
      httpck[nuidx][jj]->Fill(nuenergy,fluxWGT*vwgt_ttpCk[jj]); 
      httncpi[nuidx][jj]->Fill(nuenergy,fluxWGT*vwgt_ttnCpi[jj]);
      httpcnu[nuidx][jj]->Fill(nuenergy,fluxWGT*vwgt_ttpCnu[jj]);
      httmesinc[nuidx][jj]->Fill(nuenergy,fluxWGT*vwgt_ttmesinc[jj]);
      httnua[nuidx][jj]->Fill(nuenergy,fluxWGT*vwgt_ttnua[jj]);
      // std::cout<<"NA61 weights are:"<< vwgt_ttpipinc[jj];    
      //if(wgt_att>2)std::cout<<" hatt entry is "<<nuenergy<<" "<<vwgt_att[jj]<<std::endl;
      //std::cout<<"Target Attenuation Weights "<<wgt_att<<" vwgt quantity "<<vwgt_oth[jj]<<std::endl;
      
      //std::cout<<"Target attenuation "<<vwgt_att[jj]<<" Total Absorption "<<vwgt_abs[jj]<<std::endl;
     // if(vwgt_abs[jj]<1)std::cout<<vz<<" "<<startz<<" "<<vwgt_att[jj]<<std::endl;
    }  
  }
  
  std::cout<<"storing general histos"<<std::endl;
  fOut->cd();
  for(int ii=0;ii<Nnuhel;ii++){
    fOut->cd("nom");
    hnom[ii]->Write();
    hcv[ii]->Write();
    for(int jj=0;jj< Nuniverses; jj++){
      fOut->cd(Form("%s_thintarget",nuhel[ii]));  hthin[ii][jj]->Write();
      fOut->cd(Form("%s_na61pionInc",nuhel[ii]));  hpipCpipinc[ii][jj]->Write();
     // fOut->cd(Form("%s_na61pionInc",nuhel[ii]));  hpipCpipinc[ii][jj]->Write();
   //   fOut->cd(Form("%s_mippnumi",nuhel[ii]));    hmipp[ii][jj]->Write();
      fOut->cd(Form("%s_attenuation",nuhel[ii])); hatt[ii][jj]->Write();
      fOut->cd(Form("%s_others",nuhel[ii]));      hothers[ii][jj]->Write();
      fOut->cd(Form("%s_total",nuhel[ii]));       htotal[ii][jj]->Write();  
      fOut->cd(Form("%s_totabsorption",nuhel[ii]));htotab[ii][jj]->Write();
      fOut->cd(Form("%s_tgtabsorption",nuhel[ii]));htgtab[ii][jj]->Write();
      fOut->cd(Form("%s_targetpcpi",nuhel[ii]));httpcpi[ii][jj]->Write();
      fOut->cd(Form("%s_targetpck",nuhel[ii]));httpck[ii][jj]->Write();      
      fOut->cd(Form("%s_targetncpi",nuhel[ii]));httncpi[ii][jj]->Write();     
      fOut->cd(Form("%s_targetpcnu",nuhel[ii]));httpcnu[ii][jj]->Write();      
      fOut->cd(Form("%s_targetmesonInc",nuhel[ii]));httmesinc[ii][jj]->Write();      
      fOut->cd(Form("%s_targetnuA",nuhel[ii]));httnua[ii][jj]->Write();      
           
    }
  }
  fOut->Close();
  makerew->resetInstance(); //Releasing memory: 
  std::cout<<"End of run()"<<std::endl;

} // end doReweight

std::string change_to_xrootd_path(std::string temp){
//  pnfs access via root streaming is now not OK. Must switch to xrootd. Dec 2021 -Pierce Weatherly
//  std::string input_flux_dir_fe =  "/pnfs/dune/"+disk+"/users/"+input_user+"/fluxfiles/g4lbne/"+version+"/"+physics_list+"/"+macro+"/"+current+"/flux/";
//  std::string rootxdstr = "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/";
//  std::string input_flux_dir =  rootxdstr+disk+"/users/"+input_user+"/fluxfiles/g4lbne/"+version+"/"+physics_list+"/"+macro+"/"+current+"/flux/";
//  std::cout<<"\n\t Using input directory: "<<input_flux_dir<<std::endl;
//  if(on_grid) {
//    input_flux_dir = getenv("_CONDOR_SCRATCH_DIR");
//    input_flux_dir += "/";
//    std::cout<<"Running on grid, so getting flux files from local disk: "<<input_flux_dir<<std::endl;
//  }

  // pnfs - dCache access via XROOTD logic
  bool enable_xrootd = true;
  std::string temp_inputFile = temp;
  std::string rootxdstr = "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/";
  std::size_t index = temp_inputFile.find("/pnfs/");
  if( index == std::string::npos || index != 0) enable_xrootd = false; // cannot find /pnfs/
  else if(index == 0) enable_xrootd = true; // if the path starts with /pnfs/ , then go ahead and use xrootd.
//  if(enable_xrootd && temp_inputFile.find(grid_dir) != string::npos) enable_xrootd = false;  // if on-grid, leave alone
//  if(enable_xrootd && temp_inputFile.find("/data/") != string::npos) enable_xrootd = false;  // if grabing from dune/data or whatever, leave along; make sure not in pnfs
 // if(enable_xrootd && temp_inputFile.find("root:/") != string::npos) enable_xrootd = false;  // if it looks like it is xrootd already
  if(enable_xrootd){
        temp_inputFile.replace(0,6,rootxdstr);
        printf("\n Replacing string for inputFIle %s\n\t with XROOTD path: %s\n", temp.c_str(), temp_inputFile.c_str() );
  }else printf("\n Input File not on /pnfs/ (dCache); File path: %s\n", temp_inputFile.c_str() );
  return temp_inputFile;
}

  
int idx_hel(int pdgcode){
  int idx = -1;
  if(pdgcode ==  14)idx = 0;
  if(pdgcode == -14)idx = 1;
  if(pdgcode ==  12)idx = 2;
  if(pdgcode == -12)idx = 3;
  return idx;
}

//////////////

#ifndef __CINT__
int main(int argc, const char* argv[]){
  
  doReweight_dk2nu(argv[1],argv[2],argv[3]);
  return 0;
}
#endif
