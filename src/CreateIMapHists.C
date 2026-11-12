//#include "NeutrinoEventChain.h"
//#include "CreateHists.h"

#include "CommonIMapIncludes.h"
#include "FillIMapHists.h"
#include <iostream>
#include <getopt.h>

using namespace std;

// declarations of functions found further down in this file
int CreateHists(const char* output_filename, const char* filename, int elow, int ehigh, int nu_id, bool NA49cuts, bool NA61cuts); 
 //Replaced MIPPcuts by NA6cuts in the CreateHists function
void make_directories(TFile *f);
void name_hists(HistList * hists, TFile * out_file);
void scale_hists(HistList * hists, double total_weight=1);
void write_hists(HistList * hists, TFile * out_file);

int main( int argc, char *argv[]) 
{

  //! Default parameters
  std::vector<std::string> par;
  par.push_back("CreateHists.cpp");
  par.push_back("limiteduni.root");
  par.push_back("/pnfs/nova/persistent/stash/flux/g4numi/v6r1/me000z200i/*_0006.root");
  
//  par.push_back("");   BHUMIKA
  par.push_back("0");
  par.push_back("120");
  par.push_back("14");
  par.push_back("0");//NA49
  par.push_back("0");//MIPP, now it is na61 cut

  //! Set user parameters
  for( int i=0; i<argc; ++i){
    if (string(argv[i]) != string("-1")){
      par.at(i) = argv[i];
    }
  }
  
  int min = atoi(par[3].c_str());
  int max = atoi(par[4].c_str());
  int nu_id = atoi(par[5].c_str());
  bool NA49cuts = atoi(par[6].c_str());
  //bool MIPPcuts = atoi(par[7].c_str());
  bool NA61cuts = atoi(par[7].c_str());


  std::cout<<"min E   : "<<min<<"GeV"<<std::endl;
  std::cout<<"max E   : "<<max<<"GeV"<<std::endl;
  std::cout<<"nuhel   : "<<nu_id<<std::endl;  
  std::cout<<"TT cut  : "<<NA49cuts<<std::endl;
  std::cout<<"NA61 cut: "<<NA61cuts<<std::endl;
    
    //int CreateHists(const char* output_filename, const char* filename, int elow, int ehigh, int nu_id=56, bool NA49cuts=false, bool MIPPcuts=false))
  return CreateHists(par[1].c_str(), par[2].c_str(), min, max, nu_id, NA49cuts, NA61cuts);
}


int CreateHists(const char* output_filename, const char* filename, int elow, int ehigh, int nu_id=14, bool NA49cuts=false, bool NA61cuts=false)
{
        //single array histos	
	int nvol = IMap::nvol;
	//if(getenv("MODE")=="REF"||getenv("MODE")=="OPT") nvol = IMap::nvoldune;  BHUMIKA, as using nova
	vector<TH2D *> hmat(IMap::npop);
	vector<TH2D *> hvol(IMap::npop);
	vector<TH2D *> hmatbkw(IMap::npop);
	vector<TH2D *> hxfpt_tot(IMap::npop);   //This has a different phase space
	vector<TH1F *> henergytotal(IMap::npop);
	vector<TH1F *> hkepop_tot(IMap::npop);   //KEPOP AND TMPOP, WHAT ARE THESE EXACTLY?
	vector<TH1F *> htmpop_tot(IMap::npop);

	//double array histos
	vector< vector<TH2D *> > hxfpt(IMap::npop, vector<TH2D*>(IMap::npop));
	vector< vector<TH1F *> > henergymaterial(IMap::nspecialmat, vector<TH1F*>(IMap::npop));
	vector< vector<TH1F *> > henergyvolume(nvol, vector<TH1F*>(IMap::npop));
	vector< vector<TH1F *> > hkepop(IMap::npop, vector<TH1F*>(IMap::npop));
	vector< vector<TH1F *> > htmpop(IMap::npop, vector<TH1F*>(IMap::npop));
  
	//individual histos
	TH2D * h_in_vs_mat = NULL;

	//collect histos in a list
	HistList hists = {hmat, hvol, hmatbkw, hxfpt_tot, henergytotal, 
			  hkepop_tot, htmpop_tot, hxfpt, henergymaterial, 
			  henergyvolume, hkepop, htmpop, h_in_vs_mat};

	//make root file
	cout << "Making output file " << output_filename << endl;
	TFile *out_file= new TFile(output_filename,"RECREATE");
	
	//make directories in the output root file.	
	cout << "Making directories in " << output_filename << endl;
	make_directories(out_file);
	
	//name hists
	cout << "Naming Hists" << endl;
	name_hists(& hists, out_file);

	//Make an object of many neutrino ancestral lines	from the ntuples
	cout << "Opening flux ntuples" << endl;

	TChain* tdk2nu   = new TChain("dk2nuTree");  
	TChain* tdkmeta   = new TChain("dkmetaTree");  
	tdk2nu->Add(filename);
	tdkmeta->Add(filename);
	//Loop over the neutrinos and fill histos
	cout << "Filling histograms" << endl;	
	FillIMapHistsOpts opts;
	opts.elow=elow; opts.ehigh=ehigh; opts.nuid=nu_id;
	opts.cut_thintarget=NA49cuts; opts.cut_na61=NA61cuts; opts.cut_mipp=NA61cuts; // same 7th CLI arg drives MIPP cut again, restored for MINERvA/MIPP comparison
	double total_weight=FillIMapHists(tdk2nu, tdkmeta, &hists, &opts);
	
	cout<< "Total weight: "<<total_weight<<endl;
	//Scale the histos
	cout << "Scaling Histos" << endl;
	scale_hists(& hists,total_weight);
 
	//Write hists
 	cout << "Writing Histos" << endl;	
	write_hists(& hists, out_file);
    
	out_file->Close();
	
	return 0;

} //END INT MAP


//Function definitions for:
//make_directories
//name_hists
//scale_hists
//write_hists







// ----------------------------------------------
// Make directories in the output root file
// ----------------------------------------------
void make_directories(TFile *f){
//Projectile
//---particle_i
//------Energy_Materials
//------Energy_Volumes
//Produced
//---particle_i
  using namespace IMap;
    f->mkdir("Projectile");
    f->mkdir("Produced");

    for(int j=0;j<IMap::npop;j++){
      f->cd("Projectile");
      // Projectile/{particle}
      gDirectory->mkdir(popparticle[j].c_str());
      gDirectory->cd(popparticle[j].c_str());
      // Projectile/{particle}/Energy_Materials
      gDirectory->mkdir("Energy_Materials");
      // Projectile/{particle}/Energy_Volumes
      gDirectory->mkdir("Energy_Volumes");
    }
    
    for(int j=0;j<IMap::npop;j++){
      f->cd("Produced");
      // Produced/{particle}
      gDirectory->mkdir(popparticle[j].c_str());
    }
}

// ----------------------------------------------
// Name all the histograms
// ----------------------------------------------
void name_hists(HistList * hists, TFile * out_file){
		
  char namefile[100];
  using namespace IMap;
    // ----------------------------------------------
    // Make Histos in Projectile/
    // ----------------------------------------------
    // list of histos in Projectile
    // hmat[9]
    // hvol[9]
    // henergytotal[9]
    // henergymaterial[5][9]
    // henergyvolume[193][9]
    
    
    // Will be in Projectile/h_in_vs_mat
    hists->_h_in_vs_mat = new TH2D("h_in_vs_mat",";material ; incident particle",50,0,50,50,0,50);

    hists->_h_nint_vs_enu = new TH2D("h_nint_vs_enu","all interactions;neutrino energy (GeV); number of interactions",120,0,120,20,-0.5,19.5);
    hists->_h_nint_vs_enu_cuts=new TH2D("h_nint_vs_enu_cuts", "interactions not covered by HP;neutrino energy (GeV); number of interactions",120,0,120,20,-0.5,19.5) ;
    
    hists->_h_aveint_vs_enu_thin_pCpion     = new TH1D("h_aveint_vs_enu_thin_pCpion","thin target pC->pion;neutrino energy (GeV); number of interactions",120,0,120);
    hists->_h_aveint_vs_enu_thin_pCkaon     = new TH1D("h_aveint_vs_enu_thin_pCkaon","thin target pC->kaon;neutrino energy (GeV); number of interactions",120,0,120);
    hists->_h_aveint_vs_enu_thin_nCpion     = new TH1D("h_aveint_vs_enu_thin_nCpion","thin target nC->pion;neutrino energy (GeV); number of interactions",120,0,120);
    hists->_h_aveint_vs_enu_thin_pCnucleon  = new TH1D("h_aveint_vs_enu_thin_pCnucleon","thin target pC->nucleon;neutrino energy (GeV); number of interactions",120,0,120);
    hists->_h_aveint_vs_enu_thin_mesoninc   = new TH1D("h_aveint_vs_enu_thin_mesoninc","thin target meson incident;neutrino energy (GeV); number of interactions",120,0,120);
    hists->_h_aveint_vs_enu_pipCpip   = new TH1D("h_aveint_vs_enu_pipCpip","thin target Pi plus incident;neutrino energy (GeV); number of interactions",120,0,120);    //BHUMIKA  
  hists->_h_aveint_vs_enu_thin_nucleona   = new TH1D("h_aveint_vs_enu_thin_nucleona","thin target nuclen-A;neutrino energy (GeV); number of interactions",120,0,120);
    hists->_h_aveint_vs_enu_others          = new TH1D("h_aveint_vs_enu_others","Others;neutrino energy (GeV); number of interactions",120,0,120);
    hists->_h_aveint_vs_enu_tot             = new TH1D("h_aveint_vs_enu_tot","total;neutrino energy (GeV); number of interactions",120,0,120);
    hists->_h_nuflux   = new TH1D("h_nuflux","#nu flux;neutrino energy (GeV); #nu/m^{2}",120,0,120);

    hists->_h_occ_xfpt_pc_pip = new TH2D("h_occ_xfpt_pc_pip",";xF; pT (GeV/c)",100,0,1.0,50,0,1.0);
    hists->_h_occ_xfpt_pc_kp = new TH2D("h_occ_xfpt_pc_kp",";xF; pT (GeV/c)",100,0,1.0,50,0,1.0);
    hists->_h_hpwgt_xfpt_pc_pip = new TH2D("h_hpwgt_xfpt_pc_pip",";xF; pT (GeV/c)",100,0,1.0,50,0,1.0);
    hists->_h_hpwgt_xfpt_pc_kp = new TH2D("h_hpwgt_xfpt_pc_kp",";xF; pT (GeV/c)",100,0,1.0,50,0,1.0);
    
    for(int k=0;k<IMap::npop;k++)
    {
        sprintf(namefile,"Projectile/%s",popparticle[k].c_str());
        out_file->cd(namefile);
        
        sprintf(namefile,"hmats%s",popparticle[k].c_str());
        // Projectile/{particle}/hmats{particle}
        hists->_hmat[k]=new TH2D(namefile,popparticle[k].c_str(),50,0,50,50,0,50);
			
        sprintf(namefile,"hvols%s",popparticle[k].c_str());
        // Projectile/{particle}/hvols{particle}
        hists->_hvol[k]=new TH2D(namefile,popparticle[k].c_str(),50,0,50,50,0,50);
        
        
        sprintf(namefile,"henergy%s",popparticle[k].c_str());
        // Projectile/{particle}/henergy{particle}
        hists->_henergytotal[k] = new TH1F(namefile,popparticle[k].c_str(),1000,0,200);
        
        hists->_henergytotal[k]->GetXaxis()->SetTitle("Energy (GeV)");
        hists->_henergytotal[k]->GetYaxis()->SetTitle("Fraction of ancestry");
        
        sprintf(namefile,"Projectile/%s/Energy_Materials",popparticle[k].c_str());
        out_file->cd(namefile);
        
        for(int i=0;i<IMap::nspecialmat;i++)
	 {
            sprintf(namefile,"henergy%s%s",popparticle[k].c_str(),matlist[i].c_str());
            // Projectile/{particle}/Energy_Materials/henergy{particle}{material}
            hists->_henergymaterial[i][k] = new TH1F(namefile,"",1000,0,200);
            sprintf(namefile,"%s in %s",popparticle[k].c_str(),matlist[i].c_str());
						hists->_henergymaterial[i][k]->SetTitle(namefile);
            hists->_henergymaterial[i][k]->GetXaxis()->SetTitle("Energy (GeV)");
            hists->_henergymaterial[i][k]->GetYaxis()->SetTitle("Event Counts");
	
	 }
        
        sprintf(namefile,"Projectile/%s/Energy_Volumes",popparticle[k].c_str());
        out_file->cd(namefile);
        
        for(int i=0;i<nvol;i++)
        {
            
	  sprintf(namefile,"henergy%s%s",popparticle[k].c_str(),IMap::volume[i].c_str());
            // Projectile/{particle}/Energy_Volumes/henergy{particle}{volume}
            hists->_henergyvolume[i][k]=new TH1F(namefile,"",1000,0,200);
            sprintf(namefile,"%s in %s",popparticle[k].c_str(),IMap::volume[i].c_str());
						hists->_henergyvolume[i][k]->SetTitle(namefile);
            hists->_henergyvolume[i][k]->GetXaxis()->SetTitle("Energy (GeV)");
            hists->_henergyvolume[i][k]->GetYaxis()->SetTitle("Event Counts");
				}
    }
    
    
    // ----------------------------------------------
    // Make histos in Produced/
    // ----------------------------------------------
    // list of histos in Produced:
    // hmatbkw[9]
    // hkepop_tot[9]
    // ktmpop_tot[9]
    // hxfpt_tot[9]
    // hkepop[9][9]
    // htmpop[9][9]
    // hxfpt[9][9]
    for(int j=0;j<IMap::npop;j++)
    {
        sprintf(namefile,"Produced/%s",popparticle[j].c_str());
        out_file->cd(namefile);
        sprintf(namefile,"hmatbkg%s",popparticle[j].c_str());
        // Produced/{particle}/hmatbkg{particle}
        hists->_hmatbkw[j]=new TH2D(namefile,popparticle[j].c_str(),50,0,50,50,0,50);
        sprintf(namefile,"hke%s",popparticle[j].c_str());
        // Produced/{particle}/hke{particle}
        hists->_hkepop_tot[j]=new TH1F(namefile,"",1000,0,200);
        sprintf(namefile,"htm%s",popparticle[j].c_str());
        // Produced/{particle}/htm{particle}
        hists->_htmpop_tot[j]=new TH1F(namefile,"",1000,0,200);
        sprintf(namefile,"hxfpt%s",popparticle[j].c_str());
        // Produced/{particle}/hxfpt{particle}
        hists->_hxfpt_tot[j]=new TH2D(namefile,"",500,-1,1,100,0,2);
        
        for(int k=0;k<IMap::npop;k++)
        {
            sprintf(namefile,"hke%s_%s",popparticle[j].c_str(),popparticle[k].c_str());
            // Produced/{particle}/hke{particle}_{other particle}
            hists->_hkepop[j][k]=new TH1F(namefile,"",1000,0,200);
            sprintf(namefile,"htm%s_%s",popparticle[j].c_str(),popparticle[k].c_str());
            // Produced/{particle}/htm{particle}_{other particle}
            hists->_htmpop[j][k]=new TH1F(namefile,"",1000,0,200);
            sprintf(namefile,"hxfpt%s_%s",popparticle[j].c_str(),popparticle[k].c_str());
            // Produced/{particle}/hxfpt{particle}_{other particle}
            hists->_hxfpt[j][k]=new TH2D(namefile,"",500,-1,1,100,0,2);
        }
    }
    
}


// ----------------------------------------------
// Scale all histograms by the total weight
// ----------------------------------------------
void scale_hists(HistList * hists, double total_weight){
		
    hists->_h_in_vs_mat->Scale(1./total_weight);
    hists->_h_nint_vs_enu->Scale(1./total_weight);
    hists->_h_nint_vs_enu_cuts->Scale(1./total_weight);

    hists->_h_aveint_vs_enu_thin_pCpion->Divide(hists->_h_nuflux); 
    hists->_h_aveint_vs_enu_thin_pCkaon->Divide(hists->_h_nuflux); 
    hists->_h_aveint_vs_enu_thin_nCpion->Divide(hists->_h_nuflux); 
    hists->_h_aveint_vs_enu_thin_pCnucleon->Divide(hists->_h_nuflux); 
    hists->_h_aveint_vs_enu_thin_mesoninc->Divide(hists->_h_nuflux); 
    hists->_h_aveint_vs_enu_pipCpip->Divide(hists->_h_nuflux);    //Bhumika
    hists->_h_aveint_vs_enu_thin_nucleona->Divide(hists->_h_nuflux); 
    hists->_h_aveint_vs_enu_others->Divide(hists->_h_nuflux); 
    hists->_h_aveint_vs_enu_tot->Divide(hists->_h_nuflux); 

    hists->_h_aveint_vs_enu_thin_pCpion->Scale(1./pival);
    hists->_h_aveint_vs_enu_thin_pCkaon->Scale(1./pival);
    hists->_h_aveint_vs_enu_thin_nCpion->Scale(1./pival);
    hists->_h_aveint_vs_enu_thin_pCnucleon->Scale(1./pival);
    hists->_h_aveint_vs_enu_thin_mesoninc->Scale(1./pival);
    hists->_h_aveint_vs_enu_pipCpip->Scale(1./pival);            //BHUMIKA
    hists->_h_aveint_vs_enu_thin_nucleona->Scale(1./pival);
    hists->_h_aveint_vs_enu_others->Scale(1./pival);
    hists->_h_aveint_vs_enu_tot->Scale(1./pival);

    hists->_h_occ_xfpt_pc_pip->Scale(1.0/total_weight);
    hists->_h_occ_xfpt_pc_kp->Scale(1.0/total_weight);
    hists->_h_hpwgt_xfpt_pc_pip->Scale(1.0/total_weight);
    hists->_h_hpwgt_xfpt_pc_kp->Scale(1.0/total_weight);

    hists->_h_hpwgt_xfpt_pc_pip->Divide(hists->_h_occ_xfpt_pc_pip);
    hists->_h_hpwgt_xfpt_pc_kp->Divide(hists->_h_occ_xfpt_pc_kp);

    // loop over all particle-specific histograms
    for(int j=0;j<IMap::npop;j++)
    {
      hists->_hmat[j]->Scale(1./total_weight);
      hists->_hvol[j]->Scale(1./total_weight);
      hists->_henergytotal[j]->Scale(1./total_weight);
      hists->_hmatbkw[j]->Scale(1./total_weight);
      hists->_hkepop_tot[j]->Scale(1./total_weight);
      hists->_htmpop_tot[j]->Scale(1./total_weight);
      hists->_hxfpt_tot[j]->Scale(1./total_weight);
      for(int k=0;k<IMap::npop;k++)
      {
        hists->_hkepop[j][k]->Scale(1./total_weight);
        hists->_htmpop[j][k]->Scale(1./total_weight);
        hists->_hxfpt[j][k]->Scale(1./total_weight);
      } // for(int k=0;k<9;k++)
      // Loop over materials
      for(int m=0;m<IMap::nspecialmat;m++)
      {
        hists->_henergymaterial[m][j]->Scale(1./total_weight);
      } // for(int m=0;m<5;m++)
    } // for(int j=0;j<9;j++)

}

// ----------------------------------------------
// Write all histos to file
// ----------------------------------------------
void write_hists(HistList * hists, TFile * out_file){
    
  char namefile[100];
  using namespace IMap;
  for(int j=0;j<IMap::npop;j++)
    {
        sprintf(namefile,"Projectile/%s",popparticle[j].c_str());
        out_file->cd(namefile);
        hists->_hmat[j]->Write();
        hists->_hvol[j]->Write();
        hists->_henergytotal[j]->Write();
        
        gDirectory->cd("Energy_Materials");
        for(int m=0;m<IMap::nspecialmat;m++)
        {
            
            hists->_henergymaterial[m][j]->Write();
            
        }
        gDirectory->cd("..");
        gDirectory->cd("Energy_Volumes");
        for(int m=0;m<nvol;m++)
        {
            hists->_henergyvolume[m][j]->Write();
        }
    }
    
    
  for(int j=0;j<IMap::npop;j++)
    {
        sprintf(namefile,"Produced/%s",popparticle[j].c_str());
        out_file->cd(namefile);
        hists->_hmatbkw[j]->Write();
        hists->_hkepop_tot[j]->Write();
        hists->_htmpop_tot[j]->Write();
        hists->_hxfpt_tot[j]->Write();
        
        for(int k=0;k<IMap::npop;k++)
        {
            hists->_hkepop[j][k]->Write();
            hists->_htmpop[j][k]->Write();
            hists->_hxfpt[j][k]->Write();
        }
    }

  // Save the h_in_vs_mat histogram in Projectile/h_in_vs_mat
  out_file->cd("Projectile");
  hists->_h_in_vs_mat->Write();

  hists->_h_nint_vs_enu->Write();
  hists->_h_nint_vs_enu_cuts->Write();

  hists->_h_aveint_vs_enu_thin_pCpion->Write(); 
  hists->_h_aveint_vs_enu_thin_pCkaon->Write();
  hists->_h_aveint_vs_enu_thin_nCpion->Write();
  hists->_h_aveint_vs_enu_thin_pCnucleon->Write();
  hists->_h_aveint_vs_enu_thin_mesoninc->Write(); 
  hists->_h_aveint_vs_enu_pipCpip->Write();              //Bhumika
  hists->_h_aveint_vs_enu_thin_nucleona->Write();
  hists->_h_aveint_vs_enu_others->Write();
  hists->_h_aveint_vs_enu_tot->Write();
  hists->_h_nuflux->Write();

  hists->_h_occ_xfpt_pc_pip->Write();
  hists->_h_occ_xfpt_pc_kp->Write();
  hists->_h_hpwgt_xfpt_pc_pip->Write();
  hists->_h_hpwgt_xfpt_pc_kp->Write();

}
