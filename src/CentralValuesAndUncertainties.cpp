

#include "CentralValuesAndUncertainties.h"
#include <set>
#include <string>
#include <iostream>
#include <sstream>
#include <TMath.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "TRandom3.h"
#include "TMatrixD.h"
//using namespace std;

namespace NeutrinoFluxReweight{ 
  
  CentralValuesAndUncertainties* CentralValuesAndUncertainties::instance = 0;

  CentralValuesAndUncertainties::CentralValuesAndUncertainties(){
    //! Setting the seed to zero. 
    r3=new TRandom3(0);   //Random number generator class in ROOT, seed is selected to be zero here for reproducibility 
    r3_pip=new TRandom3(300);   //New generator for the newly introduced reweighter
    baseSeed = 0;
  }
  
  void CentralValuesAndUncertainties::readFromXML(const char* filename){
    using boost::property_tree::ptree;
    ptree top;    //Object of type ptree named top
    
    read_xml(filename,top,2); // option 2 removes comment strings
    
    //! Reading the uncorrelated parameters:
    ptree& uncorrelated = top.get_child("pars.uncorrelated");    
    ptree::iterator it = uncorrelated.begin();
    for(; it!=uncorrelated.end(); ++it){
      // it->first is the name
      // it->second is the child property tree
      double cv=it->second.get<double>("cv");
      double err=it->second.get<double>("err");
      Parameter p(it->first,cv);
      CentralValuesAndUncertainties::addUncorrelated(p,err);
    }

    //! Reading the uncorrelated list:
    ptree& uncorrelated_list = top.get_child("pars.uncorrelated_list");
    it = uncorrelated_list.begin();
    for(; it!=uncorrelated_list.end(); ++it){
      // it->first is the name, for us FIRST IS THE NAME OF THE INTERACTION
      // it->second is the child property tree, SECOND IS CV/ERROR/COVMX

      std::string cvs_string=it->second.get<std::string>("cvs");
      std::string errs_string=it->second.get<std::string>("errs");
      double cv=0;
      int ii=0;
      std::string name;
      std::stringstream ss(cvs_string);
      std::vector<Parameter> tmp_par;
      while(ss >> cv){
	std::stringstream sID;
	sID << ii;
	std::string nameID = sID.str();
	name = it->first + "_" + nameID;
	Parameter p(name,cv);
	tmp_par.push_back(p);
	ii++;
      };
      
      double err=0;
      std::stringstream sserr(errs_string);
      ii=0;
      while(sserr >> err){
	CentralValuesAndUncertainties::addUncorrelated(tmp_par[ii],err);
	ii++;
      }     
      
    }
    
    //! Reading the correlated parameters:
    ptree& correlated = top.get_child("pars.correlated");
    it = correlated.begin();
    for(; it!=correlated.end(); ++it){
      // it->first is the name
      // it->second is the child property tree
      std::string cvs_string=it->second.get<std::string>("cvs");
      std::string covmx_string=it->second.get<std::string>("covmx");

      //filling the Parameter table:
      std::stringstream ss(cvs_string);
      double cv=0;
      int ii=0;
      std::string name;
      ParameterTable ptable;
      while(ss >> cv){
	std::stringstream sID;
	sID << ii;
	std::string nameID = sID.str();
	name = it->first + "_" +nameID;
	Parameter p(name,cv);
	ptable.setParameter(p);
	ii++;
      };
      TMatrixD mcov(ii,ii);
      //filling the matrix:
      std::stringstream ssmx(covmx_string);
      double err;
      int idx = 0;
      while(ssmx >> err){
	mcov(idx/ii,idx%ii) = err;
	idx++;
      }
      CentralValuesAndUncertainties::addCorrelated(ptable,mcov);
    }

    //Auxiliar parameter.
    //cv=1.0 and error 100%:
    Parameter par_aux("aux_parameter",1.0);
    CentralValuesAndUncertainties::addUncorrelated(par_aux,1.0);

  }
  
  void CentralValuesAndUncertainties::addUncorrelated(Parameter& cv_par, double uncertainty){
    uncorrelated_pars.setParameter(cv_par);
    uncorrelated_errors[cv_par.first] = uncertainty;
  }
  
  void CentralValuesAndUncertainties::addCorrelated(ParameterTable& cv_pars, TMatrixD& cov_mx){
    correlated_par_tables.push_back(cv_pars);
    covariance_matrices.push_back(cov_mx);
  } 
  
  void CentralValuesAndUncertainties::setBaseSeed(int val){
    baseSeed = val;
  }
  ParameterTable CentralValuesAndUncertainties::calculateParsForUniverse(int universe){

    // WARNING: if you add new data here it needs to be added at the end
    
    //If universe = -1, then it is the central value:
    double cvfactor = 1.0;
    if(universe==-1)cvfactor = 0.0;
    int univ_seed = baseSeed + universe;
    r3->SetSeed(univ_seed);    
    r3_pip->SetSeed(univ_seed+300);   
    ParameterTable ptable;
    
    //! We are going to use 100% correlated bin-to-bin for systematic errors in thin target data:
    double sigma_pc_pip = r3->Gaus(0.0,1.0);
    double sigma_pc_pim = r3->Gaus(0.0,1.0);
    double sigma_pc_kap = r3->Gaus(0.0,1.0);
    double sigma_pc_kam = r3->Gaus(0.0,1.0);
    double sigma_pc_p   = r3->Gaus(0.0,1.0);
    double sigma_pc_n   = r3->Gaus(0.0,1.0);
    // double sigma_pipC_pip = r3_pip->Gaus(0.0, 1.0);   // unused
    const boost::interprocess::flat_map<std::string, double>& table_uncorr_pars = uncorrelated_pars.getMap();
    boost::interprocess::flat_map<std::string, double>::const_iterator it = table_uncorr_pars.begin();

    
  for(;it!=table_uncorr_pars.end();++it){
      double sigma= r3->Gaus(0.0,1.0);
      double sigma2 = r3_pip->Gaus(0.0,1.0);   //For pipCpip uncorrelated
      
      //redefining sigma for 100% correlation:
      if((it->first).find("ThinTarget_pC_pip_sys")<10)sigma = sigma_pc_pip;
      if((it->first).find("ThinTarget_pC_pim_sys")<10)sigma = sigma_pc_pim;
      if((it->first).find("ThinTargetLowxF_pC_kap_sys")<10)sigma = sigma_pc_kap;
      if((it->first).find("ThinTargetLowxF_pC_kam_sys")<10)sigma = sigma_pc_kam;
      if((it->first).find("ThinTarget_pC_p_sys")<10)  sigma = sigma_pc_p;
      if((it->first).find("ThinTarget_pC_n_sys")<10)  sigma = sigma_pc_n;
     
     //Redefination for NA61 data
        if((it->first).find("ThinTarget_pipC_pip_stat")<10)sigma = sigma2;
  

     	if((it->first).find("ThinTarget_pipC_pim_stat")<10)sigma = sigma2;
        if((it->first).find("ThinTarget_pipC_kp_stat")<10)sigma = sigma2;
        if((it->first).find("ThinTarget_pipC_km_stat")<10)sigma = sigma2;
        if((it->first).find("ThinTarget_pipC_p_stat")<10)sigma = sigma2;
        if((it->first).find("ThinTarget_pipC_k0s_stat")<10)sigma = sigma2;
        if((it->first).find("ThinTarget_pipC_lam_stat")<10)sigma = sigma2;
        if((it->first).find("ThinTarget_pipC_alam_stat")<10)sigma = sigma2;
       
      double new_val = it->second  + cvfactor*sigma*uncorrelated_errors[it->first];
      Parameter p(it->first,new_val);
      ptable.setParameter(p);      
    
    }
    TDecompSVD *decomp;
    
    for(size_t ii=0;ii<covariance_matrices.size();++ii){
      
      decomp=new TDecompSVD(covariance_matrices[ii],0.0); 
      
      bool isDecomposed=decomp->Decompose();
      if(isDecomposed){std::cout<<"The matrix is being decomposed using SVD decomposition"<<std::endl;}
      TMatrixD U = decomp->GetU();
      TVectorD S = decomp->GetSig();
      TMatrixD Vt = decomp->GetV();


      int nmat = S.GetNoElements();
      TVectorD vsigma(nmat);
      for(int jj=0;jj<nmat;jj++){
       vsigma[jj]= cvfactor*(r3_pip->Gaus(0.0, 1.0));
      }
      TMatrixD S_matrix(nmat, nmat);
      S_matrix.Zero();
      for (Int_t i = 0; i < nmat; ++i) {
        S_matrix(i, i) = S[i];
       }
     TVectorD sqrtS(nmat);
       for (Int_t i = 0; i < nmat; ++i) {
        sqrtS[i] = TMath::Sqrt(S[i]);
   }

    //The MATRIX FOR NEW VACTOR IS
   TMatrixD Sq_matrix(nmat, nmat);
    Sq_matrix.Zero();
    for (Int_t i = 0; i < nmat; ++i) {
      Sq_matrix(i, i) = sqrtS[i];}
    


   TVectorD vecDShift =  ((U * Sq_matrix)* vsigma); 

      
      const boost::interprocess::flat_map<std::string, double>& tb = (correlated_par_tables[ii]).getMap();
     boost::interprocess::flat_map<std::string, double>::const_iterator it_tb = tb.begin();

       for(;it_tb != tb.end();++it_tb){
	std::string tmp_name = it_tb->first;
	std::string snID = tmp_name.substr((it_tb->first).rfind("_")+1,(it_tb->first).length());
	std::stringstream ssID(snID);
	int nID;
	ssID >> nID;
	double new_val = it_tb->second + vecDShift[nID];	
	Parameter p(it_tb->first,new_val); 
	if(isDecomposed)ptable.setParameter(p);
      
      
  }  
  }      // covariance matrix      
    return ptable; 
    
  }       //universe  
  
  ParameterTable CentralValuesAndUncertainties::getCVPars(){      

    ParameterTable ptableCV;
    
    const boost::interprocess::flat_map<std::string, double>& table_pars=uncorrelated_pars.getMap();
    boost::interprocess::flat_map<std::string, double>::const_iterator it = table_pars.begin();
    for(;it!=table_pars.end();it++){
      Parameter p(it->first,it->second); 
      ptableCV.setParameter(p);      
    }
    
    //Correlated:
    for(size_t ii=0;ii<covariance_matrices.size();ii++){
    const boost::interprocess::flat_map<std::string, double>& corr_table_pars = correlated_par_tables[ii].getMap();
      it = corr_table_pars.begin();
      for(;it != corr_table_pars.end();it++){
	Parameter p(it->first,it->second); 
	ptableCV.setParameter(p);
      }
      
    }
    
    return ptableCV;
    
  }
  
  CentralValuesAndUncertainties* CentralValuesAndUncertainties::getInstance(){
    if (instance == 0) instance = new CentralValuesAndUncertainties;
    return instance;
  }

}

