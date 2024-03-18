#ifndef WREMNANTS_MUONCALIB_H
#define WREMNANTS_MUONCALIB_H


#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"


namespace wrem {

// Standalone code to make multiD histograms for the sagitta fit

  using namespace ROOT;
  using namespace ROOT::VecOps;

  //dxy_significance
  RVecF dxy_significance(RVecF Muon_dxy, RVecF Muon_dxyErr){
    return abs(Muon_dxy)/Muon_dxyErr;
  }

  ROOT::VecOps::RVec<TRandom3*> gausInSlot(unsigned int nSlot, unsigned int event) {
    ROOT::VecOps::RVec<TRandom3*> rans;
    rans.emplace_back( new TRandom3(4357) );
    rans.emplace_back( new TRandom3(3951) );
    rans.emplace_back( new TRandom3(5193) );
    rans.emplace_back( new TRandom3(9361) );
    return rans;
  }

  //simply for debugging
  bool plotInSlot(unsigned int nSlot, unsigned int event, ROOT::VecOps::RVec<TRandom3*> rg) {
    TRandom3* uni = new TRandom3(event);
    std::cout << nSlot << ">>>>" << rg[0]->Gaus(uni->Uniform(-1, 1),1) << ":" << rg[1]->Gaus(uni->Uniform(24, 28),4) << ":" << rg[2]->Gaus(uni->Uniform(39, 45),4) << ":" << rg[3]->Gaus(uni->Uniform(86, 94),2) << std::endl;
    return true;
  }
  
  using PairTuple=std::tuple<int,int,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float>;

  PairTuple muPair(unsigned int nSlot, RVecF Muon_pt, RVecI Muon_charge, RVecF Muon_eta, RVecF Muon_phi, RVecF Muon_dxy, RVecF Muon_dz) {
    RVec<PairTuple> pairs;
    PairTuple temp, pair_to_return;
    float rest_mass = 0.105658; // muMass = 0.105658 GeV
    float firstPt_reco, secondPt_reco, mll_reco;
    
    for(int i=1;i<Muon_pt.size();i++){
      for(int j=0;j<i;j++){
	if(Muon_charge[i]*Muon_charge[j]==-1 && abs(Muon_dxy[i]-Muon_dxy[j])<0.1 && abs(Muon_dz[i]-Muon_dz[j])<0.6){
	  TLorentzVector firstTrack, secondTrack, mother, firstGenTrack, secondGenTrack, motherGen, firstSmearTrack, secondSmearTrack, motherSmear;
	  firstTrack.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], rest_mass);
	  secondTrack.SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], rest_mass);
	  mother = firstTrack + secondTrack;
	  mll_reco = mother.M();
	  firstPt_reco = Muon_pt[i];
	  secondPt_reco = Muon_pt[j];
	    
	  if(75.0<mll_reco && mll_reco<105.0){ //Cut in mll
	    if(Muon_charge[i]==1){
	      temp=make_tuple(i,j,mll_reco,firstPt_reco,secondPt_reco,0.,0.,0.,0.,0.,0.,0., 0.,0.,0.,0.,0.,0.);
	    } else {
	      temp=make_tuple(j,i,mll_reco,secondPt_reco,firstPt_reco,0.,0.,0.,0.,0.,0.,0., 0.,0.,0.,0.,0.,0.);
	    }
	    pairs.push_back(temp);
	  }
	}
      }
    }
    if(pairs.size()==1){
      pair_to_return=pairs.at(0);
    } else if(pairs.size()>1){
      float diff=100.0;
      int best=0;
      for(int i=0;i<pairs.size();i++){
	if(abs(get<2>(pairs.at(i))-91)<diff){
	  diff=(abs(get<2>(pairs.at(i))-91));
	  best=i;
	}
      }
      pair_to_return=pairs.at(best);
    } else {
      pair_to_return=make_tuple(0,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
    }
    return pair_to_return;
  }

  PairTuple muPair(unsigned int nSlot, RVecF Muon_pt, RVecI Muon_charge, RVecF Muon_eta, RVecF Muon_phi, RVecF Muon_dxy, RVecF Muon_dz, RVecI GenPart_status, RVecI GenPart_pdgId, RVecI GenPart_genPartIdxMother, RVecF GenPart_pt, RVecF GenPart_eta, RVecF GenPart_phi, const ROOT::VecOps::RVec<TRandom3*> rGaus, bool isData = false) {  

    RVec<PairTuple> pairs; // <pos_muon_index, neg_muon_index, mll_reco, posPt_reco, negPt_reco, mll_diff_reco, mll_smear, posPt_smear, negPt_smear, mll_diff_smear, mll_gen, posPt_gen, negPt_gen, mll_diff_squared_smear, smear_beta_weight, posPt_smear_beta_val, negPt_smear_beta_val, mll_diff_squared_reco>
    PairTuple temp, pair_to_return;
    float rest_mass = 0.105658; // muMass = 0.105658 GeV
    float firstPt_reco, secondPt_reco, mll_reco, firstPt_smear, secondPt_smear, mll_smear, firstPt_gen, secondPt_gen, mll_gen, firstPt_smear_beta_val, secondPt_smear_beta_val, smear_beta_weight;
    float smear_pt, mean, width, beta=0.999, smear_beta_weight_first_term, smear_beta_weight_second_term;
    
    for(int i=1;i<Muon_pt.size();i++){
      for(int j=0;j<i;j++){
	if(Muon_charge[i]*Muon_charge[j]==-1 && abs(Muon_dxy[i]-Muon_dxy[j])<0.1 && abs(Muon_dz[i]-Muon_dz[j])<0.6){
	  TLorentzVector firstTrack, secondTrack, mother, firstGenTrack, secondGenTrack, motherGen, firstSmearTrack, secondSmearTrack, motherSmear;
	  firstTrack.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], rest_mass);
	  secondTrack.SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], rest_mass);
	  mother = firstTrack + secondTrack;
	  mll_reco = mother.M();
	  firstPt_reco = Muon_pt[i];
	  secondPt_reco = Muon_pt[j];
	    
	  if(75.0<mll_reco && mll_reco<105.0){ //Cut in mll
	    //Gen match
	    bool firstGenMatched = false, secondGenMatched = false;
	    for (int k=0;k<GenPart_eta.size();k++){
	      if (GenPart_status[k]==1 && abs(GenPart_pdgId[k])==13 && GenPart_pdgId[GenPart_genPartIdxMother[k]]==23){ // mu(-) has PDGID 13
		if(pow(pow(GenPart_eta[k]-Muon_eta[i],2) + pow(GenPart_phi[k]-Muon_phi[i],2),0.5)<0.3){
		  firstGenTrack.SetPtEtaPhiM(GenPart_pt[k], GenPart_eta[k], GenPart_phi[k], rest_mass);
		  firstPt_gen = GenPart_pt[k];
		  
		  //smear 1st muon
		  mean = GenPart_pt[k]; //beta = 1 
		  width = (0.004*pow(GenPart_eta[k],2)+0.01)*GenPart_pt[k];
		  firstPt_smear = rGaus[0]->Gaus(mean, width);
		  firstSmearTrack.SetPtEtaPhiM(firstPt_smear, GenPart_eta[k], GenPart_phi[k], rest_mass);
		  //smear_beta_val, weight for beta != 1
		  smear_beta_weight_first_term = TMath::Gaus(firstPt_smear, mean*beta, width) / TMath::Gaus(firstPt_smear, mean, width);
		  firstPt_smear_beta_val = rGaus[1]->Gaus(mean*beta, width);
		  
		  firstGenMatched = true;
		  if(secondGenMatched == true){break;}
		} else if(pow(pow(GenPart_eta[k]-Muon_eta[j],2) + pow(GenPart_phi[k]-Muon_phi[j],2),0.5)<0.3){
		  secondGenTrack.SetPtEtaPhiM(GenPart_pt[k], GenPart_eta[k], GenPart_phi[k], rest_mass);
		  secondPt_gen = GenPart_pt[k];
		  
		  //smear 2nd muon
		  mean = GenPart_pt[k]; //beta = 1
		  width = (0.004*pow(GenPart_eta[k],2)+0.01)*GenPart_pt[k];
		  secondPt_smear = rGaus[2]->Gaus(mean, width);
		  secondSmearTrack.SetPtEtaPhiM(secondPt_smear, GenPart_eta[k], GenPart_phi[k], rest_mass);
		  //smear_beta_val, weight for beta != 1
		  smear_beta_weight_second_term = TMath::Gaus(secondPt_smear, mean*beta, width) / TMath::Gaus(secondPt_smear, mean, width); 
		  secondPt_smear_beta_val = rGaus[3]->Gaus(mean*beta, width);
		  
		  secondGenMatched = true;
		  if(firstGenMatched == true){break;}
		}
	      }
	    }
	    if(firstGenMatched == false || secondGenMatched == false){
	      continue;
	    }
	    
	    motherSmear = firstSmearTrack + secondSmearTrack;
	    mll_smear = motherSmear.M();
	    smear_beta_weight = smear_beta_weight_first_term * smear_beta_weight_second_term;
	    
	    motherGen = firstGenTrack + secondGenTrack;
	    float mll_gen = motherGen.M();
	    float mll_diff_reco = mll_reco - mll_gen;
	    float mll_diff_smear = mll_smear - mll_gen;
	    // save for jacobians
	    float mll_diff_squared_smear = (mll_smear - mll_gen)*(mll_smear - mll_gen);
	    float mll_diff_squared_reco = (mll_reco - mll_gen)*(mll_reco - mll_gen);
            
	    if(Muon_charge[i]==1){
	      temp=make_tuple(i,j,mll_reco,firstPt_reco,secondPt_reco,mll_diff_reco,mll_smear,firstPt_smear,secondPt_smear,mll_diff_smear,mll_gen,firstPt_gen,secondPt_gen,mll_diff_squared_smear,smear_beta_weight,firstPt_smear_beta_val,secondPt_smear_beta_val, mll_diff_squared_reco);
	    } else {
	      temp=make_tuple(j,i,mll_reco,secondPt_reco,firstPt_reco,mll_diff_reco,mll_smear,secondPt_smear,firstPt_smear,mll_diff_smear,mll_gen,secondPt_gen,firstPt_gen,mll_diff_squared_smear,smear_beta_weight,secondPt_smear_beta_val,firstPt_smear_beta_val, mll_diff_squared_reco);
	    }
	    pairs.push_back(temp);
	  }
	}
      }
    }


    if(pairs.size()==1){
      pair_to_return=pairs.at(0);
    } else if(pairs.size()>1){
      float diff=100.0;
      int best=0;
      for(int i=0;i<pairs.size();i++){
	if(abs(get<2>(pairs.at(i))-91)<diff){
	  diff=(abs(get<2>(pairs.at(i))-91));
	  best=i;
	}
      }
      pair_to_return=pairs.at(best);
    } else {
      pair_to_return=make_tuple(0,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
    }
    return pair_to_return;
}

};

#endif
