///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  PYTHIA+HYDJET scanning code
//  Version : response
//  Author: Clayton Bennett
//  Date : 30 May 2022
//  Notes:
//  	- to be used with 
//
//  Version Notes:
//  	~1.0 --> - First version.  Rewrite of the forest_scan_ppMC_lxplus scanner [3/8/22]
//  	~1.1 --> - Fix matching bug.  Wasnt actually cutting on dR.  Also fix bug where
//  			we weren't skipping already-matched muons.
//  	~1.2 --> - Added reco v. ref correlation histograms for measuring resolution
//  		- Added delete file command after event loop
//
//      ~2.0 --> - Reformat to be compatible with PYTHIAHYDJET_scan_V2p0 code, AnalysisSetupV2p0, etc.
//               - Switch flavorID to be just the partonFlavor variable.
//               - allRecoJets-allRecoMuons
//        +2.01 --> - allGenJets-allRecoMuons
//		  +2.11 --> - just allGenJets, no muons
//        +jetRes --> - just calculate jetResoltuion maps
//
//
//      response
//      v1 : - just calculate a TH2D with genJet pt and matchedRecoJetPt, i.e. jet response matrix (for unfolding)
//			 - add a variabe-bin size response histogram that matched the template binnings
//
//      adapt to PYTHIA+HYDJET
//      
//			- adding plots for JER/JES calculation [6-29-22]
//			- fixing filling the allJets histograms for each centrality [7-16-22]
//			- added 500 pt point [7-18-22]
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// general ROOT/C includes
#include <iostream>
#include "TFile.h"
#include "TRandom.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TRandom.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TProfile2D.h"
#include <TF1.h>
#include "assert.h"
#include <fstream>
#include "TMath.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TMath.h"
#include <TNtuple.h>
#include "TChain.h"
#include <TString.h>
#include <TLatex.h>
#include <TCut.h>
#include "TDatime.h"
#include <vector>
#include "TCanvas.h"
#include <dirent.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

// user includes (LXPLUS)
// #include "/afs/cern.ch/user/c/cbennett/CMSSW_10_3_3_patch1/src/myProcesses/hiforest/plugin/eventMap_hiForest.h"
// #include "/afs/cern.ch/user/c/cbennett/CMSSW_10_3_3_patch1/src/JetEnergyCorrections/JetCorrector.h"
// #include "/afs/cern.ch/user/c/cbennett/condorSkim/currentCode/AnalysisSetupV2p0.h"




// user includes (local)
#include "/home/clayton/Analysis/code/myProcesses/hiforest/plugin/eventMap_hiForest.h"
#include "/home/clayton/Analysis/code/JetEnergyCorrections/JetCorrector.h"
#include "/home/clayton/Analysis/code/headers/AnalysisSetupV2p0.h"




double getPtRel(double MuonPt, double MuonEta, double MuonPhi, double JetPt, double JetEta, double JetPhi)
{

	double Muon_Px = MuonPt * TMath::Cos(MuonPhi);
	double Muon_Py = MuonPt * TMath::Sin(MuonPhi);
	double Muon_Pz = MuonPt * TMath::SinH(MuonEta);

	double Jet_Px = JetPt * TMath::Cos(JetPhi);
	double Jet_Py = JetPt * TMath::Sin(JetPhi);
	double Jet_Pz = JetPt * TMath::SinH(JetEta);

	float lj_x = Jet_Px;
	float lj_y = Jet_Py;
	float lj_z = Jet_Pz;

	// absolute values squared
	float lj2 = lj_x * lj_x + lj_y * lj_y + lj_z * lj_z;

	//float lep2 = lep.px()*lep.px()+lep.py()*lep.py()+lep.pz()*lep.pz();
	float lep2 = Muon_Px * Muon_Px + Muon_Py * Muon_Py + Muon_Pz * Muon_Pz;

	// projection vec(mu) to lepjet axis
	float lepXlj = Muon_Px * lj_x + Muon_Py * lj_y + Muon_Pz * lj_z;

	// absolute value squared and normalized
	float pLrel2 = lepXlj * lepXlj / lj2;
	float pTrel2 = lep2 - pLrel2;

	return (pTrel2 > 0) ? std::sqrt(pTrel2) : 0.0;
}


const int getCentBin(int hiBin)
{
	int result = -1;
	if (hiBin >= hiBin_C0_lo && hiBin < hiBin_C0_hi)
		result = 1;
	else if (hiBin >= hiBin_C0_hi && hiBin < hiBin_C1_hi)
		result = 2;
	else if (hiBin >= hiBin_C1_hi && hiBin < hiBin_C2_hi)
		result = 3;
	else if (hiBin >= hiBin_C2_hi && hiBin < hiBin_C3_hi)
		result = 4;
	else
		;
	return result;
}	

const int getJetPtBin(double jetpt)
{
	int result = -1;
	if (jetpt >= jetpt_J0_lo && jetpt < jetpt_J0_hi)
		result = 1;
	else if (jetpt >= jetpt_J0_hi && jetpt < jetpt_J1_hi)
		result = 2;
	else if (jetpt >= jetpt_J1_hi && jetpt < jetpt_J2_hi)
		result = 3;
	else if (jetpt >= jetpt_J2_hi && jetpt < jetpt_J3_hi)
		result = 4;
	else if (jetpt >= jetpt_J3_hi && jetpt < jetpt_J4_hi)
		result = 5;
	else if (jetpt >= jetpt_J4_hi && jetpt < jetpt_J5_hi)
		result = 6;
	else
		;
	return result;
}




double getDr(double b, double c, double y, double z)
{
	// b, y = eta values
	// c, z = phi values
	double result = -1.0;
	double dEta = b - y;
	double dPhi = acos(cos(c - z));
	result = TMath::Sqrt(TMath::Power(dEta, 2) + TMath::Power(dPhi, 2));
	return result;
}

bool passesHiBinCut(int hiBin){
	bool result = true;
	if(hiBin > hiBinCut || hiBin < hiBinCut_low) result = false;
	return result;
}

bool triggerIsOn(int triggerValue, int prescaleValue){

	bool result = true;
	if(triggerValue == 0 || prescaleValue == 0) result = false;
	return result;

}

bool isQualityMuon(double muPt,
		double muEta, 
		int muChi2NDF, 
		double muInnerD0, 
		double muInnerDz, 
		int muMuonHits, 
		int muPixelHits, 
		int muIsTracker, 
		int muStations, 
		int muTrkLayers){

	bool result = true;
	if(muPt < muPtCut ||
			TMath::Abs(muEta) > 2.0 ||
			muChi2NDF == -99 ||
			muChi2NDF > 10 ||
			TMath::Abs(muInnerD0) > 0.2 ||
			TMath::Abs(muInnerDz) > 0.5 ||
			muMuonHits <= 0 ||
			muPixelHits <= 0 ||
			muIsTracker == 0 ||
			muStations  <= 1 ||
			muTrkLayers <= 5)
		result = false;

	return result;
}

bool isWDecayMuon(double muonPt, double jetPt){

	bool result = false;
	if(muonPt/jetPt > 0.75) result = true;
	return result;


}

bool passesRecoJetPthatFilter(double pthat, double jetPt){

	bool result = true;
	if( (pthat < 40 && jetPt > 93.75) || (pthat > 40 && pthat < 300 && jetPt > 93.75 + 1.37*(pthat - 40.)) ) result = false;
	return result;

}






//void PYTHIAHYDJET_scan_response_v1(TString input = "root://cmsxrootd.fnal.gov//store/user/cmbennet/PYTHIAHYDJET_forest_noRecoJetPtCut_23May22/DiJet_pThat-15_TuneCP5_HydjetDrumMB_5p02TeV_Pythia8/crab_PYTHIAHYDJET_forest_noRecoJetPtCut_23May22/220523_170508/0000/HiForestAOD_101.root", TString output = "out.root"){

void PYTHIAHYDJET_scan_response_v1(TString input = "~/Analysis/data/eos/PbPbMC/10July22/HiForestAOD_101.root", TString output = "out.root"){



// JET ENERGY CORRECTIONS
vector<string> Files;
Files.push_back("/home/clayton/Analysis/code/JetEnergyCorrections/Autumn18_HI_V8_MC_L2Relative_AK4PF.txt"); // LOCAL
//Files.push_back("/afs/cern.ch/user/c/cbennett/CMSSW_10_3_3_patch1/src/JetEnergyCorrections/Autumn18_HI_V8_MC_L2Relative_AK4PF.txt"); // LXPLUS
JetCorrector JEC(Files);

// WEIGHT FUNCTIONS

TF1 *fVzWeightFunction = new TF1("fvz", "pol1", -15, 15);
fVzWeightFunction->SetParameters(1.01655,-0.0184946);

TF1 *fCentralityWeightFunction = new TF1("fcent", "pol5", 0, 180);
fCentralityWeightFunction->SetParameters(1.24607,-0.0121385,0.000526902,-7.57068e-06,3.61464e-08,-5.11966e-11);



//  initialize histograms
// JETS

TH2D *h_matchedRecoJetPt_genJetPt[NCentralityIndices];
TH2D *h_matchedRecoJetPt_genJetPt_var[NCentralityIndices];
TH2D *h_matchedRecoJetPtOverGenJetPt_genJetPt[NCentralityIndices][7];

TH2D *h_matchedRawJetPt_genJetPt[NCentralityIndices];
TH2D *h_matchedRawJetPt_genJetPt_var[NCentralityIndices];
TH2D *h_matchedRawJetPtOverGenJetPt_genJetPt[NCentralityIndices][7];

//RooUnfoldResponse Response(NPtBins,ptMin,ptMax);
//TH2D *h_jeteta_flavor;
//TH2D *h_jetphi_flavor[NCentralityIndices];

//TH2D *h_ratio_refpt[4];



// Define histograms
	const int N1 = 7;
	double ptAxis1[N1] = {50,60,80,120,200,300,500};
	const int N2 = 32;
	double ptAxis2[N2] = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,500};



    for(int i = 0; i < NCentralityIndices; i++){

        if(i==0) h_matchedRecoJetPt_genJetPt[i] = new TH2D(Form("h_matchedRecoJetPt_genJetPt_C%i",i),Form("genJetPt vs. matchedRecoJetPt, hiBin %i - %i", centEdges[0]-10,centEdges[NCentralityIndices-1]-10),NPtBins,ptMin,ptMax,NPtBins,ptMin,ptMax);
	    else h_matchedRecoJetPt_genJetPt[i] = new TH2D(Form("h_matchedRecoJetPt_genJetPt_C%i",i),Form("genJetPt vs. matchedRecoJetPt, hiBin %i - %i", centEdges[i-1]-10,centEdges[i]-10),NPtBins,ptMin,ptMax,NPtBins,ptMin,ptMax) ;

        h_matchedRecoJetPt_genJetPt[i]->Sumw2(); 


		if(i==0) h_matchedRawJetPt_genJetPt[i] = new TH2D(Form("h_matchedRawJetPt_genJetPt_C%i",i),Form("genJetPt vs. matchedRawJetPt, hiBin %i - %i", centEdges[0]-10,centEdges[NCentralityIndices-1]-10),NPtBins,ptMin,ptMax,NPtBins,ptMin,ptMax);
	    else h_matchedRawJetPt_genJetPt[i] = new TH2D(Form("h_matchedRawJetPt_genJetPt_C%i",i),Form("genJetPt vs. matchedRawJetPt, hiBin %i - %i", centEdges[i-1]-10,centEdges[i]-10),NPtBins,ptMin,ptMax,NPtBins,ptMin,ptMax) ;

        h_matchedRawJetPt_genJetPt[i]->Sumw2(); 



        if(i==0) h_matchedRecoJetPt_genJetPt_var[i] = new TH2D(Form("h_matchedRecoJetPt_genJetPt_var_C%i",i),Form("genJetPt vs. matchedRecoJetPt, hiBin %i - %i, var bins", centEdges[0]-10,centEdges[NCentralityIndices-1]-10),N1-1,ptAxis1,N2-1,ptAxis2);
	    else h_matchedRecoJetPt_genJetPt_var[i] = new TH2D(Form("h_matchedRecoJetPt_genJetPt_var_C%i",i),Form("genJetPt vs. matchedRecoJetPt, hiBin %i - %i, var bins", centEdges[i-1]-10,centEdges[i]-10),N1-1,ptAxis1,N2-1,ptAxis2) ;

        h_matchedRecoJetPt_genJetPt_var[i]->Sumw2();  

		if(i==0) h_matchedRawJetPt_genJetPt_var[i] = new TH2D(Form("h_matchedRawJetPt_genJetPt_var_C%i",i),Form("genJetPt vs. matchedRawJetPt, hiBin %i - %i, var bins", centEdges[0]-10,centEdges[NCentralityIndices-1]-10),N1-1,ptAxis1,N2-1,ptAxis2);
	    else h_matchedRawJetPt_genJetPt_var[i] = new TH2D(Form("h_matchedRawJetPt_genJetPt_var_C%i",i),Form("genJetPt vs. matchedRawJetPt, hiBin %i - %i, var bins", centEdges[i-1]-10,centEdges[i]-10),N1-1,ptAxis1,N2-1,ptAxis2) ;

        h_matchedRawJetPt_genJetPt_var[i]->Sumw2();   

		if(i==0) {
			h_matchedRecoJetPtOverGenJetPt_genJetPt[i][0] = new TH2D(Form("h_matchedRecoJetPtOverGenJetPt_genJetPt_allJets_C%i",i),Form("matchedRecoJetPt/genJetPt vs genJetPt, all flavors, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),500,0,5,NPtBins,ptMin,ptMax);
			h_matchedRecoJetPtOverGenJetPt_genJetPt[i][1] = new TH2D(Form("h_matchedRecoJetPtOverGenJetPt_genJetPt_bJets_C%i",i),Form("matchedRecoJetPt/genJetPt vs genJetPt, bJets, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),500,0,5,NPtBins,ptMin,ptMax);
			h_matchedRecoJetPtOverGenJetPt_genJetPt[i][2] = new TH2D(Form("h_matchedRecoJetPtOverGenJetPt_genJetPt_cJets_C%i",i),Form("matchedRecoJetPt/genJetPt vs genJetPt, cJets, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),500,0,5,NPtBins,ptMin,ptMax);
			h_matchedRecoJetPtOverGenJetPt_genJetPt[i][3] = new TH2D(Form("h_matchedRecoJetPtOverGenJetPt_genJetPt_udJets_C%i",i),Form("matchedRecoJetPt/genJetPt vs genJetPt, udJets, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),500,0,5,NPtBins,ptMin,ptMax);
			h_matchedRecoJetPtOverGenJetPt_genJetPt[i][4] = new TH2D(Form("h_matchedRecoJetPtOverGenJetPt_genJetPt_sJets_C%i",i),Form("matchedRecoJetPt/genJetPt vs genJetPt, sJets, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),500,0,5,NPtBins,ptMin,ptMax);
			h_matchedRecoJetPtOverGenJetPt_genJetPt[i][5] = new TH2D(Form("h_matchedRecoJetPtOverGenJetPt_genJetPt_gJets_C%i",i),Form("matchedRecoJetPt/genJetPt vs genJetPt, gJets, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),500,0,5,NPtBins,ptMin,ptMax);
			h_matchedRecoJetPtOverGenJetPt_genJetPt[i][6] = new TH2D(Form("h_matchedRecoJetPtOverGenJetPt_genJetPt_xJets_C%i",i),Form("matchedRecoJetPt/genJetPt vs genJetPt, xJets, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),500,0,5,NPtBins,ptMin,ptMax);   
		}
		else{
			h_matchedRecoJetPtOverGenJetPt_genJetPt[i][0] = new TH2D(Form("h_matchedRecoJetPtOverGenJetPt_genJetPt_allJets_C%i",i),Form("matchedRecoJetPt/genJetPt vs genJetPt, all flavors, hiBin %i - %i",centEdges[i-1],centEdges[i]),500,0,5,NPtBins,ptMin,ptMax);
			h_matchedRecoJetPtOverGenJetPt_genJetPt[i][1] = new TH2D(Form("h_matchedRecoJetPtOverGenJetPt_genJetPt_bJets_C%i",i),Form("matchedRecoJetPt/genJetPt vs genJetPt, bJets, hiBin %i - %i",centEdges[i-1],centEdges[i]),500,0,5,NPtBins,ptMin,ptMax);
			h_matchedRecoJetPtOverGenJetPt_genJetPt[i][2] = new TH2D(Form("h_matchedRecoJetPtOverGenJetPt_genJetPt_cJets_C%i",i),Form("matchedRecoJetPt/genJetPt vs genJetPt, cJets, hiBin %i - %i",centEdges[i-1],centEdges[i]),500,0,5,NPtBins,ptMin,ptMax);
			h_matchedRecoJetPtOverGenJetPt_genJetPt[i][3] = new TH2D(Form("h_matchedRecoJetPtOverGenJetPt_genJetPt_udJets_C%i",i),Form("matchedRecoJetPt/genJetPt vs genJetPt, udJets, hiBin %i - %i",centEdges[i-1],centEdges[i]),500,0,5,NPtBins,ptMin,ptMax);
			h_matchedRecoJetPtOverGenJetPt_genJetPt[i][4] = new TH2D(Form("h_matchedRecoJetPtOverGenJetPt_genJetPt_sJets_C%i",i),Form("matchedRecoJetPt/genJetPt vs genJetPt, sJets, hiBin %i - %i",centEdges[i-1],centEdges[i]),500,0,5,NPtBins,ptMin,ptMax);
			h_matchedRecoJetPtOverGenJetPt_genJetPt[i][5] = new TH2D(Form("h_matchedRecoJetPtOverGenJetPt_genJetPt_gJets_C%i",i),Form("matchedRecoJetPt/genJetPt vs genJetPt, gJets, hiBin %i - %i",centEdges[i-1],centEdges[i]),500,0,5,NPtBins,ptMin,ptMax);
			h_matchedRecoJetPtOverGenJetPt_genJetPt[i][6] = new TH2D(Form("h_matchedRecoJetPtOverGenJetPt_genJetPt_xJets_C%i",i),Form("matchedRecoJetPt/genJetPt vs genJetPt, xJets, hiBin %i - %i",centEdges[i-1],centEdges[i]),500,0,5,NPtBins,ptMin,ptMax); 
		}  

		h_matchedRecoJetPtOverGenJetPt_genJetPt[i][0]->Sumw2();
		h_matchedRecoJetPtOverGenJetPt_genJetPt[i][1]->Sumw2();
		h_matchedRecoJetPtOverGenJetPt_genJetPt[i][2]->Sumw2();
		h_matchedRecoJetPtOverGenJetPt_genJetPt[i][3]->Sumw2();
		h_matchedRecoJetPtOverGenJetPt_genJetPt[i][4]->Sumw2();
		h_matchedRecoJetPtOverGenJetPt_genJetPt[i][5]->Sumw2();
		h_matchedRecoJetPtOverGenJetPt_genJetPt[i][6]->Sumw2();

		if(i==0) {
			h_matchedRawJetPtOverGenJetPt_genJetPt[i][0] = new TH2D(Form("h_matchedRawJetPtOverGenJetPt_genJetPt_allJets_C%i",i),Form("matchedRawJetPt/genJetPt vs genJetPt, all flavors, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),500,0,5,NPtBins,ptMin,ptMax);
			h_matchedRawJetPtOverGenJetPt_genJetPt[i][1] = new TH2D(Form("h_matchedRawJetPtOverGenJetPt_genJetPt_bJets_C%i",i),Form("matchedRawJetPt/genJetPt vs genJetPt, bJets, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),500,0,5,NPtBins,ptMin,ptMax);
			h_matchedRawJetPtOverGenJetPt_genJetPt[i][2] = new TH2D(Form("h_matchedRawJetPtOverGenJetPt_genJetPt_cJets_C%i",i),Form("matchedRawJetPt/genJetPt vs genJetPt, cJets, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),500,0,5,NPtBins,ptMin,ptMax);
			h_matchedRawJetPtOverGenJetPt_genJetPt[i][3] = new TH2D(Form("h_matchedRawJetPtOverGenJetPt_genJetPt_udJets_C%i",i),Form("matchedRawJetPt/genJetPt vs genJetPt, udJets, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),500,0,5,NPtBins,ptMin,ptMax);
			h_matchedRawJetPtOverGenJetPt_genJetPt[i][4] = new TH2D(Form("h_matchedRawJetPtOverGenJetPt_genJetPt_sJets_C%i",i),Form("matchedRawJetPt/genJetPt vs genJetPt, sJets, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),500,0,5,NPtBins,ptMin,ptMax);
			h_matchedRawJetPtOverGenJetPt_genJetPt[i][5] = new TH2D(Form("h_matchedRawJetPtOverGenJetPt_genJetPt_gJets_C%i",i),Form("matchedRawJetPt/genJetPt vs genJetPt, gJets, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),500,0,5,NPtBins,ptMin,ptMax);
			h_matchedRawJetPtOverGenJetPt_genJetPt[i][6] = new TH2D(Form("h_matchedRawJetPtOverGenJetPt_genJetPt_xJets_C%i",i),Form("matchedRawJetPt/genJetPt vs genJetPt, xJets, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),500,0,5,NPtBins,ptMin,ptMax);   
		}
		else{
			h_matchedRawJetPtOverGenJetPt_genJetPt[i][0] = new TH2D(Form("h_matchedRawJetPtOverGenJetPt_genJetPt_allJets_C%i",i),Form("matchedRawJetPt/genJetPt vs genJetPt, all flavors, hiBin %i - %i",centEdges[i-1],centEdges[i]),500,0,5,NPtBins,ptMin,ptMax);
			h_matchedRawJetPtOverGenJetPt_genJetPt[i][1] = new TH2D(Form("h_matchedRawJetPtOverGenJetPt_genJetPt_bJets_C%i",i),Form("matchedRawJetPt/genJetPt vs genJetPt, bJets, hiBin %i - %i",centEdges[i-1],centEdges[i]),500,0,5,NPtBins,ptMin,ptMax);
			h_matchedRawJetPtOverGenJetPt_genJetPt[i][2] = new TH2D(Form("h_matchedRawJetPtOverGenJetPt_genJetPt_cJets_C%i",i),Form("matchedRawJetPt/genJetPt vs genJetPt, cJets, hiBin %i - %i",centEdges[i-1],centEdges[i]),500,0,5,NPtBins,ptMin,ptMax);
			h_matchedRawJetPtOverGenJetPt_genJetPt[i][3] = new TH2D(Form("h_matchedRawJetPtOverGenJetPt_genJetPt_udJets_C%i",i),Form("matchedRawJetPt/genJetPt vs genJetPt, udJets, hiBin %i - %i",centEdges[i-1],centEdges[i]),500,0,5,NPtBins,ptMin,ptMax);
			h_matchedRawJetPtOverGenJetPt_genJetPt[i][4] = new TH2D(Form("h_matchedRawJetPtOverGenJetPt_genJetPt_sJets_C%i",i),Form("matchedRawJetPt/genJetPt vs genJetPt, sJets, hiBin %i - %i",centEdges[i-1],centEdges[i]),500,0,5,NPtBins,ptMin,ptMax);
			h_matchedRawJetPtOverGenJetPt_genJetPt[i][5] = new TH2D(Form("h_matchedRawJetPtOverGenJetPt_genJetPt_gJets_C%i",i),Form("matchedRawJetPt/genJetPt vs genJetPt, gJets, hiBin %i - %i",centEdges[i-1],centEdges[i]),500,0,5,NPtBins,ptMin,ptMax);
			h_matchedRawJetPtOverGenJetPt_genJetPt[i][6] = new TH2D(Form("h_matchedRawJetPtOverGenJetPt_genJetPt_xJets_C%i",i),Form("matchedRawJetPt/genJetPt vs genJetPt, xJets, hiBin %i - %i",centEdges[i-1],centEdges[i]),500,0,5,NPtBins,ptMin,ptMax); 
		}  

		h_matchedRawJetPtOverGenJetPt_genJetPt[i][0]->Sumw2();
		h_matchedRawJetPtOverGenJetPt_genJetPt[i][1]->Sumw2();
		h_matchedRawJetPtOverGenJetPt_genJetPt[i][2]->Sumw2();
		h_matchedRawJetPtOverGenJetPt_genJetPt[i][3]->Sumw2();
		h_matchedRawJetPtOverGenJetPt_genJetPt[i][4]->Sumw2();
		h_matchedRawJetPtOverGenJetPt_genJetPt[i][5]->Sumw2();
		h_matchedRawJetPtOverGenJetPt_genJetPt[i][6]->Sumw2();


    }



TFile *f = TFile::Open(input);
cout << "File opened!" << endl;
auto em = new eventMap(f);
em->isMC = 1;
em->AASetup = 0;
cout << "Initializing variables ... " << endl;
em->init();
cout << "Loading jet..." << endl;
em->loadJet("akFlowPuCs4PFJetAnalyzer");
cout << "Loading muon..." << endl;
em->loadMuon("ggHiNtuplizerGED");
cout << "Loading muon triggers..." << endl;
em->loadMuonTrigger("hltanalysis");
cout << "Loading tracks..." << endl;
em->loadTrack();
cout << "Loading gen particles..." << endl;
em->loadGenParticle();
cout << "Variables initilized!" << endl << endl ;

int NEvents = em->evtTree->GetEntries();
cout << "Number of events = " << NEvents << endl;


// define event filters
std::string filters[5] = {"pprimaryVertexFilter", "HBHENoiseFilterResultRun2Loose", "collisionEventSelectionAODv2", "phfCoincFilter3Th4", "pclusterCompatibilityFilter"};
em->regEventFilter(5, filters);

// event loop
int evi_frac = 0;
for(int evi = 0; evi < NEvents ; evi++){
//for(int evi = 0; evi < 10 ; evi++){
	if(evi==0) cout << "Processing events..." << endl;

	em->getEvent(evi);

	if((100*evi / NEvents) % 5 == 0 && 100*evi / NEvents > evi_frac) cout << "evt frac: " << evi_frac << "%" << endl;
	evi_frac = 100 * evi/NEvents;

	// global event cuts
	//cout << "Applying global event cuts..." << endl;
	if(em->pthat <= pthatcut) continue;
	if(em->vz > 15.0) continue;
    if(!passesHiBinCut(em->hiBin)) continue;
	if(em->checkEventFilter()) continue;
	//cout << "Event #" << evi << " passed the global cuts!" << endl;

	// event weights
	double w_hiBin = fCentralityWeightFunction->Eval(em->hiBin-10);		
	double w_vz = fVzWeightFunction->Eval(em->vz);

	double w = em->weight * w_vz * w_hiBin;

	// RECO VARIABLES
	
	int matchFlag[10] = {0,0,0,0,0,0,0,0,0,0};

	int CentralityIndex = getCentBin(em->hiBin);

	if(CentralityIndex < 0) continue;


	// RECO JET LOOP
	for(int i = 0; i < em->ngj ; i++){

		// JET VARIABLES
		double x = em->genjetpt[i];
		double y = em->genjeteta[i];
		double z = em->genjetphi[i];

        double matchedRecoJetPt = 0.0;
		double matchedRawJetPt = 0.0;

		
		if(TMath::Abs(y) > etaMax || x < jetPtCut) continue;

	
        // GET FLAVOR FROM RECO MATCH
		bool hasRecoJetMatch = false;
		double minDr = 100.0;
		int recoJetFlavorFlag = 0;
		int jetFlavorInt = 19;
	

		

		for(int k = 0; k < em->njet; k++){
		
			double dr = getDr(em->jeteta[k],em->jetphi[k],y,z);

			if(dr < minDr){ 

				minDr = dr;

				if(minDr < epsilon_mm){
				
					hasRecoJetMatch = true;
					recoJetFlavorFlag = k;

                    JEC.SetJetPT(em->rawpt[k]);
		            JEC.SetJetEta(em->jeteta[k]);
		            JEC.SetJetPhi(em->jetphi[k]);

		            matchedRecoJetPt = JEC.GetCorrectedPT();
					//matchedRecoJetPt = em->jetpt[k];
					matchedRawJetPt = em->rawpt[k];

					if(x>60) cout << "dR(reco,gen) = " << minDr << " | genPt = " << x << " | rawPt = " << em->rawpt[k] << " | jetPt = " << em->jetpt[k] << " | corrPt = " << JEC.GetCorrectedPT() << endl;
			
				}	
			}

		}

		jetFlavorInt = em->refparton_flavorForB[recoJetFlavorFlag];
     

	


	

			

			

			
			
        // fill response matrix
        if(hasRecoJetMatch) {
			h_matchedRecoJetPt_genJetPt[0]->Fill(matchedRecoJetPt,x,w);
            h_matchedRecoJetPt_genJetPt[CentralityIndex]->Fill(matchedRecoJetPt,x,w);

			h_matchedRawJetPt_genJetPt[0]->Fill(matchedRawJetPt,x,w);
            h_matchedRawJetPt_genJetPt[CentralityIndex]->Fill(matchedRawJetPt,x,w);

			h_matchedRecoJetPt_genJetPt_var[0]->Fill(matchedRecoJetPt,x,w);
            h_matchedRecoJetPt_genJetPt_var[CentralityIndex]->Fill(matchedRecoJetPt,x,w);

			h_matchedRawJetPt_genJetPt_var[0]->Fill(matchedRawJetPt,x,w);
            h_matchedRawJetPt_genJetPt_var[CentralityIndex]->Fill(matchedRawJetPt,x,w);

			h_matchedRecoJetPtOverGenJetPt_genJetPt[0][0]->Fill(matchedRecoJetPt/x,x,w);
			h_matchedRecoJetPtOverGenJetPt_genJetPt[CentralityIndex][0]->Fill(matchedRecoJetPt/x,x,w);

			h_matchedRawJetPtOverGenJetPt_genJetPt[0][0]->Fill(matchedRawJetPt/x,x,w);
			h_matchedRawJetPtOverGenJetPt_genJetPt[CentralityIndex][0]->Fill(matchedRawJetPt/x,x,w);

			if(fabs(jetFlavorInt) == 5){
				h_matchedRecoJetPtOverGenJetPt_genJetPt[0][1]->Fill(matchedRecoJetPt/x,x,w);
				h_matchedRecoJetPtOverGenJetPt_genJetPt[CentralityIndex][1]->Fill(matchedRecoJetPt/x,x,w);

				h_matchedRawJetPtOverGenJetPt_genJetPt[0][1]->Fill(matchedRawJetPt/x,x,w);
				h_matchedRawJetPtOverGenJetPt_genJetPt[CentralityIndex][1]->Fill(matchedRawJetPt/x,x,w);
			} 
			if(fabs(jetFlavorInt) == 4){
				h_matchedRecoJetPtOverGenJetPt_genJetPt[0][2]->Fill(matchedRecoJetPt/x,x,w);
				h_matchedRecoJetPtOverGenJetPt_genJetPt[CentralityIndex][2]->Fill(matchedRecoJetPt/x,x,w);

				h_matchedRawJetPtOverGenJetPt_genJetPt[0][2]->Fill(matchedRawJetPt/x,x,w);
				h_matchedRawJetPtOverGenJetPt_genJetPt[CentralityIndex][2]->Fill(matchedRawJetPt/x,x,w);
			} 
			if(fabs(jetFlavorInt) == 1 || fabs(jetFlavorInt) == 2){
				h_matchedRecoJetPtOverGenJetPt_genJetPt[0][3]->Fill(matchedRecoJetPt/x,x,w);
				h_matchedRecoJetPtOverGenJetPt_genJetPt[CentralityIndex][3]->Fill(matchedRecoJetPt/x,x,w);

				h_matchedRawJetPtOverGenJetPt_genJetPt[0][3]->Fill(matchedRawJetPt/x,x,w);
				h_matchedRawJetPtOverGenJetPt_genJetPt[CentralityIndex][3]->Fill(matchedRawJetPt/x,x,w);
			} 
			if(fabs(jetFlavorInt) == 3){
				h_matchedRecoJetPtOverGenJetPt_genJetPt[0][4]->Fill(matchedRecoJetPt/x,x,w);
				h_matchedRecoJetPtOverGenJetPt_genJetPt[CentralityIndex][4]->Fill(matchedRecoJetPt/x,x,w);

				h_matchedRawJetPtOverGenJetPt_genJetPt[0][4]->Fill(matchedRawJetPt/x,x,w);
				h_matchedRawJetPtOverGenJetPt_genJetPt[CentralityIndex][4]->Fill(matchedRawJetPt/x,x,w);
			}  
			if(jetFlavorInt == 21){
				h_matchedRecoJetPtOverGenJetPt_genJetPt[0][5]->Fill(matchedRecoJetPt/x,x,w);
				h_matchedRecoJetPtOverGenJetPt_genJetPt[CentralityIndex][5]->Fill(matchedRecoJetPt/x,x,w);

				h_matchedRawJetPtOverGenJetPt_genJetPt[0][5]->Fill(matchedRawJetPt/x,x,w);
				h_matchedRawJetPtOverGenJetPt_genJetPt[CentralityIndex][5]->Fill(matchedRawJetPt/x,x,w);
			}  
			if(jetFlavorInt == 0){
				h_matchedRecoJetPtOverGenJetPt_genJetPt[0][6]->Fill(matchedRecoJetPt/x,x,w);
				h_matchedRecoJetPtOverGenJetPt_genJetPt[CentralityIndex][6]->Fill(matchedRecoJetPt/x,x,w);

				h_matchedRawJetPtOverGenJetPt_genJetPt[0][6]->Fill(matchedRawJetPt/x,x,w);
				h_matchedRawJetPtOverGenJetPt_genJetPt[CentralityIndex][6]->Fill(matchedRawJetPt/x,x,w);
			}
		}
			



		

		




	}
	// END GEN JET LOOP

	

} // END EVENT LOOP
delete f;
// WRITE
auto wf = TFile::Open(output,"recreate");

for(int j = 0; j < NCentralityIndices; j++){
    h_matchedRecoJetPt_genJetPt[j]->Write();
    h_matchedRecoJetPt_genJetPt_var[j]->Write();
	h_matchedRecoJetPtOverGenJetPt_genJetPt[j][0]->Write();
	h_matchedRecoJetPtOverGenJetPt_genJetPt[j][1]->Write();
	h_matchedRecoJetPtOverGenJetPt_genJetPt[j][2]->Write();
	h_matchedRecoJetPtOverGenJetPt_genJetPt[j][3]->Write();
	h_matchedRecoJetPtOverGenJetPt_genJetPt[j][4]->Write();
	h_matchedRecoJetPtOverGenJetPt_genJetPt[j][5]->Write();
	h_matchedRecoJetPtOverGenJetPt_genJetPt[j][6]->Write();

	h_matchedRawJetPt_genJetPt[j]->Write();
    h_matchedRawJetPt_genJetPt_var[j]->Write();
	h_matchedRawJetPtOverGenJetPt_genJetPt[j][0]->Write();
	h_matchedRawJetPtOverGenJetPt_genJetPt[j][1]->Write();
	h_matchedRawJetPtOverGenJetPt_genJetPt[j][2]->Write();
	h_matchedRawJetPtOverGenJetPt_genJetPt[j][3]->Write();
	h_matchedRawJetPtOverGenJetPt_genJetPt[j][4]->Write();
	h_matchedRawJetPtOverGenJetPt_genJetPt[j][5]->Write();
	h_matchedRawJetPtOverGenJetPt_genJetPt[j][6]->Write();
}






wf->Close();
return;
// END WRITE



}
