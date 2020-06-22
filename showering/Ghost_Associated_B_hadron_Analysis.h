#ifndef  Ghost_Associated_B_hadron_Analysis_H
#define  Ghost_Associated_B_hadron_Analysis_H

#include <vector>
#include <math.h>
#include <string>

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"  
#include "fastjet/tools/Filter.hh"
#include "fastjet/Selector.hh"

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TParticle.h"

#include "TH1F.h"

#include "myFastJetBase.h"
#include "Pythia8/Pythia.h"

using namespace std;

class Ghost_Associated_B_hadron_Analysis{
    private:

        string fOutName;
	float mu;
	float ycut;
	float jetRad = 1.2;
	float rho = 0;


	//We have a series of ints that give information on levels of cuts
	//This cuts are listed sequentially, so the top cut is applied first
	int numPTMISS = 0;
	int numPTCUT = 0;
	int numBTAG = 0; //This already kills like... a lot.
	int numISOLATED = 0;
	int numRHO = 0;
	int numN2 = 0;

        TFile *tF;
        TTree *tT;

	int              fTEventNumber;
        static const int MaxNJet = 2;
        int              fTNJets;

	int 		fTJetChargedMultiplicity [MaxNJet];
	int 		fTJetMultiplicity [MaxNJet];
	float 		fTJetEta[MaxNJet];
	float 		fTJetPhi[MaxNJet];
	float 		fTJetPt[MaxNJet];
	float 		fTJetE[MaxNJet];
	float 		fTHiggsPt;
	float 		fTHiggsPhi;
	float 		fTHiggsEta;
	int 		fTHiggsChargedMultiplicity;
	
	int fTmode;

        fastjet::JetDefinition     *m_jet_def;
        fastjet::JetDefinition     *m_jet_def_largeR_ALTAS;
    

    public:

        Ghost_Associated_B_hadron_Analysis ();
        ~Ghost_Associated_B_hadron_Analysis ();
        
        void Begin();
        void AnalyzeEvent(Pythia8::Pythia* pythia8, float *target);
        void End();
        void DeclareBranches();
        void ResetBranches();
        void SetOutName(string outname){
            fOutName = outname;
        }
};

#endif
