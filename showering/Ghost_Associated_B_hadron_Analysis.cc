#include <math.h>
#include <vector>
#include <string>
#include <sstream>
#include <set>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TMath.h"
#include "TVector3.h"
#include "TRandom3.h"

#include "Pythia8/Pythia.h"
#include "Pythia8/ParticleData.h"

#include "Ghost_Associated_B_hadron_Analysis.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"  
#include "fastjet/tools/Filter.hh"
#include "fastjet/Selector.hh"

#include "fastjet/tools/Recluster.hh"

#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/contrib/DistanceMeasure.hh"
#include "fastjet/contrib/QCDAwarePlugin.hh"


//#include "Pythia8/Pythia.h"

using namespace std;

// Constructor 
Ghost_Associated_B_hadron_Analysis::Ghost_Associated_B_hadron_Analysis(){

	//0.8AK jets
	jetRad = 0.8;
	m_jet_def	= new fastjet::JetDefinition(fastjet::antikt_algorithm, jetRad);
	
	//Constants for use in our algorithm
	//Following paper http://xxx.lanl.gov/pdf/0802.2470v2
	mu = 0.67;
	ycut = 0.15;
    
}

// Destructor 
Ghost_Associated_B_hadron_Analysis::~Ghost_Associated_B_hadron_Analysis(){
	delete m_jet_def;
}

// End
void Ghost_Associated_B_hadron_Analysis::End(){

	cout << "The number that were cut by PTMISS is " << numPTMISS << endl;
	cout << "The number that were cut by PTCUT is " << numPTCUT << endl;
	cout << "The number that were cut by BTAG is " << numBTAG << endl;
	cout << "The number that were cut by ISOLATED is " << numISOLATED << endl;
	cout << "The number that were cut by RHO is " << numRHO << endl;
	cout << "The number that were cut by N2 is " << numN2 << endl;

	return;
}

fastjet::ClusterSequence *largeRClusterSequence; 
std::vector<fastjet::PseudoJet> *particlesForJets; 
fastjet::PseudoJet *p;
std::vector <fastjet::PseudoJet> *largeRJets; 
std::vector <fastjet::PseudoJet> *bs;
std::vector <fastjet::PseudoJet> *pieces;
int idmod; //For tagging b-hadrons
float success = 0;
fastjet::PseudoJet trimmedJet;
vector <fastjet::PseudoJet> constit;
std::list<int> jetIndices;
std::list<int>::iterator it;
int cSize;
int totaBTag;
int thisId;
int Vtype;

int failed = 0;


//Reculster with AK0.2 jets
fastjet::Recluster recluster_ca_inf = fastjet::Recluster(fastjet::antikt_algorithm, 0.2, fastjet::Recluster::keep_all);



double DeltaPhi(double phi1, double phi2){
  if (abs(phi1 - phi2) > M_PI) return 2 * M_PI - abs(phi1 - phi2);
  else return abs(phi1 - phi2);
}

double DeltaR(double phi1, double phi2, double eta1, double eta2){
  return sqrt(pow(DeltaPhi(phi1, phi2), 2) + pow(eta1 - eta2, 2));
}


int digit(const int & loc, const int & pid){
    //  PID digits (base 10) are: n nr nl nq1 nq2 nq3 nj
    //   nj = 1, nq3=2 , nq2=3, nq1, nl, nr, n, n8, n9, n10 
    //  the location enum provides a convenient index into the PID
    int numerator = (int) std::pow(10.0,(loc-1));
    return (abs(pid)/numerator)%10;
}


int fundamentalID(const int & pid){
        if( fabs(pid)/10000000 > 0 ) return 0; //to avoid the error message from Pythia
        if( digit(3,pid) == 0 && digit(4,pid) == 0) {
            return abs(pid)%10000;
        } else if( abs(pid) <= 100 ) {
            return abs(pid);
        } else {
            return 0;
        }
}


//get B hadron
//PID for B hadron are 5XX, 5XXX
//https://gitlab.com/hepcedar/rivet/-/blob/release-3-1-x/analyses/pluginCMS/CMS_2015_I1370682.cc#L390
//https://rivet.hepforge.org/code/2.1.0/a00827.html#ad4c917595339ea52152c2950ce1225e7
bool hasBottom(const int & pid){
        if( fundamentalID(pid) > 0 ) { return false; } //reject "Fundamental particle"
        if( digit(2,pid) == 5 || digit(3,pid) == 5 || digit(4,pid) == 5 ) { return true; } 
        return false;
}


void destroy(){
	delete largeRClusterSequence;
	delete largeRJets;
	delete particlesForJets;
	delete bs;
	delete pieces;
}



// Perform clustering and pre-selection. The first argument is the showered event. 
// The second argument is a target array to write higgs jet information to.
// This target must be a float[4]. It will contain [pt, eta, phi, m] after the method is executed.
void Ghost_Associated_B_hadron_Analysis::AnalyzeEvent(Pythia8::Pythia* pythia8, float *target){

    double Ghostparam = 1e-20;

	if (!pythia8->next()){
		cout << "Failed in Analysis" << endl;
		failed++;
		target[0] = 0; //pT
		target[1] = 0; //eta
		target[2] = 0; //phi
        target[3] = -10; //mass
        target[4] = 0; //VH V type
        target[5] = 0; //DR_Bhadron//temporary
		return;
	}
	else{
	}

	particlesForJets = new std::vector<fastjet::PseudoJet>;
    
    Vtype = 0;
    
    // Particle jet clustering R=0.8 -> Large-R Jet
    // Particle loop
    for (unsigned int ip=0; ip<pythia8->event.size(); ip++){
        p = new fastjet::PseudoJet(pythia8->event[ip].px(), pythia8->event[ip].py(), pythia8->event[ip].pz(),pythia8->event[ip].e());
      
        
        if (pythia8->event[ip].idAbs() == 23 && pythia8->event[pythia8->event[ip].daughter1()].idAbs() != 23 && pythia8->event[pythia8->event[ip].daughter2()].idAbs() != 23){
          Vtype = 23;
            cout << "Vtype : =   " << Vtype << endl;  
            }
        if (pythia8->event[ip].idAbs() == 24 && pythia8->event[pythia8->event[ip].daughter1()].idAbs() != 24 && pythia8->event[pythia8->event[ip].daughter2()].idAbs() != 24){
          Vtype = 24;
          cout << "Vtype : =   " << Vtype << endl;
            }
        
        // Final state only.
      if (!pythia8->event[ip].isFinal()) continue;
      // No neutrinos or DM.
      // 12 for electron neutrino, 14 for muon neutrino, 16 for tauon neutrino, 52 for dark matter with spin 1 / 2
      // Pdgid can be accessed in https://twiki.cern.ch/twiki/bin/view/Main/PdgId
      if ( pythia8->event[ip].idAbs() == 12 || pythia8->event[ip].idAbs() == 14
          || pythia8->event[ip].idAbs() == 16)
          continue;
        
        //here let set_user_index(0) be non-ghost-associated B hadrons
        (*p).set_user_index(0); 
        (*p).set_user_info(new MyUserInfo(pythia8->event[ip].id() , ip , pythia8->event[ip].charge()));
        (*particlesForJets).push_back((*p));     
        delete p;
    } 
    
    
    
    std::vector <fastjet::PseudoJet> forBhadron_dr;//temporary
    //Prepare Ghost-Association B-hadrons and Add ghost-association B-hadrons in to particle jet
    for (int i = 0; i < pythia8->event.size(); i++){
        fastjet::PseudoJet findBhadron(pythia8->event[i].px(), pythia8->event[i].py(), pythia8->event[i].pz(), pythia8->event[i].e());
        fastjet::PseudoJet findBhadron_dr(pythia8->event[i].px(), pythia8->event[i].py(), pythia8->event[i].pz(), pythia8->event[i].e());//temporary
    
        if (!pythia8->event[i].isHadron()) continue;  
        if (!hasBottom(pythia8->event[i].idAbs())) continue;
//         cout << "B Hadron ID: " << pythia8->event[i].idAbs() << endl;
        
        float PT = pow(pow(pythia8->event[i].px(), 2) + pow(pythia8->event[i].py(), 2), 0.5);
        if ( PT < 5) continue;
            
        bool isLast = true;

        for (int daughter: pythia8->event[i].daughterList()){
            if (hasBottom(pythia8->event[daughter].idAbs())){
                isLast = false;
                cout << "not last B Hadron and it's daughter is " << pythia8->event[daughter].id() << endl;
                break;
                }
            }
        if (!isLast) continue;
        //here let set_user_index(1) be non-ghost-associated B hadrons
        findBhadron.set_user_index(1);
        findBhadron.set_user_info(new MyUserInfo(pythia8->event[i].id() , i , pythia8->event[i].charge()));
        
        //make ghost-association B-hadrons 
        findBhadron.reset_momentum(findBhadron.px() * Ghostparam, findBhadron.py() * Ghostparam, findBhadron.pz() * Ghostparam, findBhadron.e() * Ghostparam);
        
        //Add ghost-association B-hadrons in to particle jet
        (*particlesForJets).push_back(findBhadron);
        forBhadron_dr.push_back(findBhadron_dr);//temporary
        }


    
     
     //Create a R = 0.8 fat jet
    largeRClusterSequence =  new fastjet::ClusterSequence (*particlesForJets, *m_jet_def);
    // Collection of R=0.8 jets
    largeRJets = new std::vector<fastjet::PseudoJet>;
    //let minPT = 20 GeV
    *largeRJets =  fastjet::sorted_by_pt((*largeRClusterSequence).inclusive_jets(20.0)); 

    cout << "largeRJets->size():   " << largeRJets->size() << endl;

    pieces = new std::vector<fastjet::PseudoJet>;
    fastjet::contrib::SoftDrop sd(0.0,0.1);
    int leadingCentralJetId = 0;

    float DR; //temporary
    if (forBhadron_dr.size() >= 2){
     std::vector<fastjet::PseudoJet> bhadron_dr;//temporary
            bhadron_dr =  fastjet::sorted_by_pt(forBhadron_dr);//temporary
        //     float DR = DeltaR(bhadron_dr[0].phi(), bhadron_dr[1].phi(), bhadron_dr[0].eta(), bhadron_dr[1].eta());//temporary
            DR = bhadron_dr[0].delta_R(bhadron_dr[1]);//temporary
            cout << "DR  :  " << DR << endl;//temporary
    }
    
    
    for (int iJet = 0; iJet < (*largeRJets).size(); iJet++){
    totaBTag = 0;
        for (fastjet::PseudoJet con: (*largeRJets)[iJet].constituents()){
               if (fabs(con.user_index() == 1)){
                totaBTag ++;
                }
        }
        
        if (totaBTag >= 2 && fabs((*largeRJets)[iJet].eta()) < 2){
            leadingCentralJetId = iJet;
            cout << "Leading central jet chosen with eta=" << (*largeRJets)[iJet].eta() << endl;
        
            
            break; 
            }
        else if(iJet == (*largeRJets).size()-1){
            cout << "No leading central b-tagged jet!" << endl;
            numBTAG++;
            destroy();
            target[0] = 0; //pT
            target[1] = 0; //eta
            target[2] = 0; //phi
            target[3] = 0; //mass
            target[4] = Vtype; //VH V type
            target[5] = 0; //DR_Bhadron//temporary
            return;
            }
        
        }

    trimmedJet = sd((*largeRJets)[leadingCentralJetId]);
    constit = trimmedJet.constituents();
    cSize = constit.size();
    thisId = fabs(constit[0].user_info<MyUserInfo>().pdg_id());


	// pt cut - For now on the base jet
    if (trimmedJet.pt() < 400){
		numPTCUT++;
		destroy();
		target[0] = 0; //pT
		target[1] = 0; //eta
		target[2] = 0; //phi
		target[3] = 0; //mass
        target[4] = Vtype; //VH V type
        target[5] = 0; //DR_Bhadron//temporary
		return;
	}	


	//Here, we want to check over all the jets
	//And see which are just a single particle, and check if they are isolated or not

	if (cSize == 1){
	if (thisId == 11){
		//we have an ELECTRON 
		if (constit[0].pt() > 10 && fabs(constit[0].eta()) < 2.5){ 
			//we have an isolated electron/muon 
			//cout << "Particle Level : Isolated Electron" << endl;
			numISOLATED++;
			destroy();
			target[0] = 0; //pT
			target[1] = 0; //eta
			target[2] = 0; //phi
			target[3] = 0; //mass
            target[4] = Vtype; //VH V type
            target[5] = 0; //DR_Bhadron//temporary
			return;
		}
	}
	if (thisId == 13){
		//we have an muon
		if (constit[0].pt() > 10 && fabs(constit[0].eta()) < 2.4){
			//we have an isolated electron/muon
			//cout << "Particle Level : Isolated Muon" << endl;
			numISOLATED++;
			destroy();
			target[0] = 0; //pT
			target[1] = 0; //eta
			target[2] = 0; //phi
			target[3] = 0; //mass
            target[4] = Vtype; //VH V type
            target[5] = 0; //DR_Bhadron//temporary
			return;
		}
	}
	if (thisId ==15){
		//we have a tau
		if (constit[0].pt() > 18 && fabs(constit[0].eta()) < 2.3){
			//we have an isolated tau
			//cout << "Particle Level : Isolated Tau" << endl;
			numISOLATED++;
			destroy();
			target[0] = 0; //pT
			target[1] = 0; //eta
			target[2] = 0; //phi
			target[3] = 0; //mass
            target[4] = Vtype; //VH V type
            target[5] = 0; //DR_Bhadron//temporary
			return;
		}
	}
	}

	//Now we do a cut on rho, but we need to get the softdrop working first

	rho = log(pow(trimmedJet.m(),2)/pow((*largeRJets)[leadingCentralJetId].pt(),2)); //Use the ungroomed pt
	if (rho <= -6){
		numRHO++;
		destroy();
		target[0] = 0; //pT
		target[1] = 0; //eta
		target[2] = 0; //phi
		target[3] = 0; //mass
        target[4] = Vtype; //VH V type
        target[5] = 0; //DR_Bhadron//temporary
		return;
	}
	if (rho >= -2.1){
		numRHO++;
		destroy();
		target[0] = 0; //pT
		target[1] = 0; //eta
		target[2] = 0; //phi
		target[3] = 0; //mass
        target[4] = Vtype; //VH V type
        target[5] = 0; //DR_Bhadron//temporary
		return;
	}

	//Now, lets actually just implement an N12 and see what kind of shenanigans we get
	//Here, beta is set to zero
	//https://arxiv.org/pdf/1609.07483.pdf
	
	float e2 = 0;
	float e3 = 0;
	float sumpt = 0;

	float rs = 0;
	float rt = 0;
	float st = 0;

	for (int r = 0; r < cSize; r++){
		sumpt = sumpt + constit[r].pt();
	for (int s = 0; s < cSize; s++){
		if (r == s) continue;	
		e2 = e2 + constit[r].pt()*constit[s].pt()*constit[r].delta_R(constit[s]);
	for (int t = 0; t < cSize; t++){
		if (s == t) continue;
		rs = constit[r].delta_R(constit[s]);
		rt = constit[r].delta_R(constit[t]);
		st = constit[s].delta_R(constit[t]);
		e3 = e3 + constit[r].pt()*constit[s].pt()*constit[t].pt()*std::min({rs*rt,rs*st,rt*st});
	}
	}
	}

	e2 = e2/(pow(sumpt,2));
	e3 = e3/(pow(sumpt,3));

	float N2 = e3/pow(e2,2);
	
	//N2 cut
	if (N2 > 0.45){
		numN2++;
		destroy();
		target[0] = 0; //pT
		target[1] = 0; //eta
		target[2] = 0; //phi
		target[3] = 0; //mass
        target[4] = Vtype; //VH V type
        target[5] = 0; //DR_Bhadron//temporary
		return;
	}
	
	/*
	delete largeRClusterSequence;
	delete largeRJets;
	delete particlesForJets;
	*/
	destroy();	

	//storing higgs jet information to target array
	target[0] = pow(pow(trimmedJet.px(), 2) + pow(trimmedJet.py(), 2), 0.5); //pT
	target[1] = trimmedJet.eta(); //eta
	target[2] = trimmedJet.phi(); //phi
	target[3] = trimmedJet.m(); //mass
    target[4] = Vtype; //VH V type
    target[5] = DR; //DR_Bhadron//temporary
	cout << "higgs output test - analysis ";
	cout << target[0] << "," << target[1] << ",";
	cout << target[2] << "," << target[3] << endl;
//     cout << target[4] << endl;
	cout << target << endl;
}
