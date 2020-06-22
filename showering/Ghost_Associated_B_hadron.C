// main89.cc is a part of the PYTHIA event generator.
// Copyright (C) 2017 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This program is written by Stefan Prestel.
// It illustrates how to do run PYTHIA with LHEF input, allowing a
// sample-by-sample generation of
// a) Non-matched/non-merged events
// b) MLM jet-matched events (kT-MLM, shower-kT, FxFx)
// c) CKKW-L and UMEPS-merged events
// d) UNLOPS NLO merged events
// see the respective sections in the online manual for details.

#include "Ghost_Associated_B_hadron_Analysis.h"

	
#include "Pythia8/Pythia.h"
#include "Pythia8/ParticleData.h"

#include <unistd.h>

// Include UserHooks for Jet Matching.
#include "Pythia8Plugins/CombineMatchingInput.h"
 // Include UserHooks for randomly choosing between integrated and
// non-integrated treatment for unitarised merging.
#include "Pythia8Plugins/aMCatNLOHooks.h"

using namespace Pythia8;

//==========================================================================

// Example main programm to illustrate merging.
//
//
//

//Random Seed Generator using times
int getSeed(int seed) {
	if (seed > -1) return seed;
	int timeSeed = time(NULL);
	return abs(((timeSeed * 181) * ((getpid() - 83) * 359)) % 104729);
}

int main(int argc, char* argv[]) {

	// Check that correct number of command-line arguments
	if (argc != 5) {
		cerr << " Unexpected number of command-line arguments (" << argc << "). \n" <<
			" You are expected to provide the arguments" << endl <<
			" 1. Input file for settings" << endl <<
			" 2. Output file (csv)" << endl <<
			" 3. Location of LHE file" << endl <<
            " 4. Output Truth Record(bb~) file (csv)" << endl <<
			" Program stopped. " << endl;
		return 1;
	}

	Pythia pythia;
	//Generate our seed
	int seed = -1;
	seed = getSeed(seed);
	pythia.readString("Random:setSeed = on");
	std::stringstream ss;
	ss << "Random:seed = " << seed;
	cout << ss.str() << endl;
	pythia.readString(ss.str());

	ofstream myfile;
	std::string str_name(argv[2]);
	//  std::string str_seed(seed);
	myfile.open(str_name + "_seed_" + std::to_string(seed) + ".csv");
    
    
    
    //Create a File to restore b b~ Information
    ofstream truthrecord;
    std::string str_name_2(argv[4]);
    truthrecord.open(str_name_2 + "_seed_" + std::to_string(seed) + ".csv");
    //ofstream truthrecord_2;//temporary
    //truthrecord_2.open("/home/alan/FrankShower/B_hadron_Test_2/fordeltaR.csv");//temporary
    
    truthrecord << "b_E";
    truthrecord << ",";
    truthrecord << "b_Px";
    truthrecord << ",";
    truthrecord << "b_Py";
    truthrecord << ",";
    truthrecord << "b_Pz";
    truthrecord << ",";
    truthrecord << "bbar_E";
    truthrecord << ",";
    truthrecord << "bbar_Px";
    truthrecord << ",";
    truthrecord << "bbar_Py";
    truthrecord << ",";
    truthrecord << "bbar_Pz";
    truthrecord << ",";
    truthrecord << "evtweight";
    truthrecord << "\n";

    //truthrecord_2 << "dbb" ;//temporary
    //truthrecord_2 << "\n";//temporary


	// New setting to allow processing of multiple input LHEFs.
	pythia.settings.addMode("LHEFInputs:nSubruns", 0, true, false, 0, 100);

	// Input parameters:
	pythia.readFile(argv[1], 0);
	const char *t = argv[3];
	std::string str(t);
	pythia.readString("Beams:LHEF = " + str);

	// Interface for conversion from Pythia8::Event to HepMC one.
	//HepMC::Pythia8ToHepMC ToHepMC;
	// Specify file where HepMC events will be stored.
	//HepMC::IO_GenEvent ascii_io(argv[2], std::ios::out);
	// Switch off warnings for parton-level events.
	//ToHepMC.set_print_inconsistency(false);
	//ToHepMC.set_free_parton_exception(false);
	// Do not store cross section information, as this will be done manually.
	//ToHepMC.set_store_pdf(false);
	//ToHepMC.set_store_proc(false);
	//ToHepMC.set_store_xsec(false);

	// Check if jet matching should be applied.
	bool doMatch = pythia.settings.flag("JetMatching:merge");

	// Check if internal merging should be applied.
	bool doMerge = !(pythia.settings.word("Merging:Process").compare("void") ==
		0);

	cout << "doMatch : " << doMatch << "doMerge : " << doMerge << endl;

	// Currently, only one scheme at a time is allowed.
	if (doMatch && doMerge) {
		cerr << " Jet matching and merging cannot be used simultaneously.\n" <<
			" Program stopped.";
	}

	// Get number of subruns.
	int nMerge = pythia.mode("LHEFInputs:nSubruns");

	// Number of events. Negative numbers mean all events in the LHEF will be
	// used.
	long nEvent = pythia.settings.mode("Main:numberOfEvents");
	if (nEvent < 1) nEvent = 1000000000000000;

	// For jet matching, initialise the respective user hooks code.
	UserHooks* matching = NULL;

	// Allow to set the number of addtional partons dynamically.
	amcnlo_unitarised_interface* setting = NULL;
	if (doMerge) {
		// Store merging scheme.
		int scheme = (pythia.settings.flag("Merging:doUMEPSTree") ||
				pythia.settings.flag("Merging:doUMEPSSubt")) ?
				1 :
					((pythia.settings.flag("Merging:doUNLOPSTree") ||
					pythia.settings.flag("Merging:doUNLOPSSubt") ||
					pythia.settings.flag("Merging:doUNLOPSLoop") ||
					pythia.settings.flag("Merging:doUNLOPSSubtNLO")) ?
				2 :
				0);
		setting = new amcnlo_unitarised_interface(scheme);
		pythia.setUserHooksPtr(setting);
	}

	// For jet matching, initialise the respective user hooks code.
	if (doMatch) {
		CombineMatchingInput combined;
		matching = combined.getHook(pythia);
		if (!matching) {
			cerr << " Failed to initialise jet matching structures.\n" <<
				" Program stopped.";
			return 1;
		}
		pythia.setUserHooksPtr(matching);
	}

	// Cross section and error.
	double sigmaTotal = 0.;
	double errorTotal = 0.;

	// Allow abort of run if many errors.
	int nAbort = pythia.mode("Main:timesAllowErrors");
	int iAbort = 0;
	bool doAbort = false;

	Ghost_Associated_B_hadron_Analysis *analysis1 = new Ghost_Associated_B_hadron_Analysis();

	cout << endl << endl << endl;
	cout << "Start generating events" << endl;

	int numShowered = 0;
	int numSelected = 0; // total number of events stored
	int numTotal = 0; // total number of events analyzed
	int numFiltered = 0; // total number of events filtered by pre-selection
	int numFailed = 0; // total number of events that failed analysis
    
    int total = 0;
        
	// Loop over subruns with varying number of jets.
	for (int iMerge = 0; iMerge < nMerge; ++iMerge) {

		double sigmaSample = 0., errorSample = 0.;

		// Read in name of LHE file for current subrun and initialize.
		pythia.readFile(argv[1], iMerge);
		// Initialize.
		pythia.init();

		// Get the inclusive x-section by summing over all process x-sections.
		double xs = 0.;
		for (int i = 0; i < pythia.info.nProcessesLHEF(); ++i)
			xs += pythia.info.sigmaLHEF(i);

		int iev = 0;
		float trimmedJetMass = 0;
		float higgsJet[6]; // The target array. It will contain [pt, eta, phi, m] after the method is executed.
						   // It will be [0, 0, 0, -10] if there is an error
						   // It will be [0, 0, 0, 0] if an event is filtered

                           
		// Start generation loop
		// This will result in a .csv file with format:
		// particle pT, eta, phi, m, id, isCharged \n  \
		// particle pT, eta, phi, m, id, isCharged \n  |
		// ...                                         > an event
		// particle pT, eta, phi, m, id, isCharged \n  |
		// higgs jet pT, eta, phi, m.\n                /
		// \n                                            an empty line
		// particle pT, eta, phi, m, id, isCharged \n  \
		// particle pT, eta, phi, m, id, isCharged \n  |
		// ...                                         > an event
		// particle pT, eta, phi, m, id, isCharged \n  |
		// higgs jet pT, eta, phi, m, evtweight \n     /
		// \n                                            an empty line
		// ...
        
        
		while (pythia.info.nSelected() < nEvent) {

			numTotal++;
            
			// Get event weight(s).
			double evtweight = pythia.info.weight();
			// Additional PDF/alphaS weight for internal merging.

			//we are not merging so this is tame
			if (doMerge) evtweight *= pythia.info.mergingWeightNLO()
				// Additional weight due to random choice of reclustered/non-reclustered
				// treatment. Also contains additional sign for subtractive samples.
				*setting -> getNormFactor();

			// Do not print zero-weight events.
			// if (evtweight == 0.) continue;
            
            
            
            
            
            //Load Truth b b~ information 
            int nB = 0;
            int pB = 0;
            int nH = 0;
            cout << "event size " << pythia.event.size() << endl;
            cout << "numTotal " << numTotal << endl;
            
            if (pythia.event.size() > 0){
            
            
            for (int h_index = 0; h_index < pythia.event.size(); ++h_index) {
                 if (pythia.event[h_index].id() == 25 ) nH = h_index;
            }
            
            if (pythia.event[pythia.event[nH].daughter1()].id() == 5) {
            nB = pythia.event[nH].daughter1();
            pB = pythia.event[nH].daughter2();
            }
            
            if (pythia.event[pythia.event[nH].daughter1()].id() == -5) {
            pB = pythia.event[nH].daughter1();
            nB = pythia.event[nH].daughter2();
            }
            
            total = total + 1;
            cout << "# of event: " << total << endl;
            cout << "i = " << nB << ", id = " 
                  << pythia.event[nB].id() << endl;
            cout << pythia.event[nB].e() << endl;
            cout << pythia.event[nB].px() << endl;
            cout << pythia.event[nB].py() << endl;
            cout << pythia.event[nB].pz() << endl;   

            cout << "i = " << pB << ", id = " 
                  << pythia.event[pB].id() << endl; 
            cout << pythia.event[pB].e() << endl;
            cout << pythia.event[pB].px() << endl;
            cout << pythia.event[pB].py() << endl;
            cout << pythia.event[pB].pz() << endl; 
            
            
            truthrecord << pythia.event[nB].e();
            truthrecord << ",";
            truthrecord << pythia.event[nB].px();
            truthrecord << ",";
            truthrecord << pythia.event[nB].py();
            truthrecord << ",";
            truthrecord << pythia.event[nB].pz();
            truthrecord << ",";
            truthrecord << pythia.event[pB].e();
            truthrecord << ",";
            truthrecord << pythia.event[pB].px();
            truthrecord << ",";
            truthrecord << pythia.event[pB].py();
            truthrecord << ",";
            truthrecord << pythia.event[pB].pz();
            truthrecord << ",";
            truthrecord << evtweight;
            truthrecord << "\n";
            }



			higgsJet[0] = 0;
			higgsJet[1] = 0;
			higgsJet[2] = 0;
			higgsJet[3] = 0;
            higgsJet[4] = 0;
            higgsJet[5] = 0;
            
            
			analysis1 -> AnalyzeEvent(&pythia, &higgsJet[0]);			
			cout << "higgs output test ----- main ";
			cout << higgsJet[0] << "," << higgsJet[1] << ",";
			cout << higgsJet[2] << "," << higgsJet[3] << endl;	
            cout << "V type: " << higgsJet[4] << endl;	
            cout << "DR(two leading B hadron): " << higgsJet[5] << endl;	
			cout << "Event weight: " << evtweight << endl;
			trimmedJetMass = higgsJet[3];
			if (trimmedJetMass < -9) {
				numFailed++;
				if (pythia.info.atEndOfFile()) break;
				if (pythia.info.nSelected() == 0){
					cout << "End of file!" << endl;
					break; // Break if file is empty
				} 
				cout << "Failed to perform parton shower, error code:" << endl;
				cout << pythia.info.nSelected() << endl;
				continue;
			}
			if (trimmedJetMass == 0) {
				numFiltered++;
			}
			if (trimmedJetMass > 50) {
				//cout << "writing" << endl;
				for (unsigned int ip = 0; ip < pythia.event.size(); ip++) {
					if (!pythia.event[ip].isFinal()) continue;
					if (fabs(pythia.event[ip].id()) == 12) continue;
					if (fabs(pythia.event[ip].id()) == 14) continue;
					if (fabs(pythia.event[ip].id()) == 16) continue;

					myfile << pythia.event[ip].pT();
					myfile << ",";
					myfile << pythia.event[ip].eta();
					myfile << ",";
					myfile << pythia.event[ip].phi();
					myfile << ",";
					myfile << pythia.event[ip].m();
					myfile << ",";
					myfile << pythia.event[ip].id();
					myfile << ",";
					if (pythia.event[ip].isCharged()) {
						myfile << "1";
					} else {
						myfile << "0";
					}
                    //myfile << "" <<","; //temporary
					myfile << "\n";

				}
				myfile << higgsJet[0] << ",";
				myfile << higgsJet[1] << ",";
				myfile << higgsJet[2] << ",";
				myfile << trimmedJetMass << ",";
				myfile << evtweight << ",";
                myfile << higgsJet[4] << "\n"; 
                //myfile << higgsJet[5] << "\n"; //temporary

				myfile << "\n";
                
                //truthrecord_2 << higgsJet[5];//temporary
                //truthrecord_2 << "\n";//temporary
                
				numSelected++;
			}

			if (iev % 200 == 0) {
				cout << "Running on " << iev << " rioght now" << endl;
			}


			// Work with weighted (LHA strategy=-4) events.
			// What is this used for?
			double normhepmc = 1.;
			if (abs(pythia.info.lhaStrategy()) == 4)
				normhepmc = 1. / double(1e9 * nEvent);
			// Work with unweighted events.
			else
				normhepmc = xs / double(1e9 * nEvent);


			// Construct new empty HepMC event.
			//HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
			// Set event weight
			//hepmcevt->weights().push_back(evtweight*normhepmc);
			// Fill HepMC event
			//ToHepMC.fill_next_event( pythia, hepmcevt );
			// Add the weight of the current event to the cross section.
			sigmaTotal += evtweight * normhepmc;
			sigmaSample += evtweight * normhepmc;
			errorTotal += pow2(evtweight * normhepmc);
			errorSample += pow2(evtweight * normhepmc);

			numShowered++;

			// Report cross section to hepmc
			//HepMC::GenCrossSection xsec;
			//xsec.set_cross_section( sigmaTotal*1e9, pythia.info.sigmaErr()*1e9 );
			//hepmcevt->set_cross_section( xsec );
			// Write the HepMC event to file. Done with it.
			//ascii_io << hepmcevt;
			//delete hepmcevt;

			iev++;

		} // end loop over events to generate.
		if (doAbort) break;

		cout << endl;
		cout << "The number that passed all cuts is " << numSelected << " out of " << numTotal-1 << endl; // numtotal counts end of file, and need to be deducted by 1
		cout << "The number that are filtered by pre-selection is " << numFiltered << endl;
		cout << "The number with failed analysis is " << numFailed << endl;
		pythia.settings.list("seed");
		//cout << pythia.settings.list("seed") << " seed";

		analysis1 -> End();
		delete analysis1;

		// print cross section, errors
		pythia.stat();

		cout << endl << " Contribution of sample " << iMerge <<
			" to the inclusive cross section : " <<
			scientific << setprecision(8) <<
			sigmaSample << "  +-  " << sqrt(errorSample) << endl;

	}

	cout << endl << endl << endl;
	if (doAbort)
		cout << " Run was not completed owing to too many aborted events" << endl;
	else
		cout << "Inclusive cross section: " << scientific << setprecision(8) <<
		sigmaTotal << "  +-  " << sqrt(errorTotal) << " mb " << endl;
	cout << endl << endl << endl;

	// Clean-up
	if (doMerge) delete setting;
	if (doMatch) delete matching;

	cout << numSelected << " events stored in final .csv file." << endl;
    
    cout << total << " events stored in truthrecord .csv file." << endl;

	// Done
	return 0;

}
