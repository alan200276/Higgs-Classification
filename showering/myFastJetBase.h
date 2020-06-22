#ifndef MYFASTJETBASE_H
#define MYFASTJETBASE_H

#include "fastjet/PseudoJet.hh"


class MyUserInfo : public fastjet::PseudoJet::UserInfoBase{
	public:
		MyUserInfo(const int & pdg_id_in,const int & pythia_id_in,  const double & charge_in) :
			_pdg_id(pdg_id_in),_pythia_id(pythia_id_in), _charge(charge_in){}
		int pdg_id() const { return _pdg_id;}
		int pythia_id() const {return _pythia_id;}
		double charge() const { return _charge;}
	protected:
		int _pdg_id;     // the associated pdg id
		int _pythia_id;  // index in pythia.event
		double _charge;  // the particle charge
};

#endif
