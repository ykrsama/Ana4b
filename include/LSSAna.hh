#ifndef _LSSAna_hh_
#define _LSSAna_hh_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>

#include "marlin/Processor.h"
#include "LCIOSTLTypes.h"
#include "EVENT/LCCollection.h"
#include "EVENT/LCParameters.h"
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/MCParticle.h"
#include "UTIL/PIDHandler.h"

#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"


class TTree;


class LSSAna : public marlin::Processor {
public:
    Processor* newProcessor() { return new LSSAna; };

    LSSAna();

    ~LSSAna();

    void init();

    void processEvent( LCEvent* evtP );

    void end();

private:

    double getInvMass(MCParticle* part1, MCParticle* part2);
    double getInvMass(ReconstructedParticle* part1, ReconstructedParticle* part2);

    TLorentzVector getTLorentzVector(MCParticle* part);
    TLorentzVector getTLorentzVector(ReconstructedParticle* part);

    double getMassjj( ReconstructedParticle* part1, ReconstructedParticle* part2 );

    //double getDeltaM( int j1, int j2, int j3, int j4) { return fabs( Mjj[j1][j2] - Mjj[j3][j4] ); };

    //double getRm( int j1, int j2, int j3, int j4 ) { return fabs( ( Mjj[j1][j2] - Mjj[j3][j4] ) / ( Mjj[j1][j2] + Mjj[j3][j4] ) ); };

    bool lessDeltaRjl(ReconstructedParticle *part1, ReconstructedParticle *part2);

    bool lessDeltaRjl(MCParticle *part1, ReconstructedParticle *part2);

    bool greaterPT(ReconstructedParticle *part1, ReconstructedParticle *part2);

    bool greaterE(ReconstructedParticle *part1, ReconstructedParticle *part2);

    bool greaterE(MCParticle* part1, MCParticle *part2);

    std::vector<ReconstructedParticle*> sortLessDeltaRjl(std::vector<ReconstructedParticle*> &ivec);
    
    std::vector<ReconstructedParticle*> sortGreaterPT(std::vector<ReconstructedParticle*> &ivec);

    std::vector<ReconstructedParticle*> sortGreaterE(std::vector<ReconstructedParticle*> &ivec);

    std::vector<MCParticle*> sortGreaterE(std::vector<MCParticle*> &ivec);

    std::vector<ReconstructedParticle*> arrangeJets(std::vector<ReconstructedParticle*> &ivec);

    struct {
        bool operator()(const DoubleVec lcfi1, const DoubleVec lcfi2) {
            return ( lcfi1.at(0) > lcfi2.at(0) );
        };
    } greaterBTag;

private:
    // ROOT related
    TFile* tree_file;
    TTree* _outputTree_truth;
    TTree* _outputTree_event;

    // Processor Parameters
    std::string _treeFileName;
    //std::string _treeName_truth;
    std::string _treeName_event;
    std::string _colName;
    std::ostream _output;
    int _overwrite;
    int fNJet;

    // Vars fill to tree
    //===============================
    // event
    //-------------------------------
    int feventNum;
    double fEvisible;
    double fMass_visible;
    double fEmiss;
    
    // particle
    double feta;
    double fcosThetaj;
    double fyij;
    double fmll; // inv mass of Z leptons
    double fEll; // Energy of Z leptons
    double fmz;
    double fmh;
    double fmrecoil;
    double fdelta_mll_mz;
    double fdelta_mrecoil_mh;
    int fMCbNum; // b quark number
    int fMCcNum;
    int fMCqNum;
    DoubleVec fEb;

    // reco jet
    double fjet0_pT;
    double fjet0_E;
    double fjet1_pT;
    double fjet1_E;
    double fjet2_pT;
    double fjet2_E;
    double fjet3_pT;
    double fjet3_E;
    DoubleVec fPTjet;
    DoubleVec fEjet;
    DoubleVec fDeltaRBJ;
    DoubleVec fEBJ;
    double fbTagNum;
    DoubleVec flcfiplus0;
    DoubleVec flcfiplus1;
    DoubleVec flcfiplus2;
    DoubleVec flcfiplus3;

    bool ffound_Zlepton;

    // jet

    //===============================
    // h1
    double Mj00j01; // invariant mass of h1
    double Mj10j11;
    double fdeltaMjj;
    double fDeltaRj00j01;
    double fDeltaRj10j11;
    StringVec fJetTags;

    //==============================
    // rech1 tree

    double fh1InvMass;
    double fDletaRjj;
    //===============================
    // constant
    int Leptons[6];
    int S;

    // run time vars
    int NMCP; // element number of MCParticle
    MCParticle* Zleptons[2]; // two leptons from Z

    int alcfiplus; // Algorithm ID of lcfiplus
    int ayth; // Algorithm ID of yth
    int ibtag; // Parameter index of BTag
    int ictag; // Parameter index of CTag
    int iltag; // value of light jet param
    int icategory; // Parameter index of Category
    int NJetsNum; // element number of Jets
    double fbTag;
    double fcTag;
    double fCategory;
    double flTag;
    double ftempTagParam;
    DoubleVec lcfiplus_temp;
    std::string ftempTag;
    ReconstructedParticle* leptonJets[2];
    ReconstructedParticle* realJets[4];
    MCParticle* MCbs[4];
    std::vector<DoubleVec> vlcfiplus;

    TLorentzVector Vl[2];
    TLorentzVector Vj[2];
    
};

#endif
