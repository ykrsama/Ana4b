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

    double getMCInvMass(MCParticle* part1, MCParticle* part2);
    double getReInvMass(ReconstructedParticle* part1, ReconstructedParticle* part2);

    TLorentzVector getMCTLorentzVector(MCParticle* part);
    TLorentzVector getReTLorentzVector(ReconstructedParticle* part);

    double getMassjj(  double j1E,
                    double j1Px,
                    double j1Py,
                    double j1Pz,
                    double j2E,
                    double j2Px,
                    double j2Py,
                    double j2Pz );


    double getDeltaM( int j1, int j2, int j3, int j4) { return fabs( Mjj[j1][j2] - Mjj[j3][j4] ); };

    double getRm( int j1, int j2, int j3, int j4 ) { return fabs( ( Mjj[j1][j2] - Mjj[j3][j4] ) / ( Mjj[j1][j2] + Mjj[j3][j4] ) ); };

    bool lessDeltaRjl(ReconstructedParticle *part1, ReconstructedParticle *part2);

    bool greaterPT(ReconstructedParticle *part1, ReconstructedParticle *part2);

    std::vector<ReconstructedParticle*> sortLessDeltaRjl(std::vector<ReconstructedParticle*> &ivec);
    
    std::vector<ReconstructedParticle*> sortGreaterPT(std::vector<ReconstructedParticle*> &ivec);

private:
    // ROOT related
    TFile* tree_file;
    TTree* _outputTree_event;
    TTree* _outputTree_rech1;

    // Processor Parameters
    std::string _treeFileName;
    std::string _treeName_event;
    std::string _treeName_rech1;
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
    double fEmiss;
    
    // particle
    double feta;
    double fcosThetaj;
    double fyij;
    double fmll;
    double fmz;
    double fmh;
    double fmrecoil;
    double fdelta_mll_mz;
    double fdelta_mrecoil_mh;
    double fjet0_pT;
    double fjet0_E;
    double fjet1_pT;
    double fjet1_E;
    double fjet2_pT;
    double fjet2_E;
    double fjet3_pT;
    double fjet3_E;

    bool ffound_Zlepton;

    // jet
    DoubleVec fvjTag;

    
    //===============================
    // h1
    double fh1InvMass;
    double fDletaRjj;
    std::string fj1Tag;
    std::string fj2Tag;
    double Mjj[4][4]; // matrix of massjj. Is a upper triangular matrix.

    // constant
    int Leptons[6];

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
    ReconstructedParticle* leptonJets[2];
    ReconstructedParticle* realJets[4];

    TLorentzVector Vl[2];
    TLorentzVector Vj[2];
    
};

#endif
