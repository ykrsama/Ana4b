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
    double getMassjj(  double j1E,
                    double j1Px,
                    double j1Py,
                    double j1Pz,
                    double j2E,
                    double j2Px,
                    double j2Py,
                    double j2Pz );


    //double getDeltaM( int j1, int j2, int j3, int j4) { return fabs( Mjj[j1][j2] - Mjj[j3][j4] ); };

    //double getRm( int j1, int j2, int j3, int j4 ) { return fabs( ( Mjj[j1][j2] - Mjj[j3][j4] ) / ( Mjj[j1][j2] + Mjj[j3][j4] ) ); };

private:
    // ROOT related
    TFile* ptree_file;
    TTree* pevent_tree;
    TTree* prech1_tree;

    // Processor Parameters
    std::string ftree_file_name;
    std::string fevent_tree_name;
    std::string frech1_tree_name;
    std::string fcol_name;
    std::ostream foutput;
    int foverwrite;
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
    double fleading_jet_pT;

    // jet
    DoubleVec fvjTag;

    
    //===============================
    // h1
    double fh1InvMass;
    double fDletaRjj;
    std::string fj1Tag;
    std::string fj2Tag;

    
    // run time vars

    double h1InvMass;
};

#endif
