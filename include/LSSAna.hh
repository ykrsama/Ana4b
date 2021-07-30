#ifndef _LSSAna_hh 
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


class LSSAna : public marlin::Processor
{
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


    double getDeltaM( int j1, int j2, int j3, int j4) { return fabs( Mjj[j1][j2] - Mjj[j3][j4] ); };

    double getRm( int j1, int j2, int j3, int j4 ) { return fabs( ( Mjj[j1][j2] - Mjj[j3][j4] ) / ( Mjj[j1][j2] + Mjj[j3][j4] ) ); };

private:
    // ROOT related
    TFile* tree_file;
    TTree* _event_Tree;
    TTree* _rech1_Tree;

    // Processor Parameters
    std::string _treeFileName;
    std::stirng _event_treeName;
    std::string _rech1_treeName;
    std::stirng _colName;
    std::ostream _output;
    int  _overwrite;
    int _NJet;

    // Vars fill to tree
    //===============================
    // event
    //-------------------------------
    int _eventNum;
    double _Evisible;
    double _Emiss;
    
    // particle
    double _eta;
    double _cosThetaj;
    double _yij;
    double mll;
    double mz;
    double mh;
    double _mrecoil;
    double _delta_mll_mz;
    double _delta_mrecoil_mh;
    double _leading_jet_pT;

    // jet
    DoubleVec _vjTag;

    
    //===============================
    // h1
    double _h1InvMass;
    double _DletaRjj;
    std::string _j1Tag;
    std::string _j2Tag;

    
    // run time vars

    double h1InvMass;
}

#endif
