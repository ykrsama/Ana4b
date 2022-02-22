#ifndef _Truth4b_hh_
#define _Truth4b_hh_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "marlin/Processor.h"
#include "LCIOSTLTypes.h"

#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"


class TTree;


class Truth4b: public marlin::Processor
{
public:
    Processor* newProcessor() {return new Truth4b ; }

    Truth4b();

    ~Truth4b();

    void init();

    void processEvent( LCEvent * evtP );

    void end();


protected:
    // ROOT related
    TFile *tree_file;
    TTree* _outputTree;

    // Processor Parameters
    std::string _treeFileName;
    std::string _treeName;
    std::string _colName;
    std::ostream _output;
    int _overwrite;

    // Variables fill to tree
    int _eventNum;
    // mass
    double _mass1;
    double _mass2;
    double _mass3;
    double _mass4;
    // energy
    double _E1;
    double _E2;
    double _E3;
    double _E4;
    // momentum
    FloatVec p1Vec;
    FloatVec p2Vec;
    FloatVec p3Vec;
    FloatVec p4Vec;
    // calculated
    double _h2E;
    double _h2Psqr;
    double _h2InvMass;
    double _h1E1;
    double _h1Psqr1;
    double _h1InvMass1;
    double _h1E2;
    double _h1Psqr2;
    double _h1InvMass2;
    double _deltaR1;
    double _deltaR2;
    FloatVec _h2DecayLength;
    //FloatVec _bbMass;

    // run time variables
    int pdgid;
    int partNum;
    int bcount;
    int h1count;
    float recoEff;
    int totalRec4b;
    float h2Px;
    float h2Py;
    float h2Pz;
    float h1Px;
    float h1Py;
    float h1Pz;
    TLorentzVector v1;
    TLorentzVector v2;
    TLorentzVector v3;
    TLorentzVector v4;

};


#endif // _Truth4b_hh_
