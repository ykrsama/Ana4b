#ifndef _Ana4b_hh_
#define _Ana4b_hh_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "marlin/Processor.h"
#include "LCIOSTLTypes.h"

#include "TTree.h"
#include "TFile.h"


class TTree;


class Ana4b_old : public marlin::Processor
{
public:
    Processor* newProcessor() {return new Ana4b_old ; };

    Ana4b_old();

    ~Ana4b_old();

    void init();

    void processEvent( LCEvent * evtP );

    void end();


protected:
    // ROOT related
    TFile *tree_file;
    TTree* _outputTree;
    TTree* _outputTree2;

    // Processor Parameters
    std::string _treeFileName;
    std::string _treeName;
    std::string _treeName2;
    std::string _colName;
    std::ostream _output;
    int _overwrite;

    // Variables fill to tree
    int _eventNum;
    // mass of b jets
    double _mass1;
    double _mass2;
    double _mass3;
    double _mass4;
    // energy of b jets
    double _E1;
    double _E2;
    double _E3;
    double _E4;
    // momentum of b jets
    FloatVec p1Vec;
    FloatVec p2Vec;
    FloatVec p3Vec;
    FloatVec p4Vec;
    // algorithm lcfiplus parameters
    FloatVec _bTags; // size = 4
    FloatVec _cTags;
    FloatVec _categories;
    // algorithm yth parameters
    FloatVec _y01s; // size=4
    FloatVec _y12s;
    FloatVec _y23s;
    FloatVec _y34s;
    FloatVec _y45s;
    FloatVec _y56s;
    FloatVec _y67s;
    FloatVec _y78s;
    FloatVec _y89s;
    FloatVec _y910s;
    // calculated
    int _recBNum;
    int _isRec4b;
    double _h2Psqr;
    double _h2E;
    double _h2InvMass;
    double _h11Psqr;
    double _h12Psqr;
    double _h11E;
    double _h12E;
    double _h11InvMass;
    double _h12InvMass;

    // run time variables
    int NJetsNum;
    FloatVec yth; // size = 10
    std::vector<FloatVec> yith; // size = 4 * 10
    float recoEff;
    int totalRec4b;
    float h2Px;
    float h2Py;
    float h2Pz;

    // PIDHandler related
    int alcfiplus; // Algorithm ID of lcfiplus
    int ayth; // Algorithm ID of yth
    int ibtag; // Parameter index of BTag
    int ictag; // Parameter index of CTagf
    int icategory; // Parameter index of Categroy
    IntVec iyth; // Parameter indedx of yth (iyth[0] = iy01, iyth[1] = iy12)


};


#endif // _Ana4b_hh_
