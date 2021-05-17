#ifndef _RecSinglet_hh_
#define _RecSinglet_hh_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "marlin/Processor.h"
#include "LCIOSTLTypes.h"

#include "TTree.h"
#include "TFile.h"


class TTree;


class RecSinglet : public marlin::Processor
{
public:
    Processor* newProcessor() { return new RecSinglet ; };

    RecSinglet();

    ~RecSinglet();

    void init();

    void processEvent( LCEvent * evtP );

    void end();

protected:
    // ROOT related
    TFile *tree_file;
    TTree *_outputTree;

    // Processor Parameters
    std::string _treeFileName;
    std::string _treeName;
    std::string _colName;
    std::ostream _output;
    int _overwrite;

    // Variables fill to tree
    int _eventNum;
    double _h1InvMass; // invariant mass of h1 (singlet)
    double _h1Psqr; // Singlet momentum square
    double _h1E; // Singlet Kinetic Energy
    double _DeltaR; // DeltaR of j1 j2
    double _j1Tag; // Tag of j1
    double _j2Tag; // Tag of j2
    double _EjCut; // Cut of jet Energy

    // run time variables
    int NJetsNum;
    int alcfiplus; // Algorithm ID of lcfiplus
    int ibtag; // Parameter index of BTag
    int ictag; // Parameter index of CTag
    int iltag; // value of light jet param
    int icategory; // Parameter index of Category 

};

#endif