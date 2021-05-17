#ifndef _Exotic4b_hh_
#define _Exotic4b_hh_

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
#include "UTIL/PIDHandler.h"

#include "TTree.h"
#include "TFile.h"


class TTree;


class Exotic4b : public marlin::Processor
{
public:
    Processor* newProcessor() { return new Exotic4b ; };

    Exotic4b();

    ~Exotic4b();

    void init();

    void processEvent( LCEvent * evtP );

    void end();

private:
    double Get_massjj(  double j1E,
                    double j1Px,
                    double j1Py,
                    double j1Pz,
                    double j2E,
                    double j2Px,
                    double j2Py,
                    double j2Pz );


    double Get_deltaM( int j1, int j2, int j3, int j4) { return fabs( Mjj[j1][j2] - Mjj[j3][j4] ); };

    double Get_Rm( int j1, int j2, int j3, int j4 ) { return fabs( ( Mjj[j1][j2] - Mjj[j3][j4] ) / ( Mjj[j1][j2] + Mjj[j3][j4] ) ); };

private:
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
    double _deltaM;
    double _Rm;
    // double _h1Psqr; // Singlet momentum square
    // double _h1E; // Singlet Kinetic Energy
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

    double bTag;
    double cTag;
    double lTag;
    double tempTagParam; // temperory tag probability
    std::string tempTag;

    double Mjj[4][4]; // matrix of massjj. Is a upper triangular matrix.
    int jetIndex[4];
};

#endif