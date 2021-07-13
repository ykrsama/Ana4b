#ifndef _BG2j_hh_
#define _BG2j_hh_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "marlin/Processor.h"
#include "LCIOSTLTypes.h"

#include "TTree.h"
#include "TFile.h"


class TTree;


class BG2j : public marlin::Processor
{
public:
    Processor* newProcessor() {return new BG2j; };

    BG2j();

    ~BG2j();

    void init();

    void processEvent( LCEvent * evtP );
};

#endif
