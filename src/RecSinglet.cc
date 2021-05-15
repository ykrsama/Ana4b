#include "RecSinglet.hh"
#include "EVENT/LCCollection.h"
#include "EVENT/LCParameters.h"
#include "EVENT/ReconstructedParticle.h"
#include "UTIL/PIDHandler.h"
#include <math.h>


RecSinglet RecSinglet_instance;


RecSinglet::RecSinglet()
    : Processor("RecSinglet"),
    _output(0)
{
    _description = "Channel: h2 > h1h1 > 4b.";

    _treeFileName = "RecSinglet.root";
    registerProcessorParameter( "TreeOutputFile",
        "The name of the file to which the ROOT tree will be written",
        _treeFileName,
        _treeFileName);

    _colName = "RefinedJets";
    registerProcessorParameter( "CollctionName",
        "The name of Jet Collection.",
        _colName,
        _colName);

    _overwrite=0;
    registerProcessorParameter( "OverwriteFile",
        "If zero an already existing file will not be overwritten.",
        _overwrite,
        _overwrite);
}


RecSinglet::~RecSinglet()
{
    delete tree_file;
}


void RecSinglet::init()
{
    printParameters();
}