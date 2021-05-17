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
    _description = "Reconstruct the light singlet h1. Channel: h2 > h1h1 > 4b.";

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

    tree_file = new TFile(_treeFileName.c_str(), (_overwrite ? "RECREATE" : "UPDATE") );

    if ( !tree_file->IsOpen() )
    {
        delete tree_file;
        tree_file = new TFile(_treeFileName.c_str(), "NEW");
    }

    _outputTree = new TTree( _treeName.c_str(), _treeName.c_str() );
    
    // Register Branch
    _outputTree->Branch("EventNum", &_eventNum, "EventNum/I");
    // Singlet information
    _outputTree->Branch("h1InvMass", &_h1InvMass);
    _outputTree->Branch("h1Psqr", &_h1Psqr);
    _outputTree->Branch("h1E", &_h1E);
    _outputTree->Branch("DeltaR", &_DeltaR);
    _outputTree->Branch("j1Tag",&_j1Tag);
    _outputTree->Branch("j2Tag",&_j2Tag);
    _outputTree->Branch("EjCut",&_EjCut);
    
}


void RecSinglet::processEvent( LCEvent *evtP )
{
    if (evtP)
    {
        try
        {
            // Get Collection
            LCCollection* colJet = evtP->getCollection(_colName);
            NJetsNum = colJet->getNumberOfElements();
            if ( NJetsNum != 4 ) return;
            // Handle PID information
            PIDHandler pidH (colJet);
            // Get algorithm IDs
            alcfiplus = pidH.getAlgorithmID("lcfiplus");
            // Get lcfiplus parameter indicies
            ibtag = pidH.getParameterIndex(alcfiplus, "BTag");
            ictag = pidH.getParameterIndex(alcfiplus, "CTag");
            // Get PID
            for (int i = 0; i < NJetsNum; i++)
            {
                ReconstructedParticle* recP = dynamic_cast<ReconstructedParticle*>(colJet->getElementAt(i));
                // mass
            }
        }catch (lcio::DataNotAvailableException err) {}
    }

}


void RecSinglet::end()
{
    if (_outputTree)
    {
        tree_file = _outputTree->GetCurrentFile(); // just in case we switched to a new file.
        tree_file->Write();
    }

    delete tree_file;
}