#include "Truth4b.hh"
#include "EVENT/LCCollection.h"
#include "EVENT/LCParameters.h"
#include "EVENT/MCParticle.h"
//#include "UTIL/PIDHandler.h"
#include <math.h>


Truth4b Truth4b_instance;


Truth4b::Truth4b()
    : Processor("Truth4b"),
    _output(0),
    totalRec4b(0)
{
    _description = "Channel: h2 > 4b. Save b-tagging efficiency and h1 invariant mass";

    _treeFileName = "Truth4b.root";
    registerProcessorParameter( "TreeOutputFile",
        "The name of the file to which the ROOT tree will be written",
        _treeFileName,
        _treeFileName);

    _treeName = "truth";
    registerProcessorParameter( "TreeName",
        "The name of the ROOT tree",
        _treeName,
        _treeName);


    _colName  = "MCParticle";
    registerProcessorParameter( "CollectionName",
        "The name of Jet Collection.",
        _colName,
        _colName);

    _overwrite=0;
    registerProcessorParameter( "OverwriteFile",
        "If zero an already existing file will not be overwritten.",
        _overwrite,
        _overwrite);
}


Truth4b::~Truth4b() {
    delete tree_file;
}


void Truth4b::init() {
    printParameters();

    tree_file = new TFile(_treeFileName.c_str(), (_overwrite ? "RECREATE" : "UPDATE") );

    if ( !tree_file->IsOpen() ) {
        delete tree_file;
        tree_file = new TFile(_treeFileName.c_str(), "NEW");
    }

    _outputTree = new TTree( _treeName.c_str(), _treeName.c_str() );
    _outputTree->Branch("EventNum", &_eventNum, "EventNum/I");

    // reconstructedParticle class
    _outputTree->Branch("Mass1", &_mass1, "Mass1/D");
    _outputTree->Branch("Mass2", &_mass2, "Mass2/D");
    _outputTree->Branch("Mass3", &_mass3, "Mass3/D");
    _outputTree->Branch("Mass4", &_mass4, "Mass4/D");
    _outputTree->Branch("E1", &_E1, "E1/D");
    _outputTree->Branch("E2", &_E2, "E2/D");
    _outputTree->Branch("E3", &_E3, "E3/D");
    _outputTree->Branch("E4", &_E4, "E4/D");
    _outputTree->Branch("h2E", &_h2E, "h2E/D");
    _outputTree->Branch("h2Psqr", &_h2Psqr, "h2Psqr/D");
    _outputTree->Branch("h2InvMass", &_h2InvMass, "h2InvMass/D");
    _outputTree->Branch("h1E1", &_h1E1, "h1E1/D");
    _outputTree->Branch("h1Psqr1", &_h1Psqr1, "h1Psqr1/D");
    _outputTree->Branch("h1InvMass1", &_h1InvMass1, "h1InvMass1/D");
    _outputTree->Branch("h1E2", &_h1E2, "h1E2/D");
    _outputTree->Branch("h1Psqr2", &_h1Psqr2, "h1Psqr2/D");
    _outputTree->Branch("h1InvMass2", &_h1InvMass2, "h1InvMass2/D");
    _outputTree->Branch("deltaR1", &_deltaR1, "deltaR1/D");
    _outputTree->Branch("deltaR2", &_deltaR2, "deltaR2/D");

}


void Truth4b::processEvent( LCEvent *evtP ) {
    if (evtP) {
        try {
            // Init variables
            _h2E = 0;
            h2Px = 0;
            h2Py = 0;
            h2Pz = 0;
            bcount = 0;
            h1count = 0;

            // Get Collection
            LCCollection* col = evtP->getCollection(_colName);
            partNum = col->getNumberOfElements();
            
            // Search truth information
                for (int i = 0; i < partNum; i++ ) {
                    MCParticle* simP = dynamic_cast<MCParticle*>(col->getElementAt(i));
                    pdgid = simP->getPDG();
                    if ( h1count == 2 ) break;
                    if ( pdgid == 9999992 ) {
                        h1count++;
                        MCParticleVec h1Daughters = simP->getDaughters();
                        // init vars
                        _h1E1 = 0;
                        _h1E2 = 0;
                        h1Px = 0;
                        h1Py = 0;
                        h1Pz = 0;
                        for (MCParticleVec::iterator itr = h1Daughters.begin(); itr != h1Daughters.end(); itr++) {
                            switch ( bcount ) {
                                case 0 :
                                    _mass1 = (*itr)->getMass();
                                    _E1 = (*itr)->getEnergy();
                                    break;

                                case 1 :
                                    _mass2 = (*itr)->getMass();
                                    _E2 = (*itr)->getEnergy();
                                    break;

                                case 2 :
                                    _mass3 = (*itr)->getMass();
                                    _E3 = (*itr)->getEnergy();
                                    break;

                                case 3 :
                                    _mass4 = (*itr)->getMass();
                                    _E4 = (*itr)->getEnergy();
                                    break;

                                default :
                                    break;
                            }
                            bcount++;
                            // calculate h1
                            switch ( h1count ) {
                                case 1 : 
                                    _h1E1 += (*itr)->getEnergy();
                                    break;
                                case 2 :
                                    _h1E2 += (*itr)->getEnergy();
                                    break;
                                default :
                                    break;
                            }
                            h1Px += (*itr)->getMomentum()[0];
                            h1Py += (*itr)->getMomentum()[1];
                            h1Pz += (*itr)->getMomentum()[2];
                            // calculate h2
                            _h2E += (*itr)->getEnergy();
                            h2Px += (*itr)->getMomentum()[0];
                            h2Py += (*itr)->getMomentum()[1];
                            h2Pz += (*itr)->getMomentum()[2];
                        }
                        // calculate h1
                        switch ( h1count ) {
                            case 1 :
                                _h1Psqr1 = h1Px * h1Px + h1Py * h1Py + h1Pz * h1Pz;
                                _h1InvMass1 = sqrt(_h1E1 * _h1E1 - _h1Psqr1);
                                break;
                            case 2 :
                                _h1Psqr2 = h1Px * h1Px + h1Py * h1Py + h1Pz * h1Pz;
                                _h1InvMass2 = sqrt(_h1E2 * _h1E2 - _h1Psqr2);
                                break;
                            default :
                                break;
                        }
                    }
                }
            
                // Calcuate
                _h2Psqr = h2Px * h2Px + h2Py * h2Py + h2Pz * h2Pz;
                _h2InvMass = sqrt( _h2E * _h2E - _h2Psqr );

                // process branch EventNum
                _eventNum = evtP->getEventNumber();

                _outputTree->Fill();

        }catch (lcio::DataNotAvailableException err) {}
        
    }
}


void Truth4b::end() {
    if (_outputTree) {
        tree_file = _outputTree->GetCurrentFile(); //just in case we switched to a new file
        tree_file->Write();
    }

    delete tree_file;
}