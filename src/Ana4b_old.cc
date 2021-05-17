#include "Ana4b.hh"
#include "EVENT/LCCollection.h"
#include "EVENT/LCParameters.h"
#include "EVENT/ReconstructedParticle.h"
#include "UTIL/PIDHandler.h"
#include <math.h>


Ana4b_old Ana4b_instance;


Ana4b_old::Ana4b_old()
    : Processor("Ana4b"),
    _output(0),
    totalRec4b(0)
{
    _description = "Channel: h2 > 4b. Save b-tagging efficiency and h1 invariant mass";

    _treeFileName = "Ana4b.root";
    registerProcessorParameter( "TreeOutputFile",
        "The name of the file to which the ROOT tree will be written",
        _treeFileName,
        _treeFileName);

    _treeName = "ana";
    registerProcessorParameter( "TreeName",
        "The name of the ROOT tree",
        _treeName,
        _treeName);

    _treeName2 = "correct";
    registerProcessorParameter( "TreeName2",
        "The name of the ROOT tree about correctly reconstructed event",
        _treeName2,
        _treeName2);
    

    _colName  = "RefinedJets";
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


Ana4b_old::~Ana4b_old() {
    delete tree_file;
}


void Ana4b_old::init() {
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
    // Algorithm lcfiplus parameters
    _outputTree->Branch("BTags", &_bTags);
    _outputTree->Branch("CTags", &_cTags);
    _outputTree->Branch("Categories", &_categories);
    // Algorithm yth parameters
    //_outputTree->Branch("yths", &yths);
    _outputTree->Branch("y01", &_y01s);
    _outputTree->Branch("y12", &_y12s);
    _outputTree->Branch("y23", &_y23s);
    _outputTree->Branch("y34", &_y34s);
    _outputTree->Branch("y45", &_y45s);
    _outputTree->Branch("y56", &_y56s);
    _outputTree->Branch("y67", &_y67s);
    _outputTree->Branch("y78", &_y78s);
    _outputTree->Branch("y89", &_y89s);
    _outputTree->Branch("y910", &_y910s);
    // Calculatepd
    _outputTree->Branch("RecBNum", &_recBNum, "RecBNum/i");
    _outputTree->Branch("IsRec4b", &_isRec4b, "IsRec4b/b");

    _outputTree2 = new TTree(_treeName2.c_str(), _treeName2.c_str() );
    _outputTree2->Branch("Mass1", &_mass1, "Mass1/D");
    _outputTree2->Branch("Mass2", &_mass2, "Mass2/D");
    _outputTree2->Branch("Mass3", &_mass3, "Mass3/D");
    _outputTree2->Branch("Mass4", &_mass4, "Mass4/D");
    _outputTree2->Branch("h2E", &_h2E, "h2E/D");
    _outputTree2->Branch("h2Psqr", &_h2Psqr, "h2Psqr/D");
    _outputTree2->Branch("h2InvMass", &_h2InvMass, "h2InvMass/D");
    _outputTree2->Branch("E1", &_E1, "E1/D");
    _outputTree2->Branch("E2", &_E2, "E2/D");
    _outputTree2->Branch("E3", &_E3, "E3/D");
    _outputTree2->Branch("E4", &_E4, "E4/D");
}


void Ana4b_old::processEvent( LCEvent *evtP ) {
    if (evtP) {
        try {
            // Clear containers
            _bTags.clear();
            _cTags.clear();
            _categories.clear();
            yith.clear();
            _y01s.clear();
            _y12s.clear();
            _y23s.clear();
            _y34s.clear();
            _y45s.clear();
            _y56s.clear();
            _y67s.clear();
            _y78s.clear();
            _y89s.clear();
            _y910s.clear();
            iyth.clear();
            // Init variables
            _recBNum = 0;
            _h2E = 0;
            h2Px = 0;
            h2Py = 0;
            h2Pz = 0;

            // Get Collection
            LCCollection* colJet = evtP->getCollection(_colName);
            NJetsNum = colJet->getNumberOfElements();

            // Handle PID information
            PIDHandler pidH (colJet);
            // get algorithm IDs
            alcfiplus = pidH.getAlgorithmID("lcfiplus"); // [id: 1]   lcfiplus - params:  BTag CTag Category
            ayth = pidH.getAlgorithmID("yth"); // [id: 2]   yth - params:  y01 y12 y23 y34 y45 y56 y67 y78 y89 y910
            // get lcfiplus parameter indices
            ibtag = pidH.getParameterIndex(alcfiplus, "BTag");
            ictag = pidH.getParameterIndex(alcfiplus, "CTag");
            icategory = pidH.getParameterIndex(alcfiplus, "Category");
            // get yth parameter indices
            for (int i = 0; i < 10; i++) {
                std::stringstream ss;
                ss << "y" << i << (i + 1);
                iyth.push_back( pidH.getParameterIndex( ayth, ss.str() ) );
            }
            
            // Search PID, read PID information
            try {
                for (int i = 0; i < NJetsNum; i++) {
                    ReconstructedParticle* recP = dynamic_cast<ReconstructedParticle*>(colJet->getElementAt(i));
                    // mass
                    switch (i) {
                        case 0 :
                            _mass1 = recP->getMass();
                            _E1 = recP->getEnergy();
                            break;
                        
                        case 1 :
                            _mass2 = recP->getMass();
                            _E2 = recP->getEnergy();
                            break;

                        case 2 :
                            _mass3 = recP->getMass();
                            _E3 = recP->getEnergy();
                            break;

                        case 3 :
                            _mass4 = recP->getMass();
                            _E4 = recP->getEnergy();
                            break;
                        
                        default :
                            break;
                    }
                    _h2E += recP->getEnergy();
                    h2Px += recP->getMomentum()[0];
                    h2Py += recP->getMomentum()[1];
                    h2Pz += recP->getMomentum()[2];
                    // lcfiplus parametes
                    const ParticleID &pid_lcfiplus = pidH.getParticleID(recP, alcfiplus);
                    _bTags.push_back(pid_lcfiplus.getParameters()[ibtag]);
                    _cTags.push_back(pid_lcfiplus.getParameters()[ictag]);
                    _categories.push_back(pid_lcfiplus.getParameters()[icategory]);
                    // yth parameters
                    const ParticleID &pid_yth = pidH.getParticleID(recP, ayth);
                    _y01s.push_back(pid_yth.getParameters()[ iyth.at(0) ]);
                    _y12s.push_back(pid_yth.getParameters()[ iyth.at(1) ]);
                    _y23s.push_back(pid_yth.getParameters()[ iyth.at(2) ]);
                    _y34s.push_back(pid_yth.getParameters()[ iyth.at(3) ]);
                    _y45s.push_back(pid_yth.getParameters()[ iyth.at(4) ]);
                    _y56s.push_back(pid_yth.getParameters()[ iyth.at(5) ]);
                    _y67s.push_back(pid_yth.getParameters()[ iyth.at(6) ]);
                    _y78s.push_back(pid_yth.getParameters()[ iyth.at(7) ]);
                    _y89s.push_back(pid_yth.getParameters()[ iyth.at(8) ]);
                    _y910s.push_back(pid_yth.getParameters()[ iyth.at(9) ]);
                    // backup
                    yth.clear();
                    for (int j = 0; j < 10; j++) {
                       yth.push_back(pid_yth.getParameters()[ iyth.at(j) ]);
                    }
                    yith.push_back(yth);
                }
            }catch (lcio::DataNotAvailableException err) {}

            // Calcuate
            for (int i = 0; i < NJetsNum; i++) {
                if ( _bTags.at(i) > 0.5 && _bTags.at(i) > _cTags.at(i) ) {
                    _recBNum++;
                }
            }
            
            if ( _recBNum == 4 ) {
                _isRec4b = 1;
                totalRec4b++;
            }
            else _isRec4b = 0;

            _h2Psqr = h2Px*h2Px + h2Py*h2Py + h2Pz*h2Pz;
            _h2InvMass = sqrt( _h2E*_h2E - _h2Psqr );


            // process branch EventNum
            _eventNum = evtP->getEventNumber();

            _outputTree->Fill();

            if ( _isRec4b ) _outputTree2->Fill();
        }catch (lcio::DataNotAvailableException err) {}
        
    }
}


void Ana4b_old::end() {
    if (_outputTree) {
        tree_file = _outputTree->GetCurrentFile(); //just in case we switched to a new file
        tree_file->Write();
    }

    recoEff = totalRec4b / (_eventNum + 1.);
    ofstream OutFile("recoEff.log");
    OutFile << recoEff << "\n";
    OutFile.close();

    delete tree_file;
}