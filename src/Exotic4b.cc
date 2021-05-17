#include "Exotic4b.hh"


Exotic4b Exotic4b_instance;


Exotic4b::Exotic4b()
    : Processor("Exotic4b"),
    _output(0)
{
    _description = "Reconstruct the light singlet h1. Channel: h2 > h1h1 > 4b.";

    _treeFileName = "Exotic4b.root";
    registerProcessorParameter( "TreeOutputFile",
        "The name of the file to which the ROOT tree will be written",
        _treeFileName,
        _treeFileName);

    _treeName = "exotic";
    registerProcessorParameter( "TreeName",
        "The name of the ROOT tree",
        _treeName,
        _treeName);

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
    
    _EjCut=5; // GeV
    registerProcessorParameter( "EjCut",
        "4 jets satisfy Energy >= EjCut",
        _EjCut,
        _EjCut);

    _RmCut=0.1;
    registerProcessorParameter( "RmCut",
        "Rm < RmCut",
        _RmCut,
        _RmCut);
}


Exotic4b::~Exotic4b()
{
    delete tree_file;
}


double Exotic4b::getMassjj(  double j1E,
                double j1Px,
                double j1Py,
                double j1Pz,
                double j2E,
                double j2Px,
                double j2Py,
                double j2Pz)
{
    double Ejj = j1E + j2E;
    double Pjjx = j1Px + j2Px;
    double Pjjy = j1Py + j2Py;
    double Pjjz = j1Pz + j2Pz;
    double Pjjsqr = Pjjx * Pjjx + Pjjy * Pjjy + Pjjz * Pjjz;
    double Mjj = sqrt( Ejj * Ejj - Pjjsqr );
    return Mjj;
}


void Exotic4b::init()
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
    _outputTree->Branch("deltaM", &_deltaM);
    _outputTree->Branch("Rm", &_Rm);
    // _outputTree->Branch("h1Psqr", &_h1Psqr);
    // _outputTree->Branch("h1E", &_h1E);
    _outputTree->Branch("DeltaR", &_DeltaR);
    _outputTree->Branch("j1Tag",&_j1Tag);
    _outputTree->Branch("j2Tag",&_j2Tag);

}


void Exotic4b::processEvent( LCEvent *evtP )
{
    if (evtP)
    {
        // declear containers
        DoubleVec vjE; // x0 of jet i of 4
        DoubleVec vjPx; // momentum x of jet i of 4
        DoubleVec vjPy;
        DoubleVec vjPz;
        StringVec vjTag;

        // read lcio data
        try
        {
            // Get Collection
            LCCollection* colJet = evtP->getCollection(_colName);
            NJetsNum = colJet->getNumberOfElements();
            _eventNum = evtP->getEventNumber();
            // check NJetsNum
            if ( NJetsNum != 4 ) return;
            // cut Energy of jets
            for (int i = 0; i < NJetsNum; i++)
            {
                // reconstructed jet particle
                ReconstructedParticle* recP = dynamic_cast<ReconstructedParticle*>(colJet->getElementAt(i));
                if ( recP->getEnergy() < _EjCut) return;
            }
            // Handle PID information
            PIDHandler pidH (colJet);
            // Get algorithm IDs
            alcfiplus = pidH.getAlgorithmID("lcfiplus");
            // Get lcfiplus parameter indicies
            ibtag = pidH.getParameterIndex(alcfiplus, "BTag");
            ictag = pidH.getParameterIndex(alcfiplus, "CTag");
            // Get jets PID
            for (int i = 0; i < NJetsNum; i++)
            {
                // reconstructed jet particle
                ReconstructedParticle* recP = dynamic_cast<ReconstructedParticle*>(colJet->getElementAt(i));
                // jets 4-momentum
                vjE.push_back( recP->getEnergy() );
                vjPx.push_back( recP->getMomentum()[0] );
                vjPy.push_back( recP->getMomentum()[1] );
                vjPz.push_back( recP->getMomentum()[2] );
                // lcfipuls parameters of recP
                const ParticleID &pid_lcfiplus = pidH.getParticleID(recP, alcfiplus);
                bTag = pid_lcfiplus.getParameters()[ibtag];
                cTag = pid_lcfiplus.getParameters()[ictag];
                lTag = 1 - bTag - cTag;
                // find most likely jet tag
                tempTagParam = bTag;
                tempTag = "b";
                if ( cTag > tempTagParam )
                {
                    tempTagParam = cTag;
                    tempTag = "c";
                }
                if ( lTag > tempTagParam )
                {
                    tempTagParam = lTag;
                    tempTag = "l";
                }
                vjTag.push_back( tempTag );
            }
        }catch (lcio::DataNotAvailableException err) {}
        
        // reconstruct singlet scaler
        for (int i = 0; i < NJetsNum; i++)
        {
            for (int j = i + 1; j < NJetsNum; j++)
            {
                Mjj[i][j] = getMassjj(vjE.at(i),
                                    vjPx.at(i),
                                    vjPy.at(i),
                                    vjPz.at(i),
                                    vjE.at(j),
                                    vjPx.at(j),
                                    vjPy.at(j),
                                    vjPz.at(j));
            }
        }
        _deltaM = getDeltaM(0, 1, 2, 3);
        jIndex[0] = 0;
        jIndex[1] = 1;
        jIndex[2] = 2;
        jIndex[3] = 3;
        if ( getDeltaM(0, 2, 1, 3) < _deltaM )
        {
            _deltaM = getDeltaM(0, 2, 1, 3);
            jIndex[0] = 0;
            jIndex[1] = 2;
            jIndex[2] = 1;
            jIndex[3] = 3;
        }
        if ( getDeltaM(0, 3, 1, 2) < _deltaM )
        {
            _deltaM = getDeltaM(0, 3, 1, 2);
            jIndex[0] = 0;
            jIndex[1] = 3;
            jIndex[2] = 1;
            jIndex[3] = 2;
        }
        _Rm = getRm(jIndex[0], jIndex[1], jIndex[2], jIndex[3]);
        if ( _Rm >= _RmCut ) return;
        
        // fill the first singlet
        for (int i = 0; i < 4; i += 2)
        {
            j1I = jIndex[i];
            j2I = jIndex[i + 1];
            
            _h1InvMass = Mjj[ j1I ][ j2I ];
            _j1Tag = vjTag[ j1I ];
            _j2Tag = vjTag[ j2I ];
            // calculate deltaR
            Vj1.SetPxPyPzE( vjPx[ j1I ],
                            vjPy[ j1I ],
                            vjPz[ j1I ],
                            vjE[ j1I ]);
            Vj2.SetPxPyPzE( vjPx[ j2I ],
                            vjPy[ j2I ],
                            vjPz[ j2I ],
                            vjE[ j2I ]);
            _DeltaR = Vj1.DeltaR(Vj2);

            _outputTree->Fill();
        }
        

    }

}


void Exotic4b::end()
{
    if (_outputTree)
    {
        tree_file = _outputTree->GetCurrentFile(); // just in case we switched to a new file.
        tree_file->Write();
    }

    delete tree_file;
}