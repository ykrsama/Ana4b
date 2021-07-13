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

    _treeName2 = "bg2j";

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

    _deltaMCut=1; // GeV
    registerProcessorParameter( "deltaMCut",
        "delta M cut. If set deltaMCut < 0, no cut.",
        _deltaMCut,
        _deltaMCut);

    _DeltaRMax=1;
    registerProcessorParameter( "DeltaRMax",
        "Delta R max. If set deltaMCut < 0, no cut.",
        _DeltaRMax,
        _DeltaRMax);
    
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

    _DeltaRjlMax=1;
    registerProcessorParameter( "DeltaRjlMax",
        "DeltaR_{jet,lepton}.If set < 0 turn off this cut.",
        _DeltaRjlMax,
        _DeltaRjlMax);

}


Exotic4b::~Exotic4b()
{
    delete tree_file;
}


double Exotic4b::getMassjj( double j1E,
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
    //_outputTree->Branch("dMsM", &_dMsM);
    // _outputTree->Branch("h1Psqr", &_h1Psqr);
    // _outputTree->Branch("h1E", &_h1E);
    _outputTree->Branch("DeltaR", &_DeltaR);
    _outputTree->Branch("j1Tag",&_j1Tag);
    _outputTree->Branch("j2Tag",&_j2Tag);

    _outputTree2 = new TTree(_treeName2.c_str(), _treeName2.c_str() );
    _outputTree2->Branch("h2InvMass", &_h2InvMass);
    _outputTree2->Branch("DeltaR", &_DeltaR);
    _outputTree2->Branch("j1Tag",&_j1Tag);
    _outputTree2->Branch("j2Tag",&_j2Tag);

    Leptons[0]=11; Leptons[1]=-11; Leptons[2]=13; Leptons[3]=-13; Leptons[4]=15; Leptons[5]=-15;

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
        std::vector<TLorentzVector> Vlepton_table;

        // read lcio data
        try
        {
            // Get Collection
            LCCollection* colJet = evtP->getCollection(_colName);
            LCCollection* MCPart = evtP->getCollection("MCParticle");
            NMCP = MCPart->getNumberOfElements();
            NJetsNum = colJet->getNumberOfElements();
            _eventNum = evtP->getEventNumber();

            // cut Energy of jets
            for (int i = 0; i < NJetsNum; i++)
            {
                // reconstructed jet particle
                ReconstructedParticle* recP = dynamic_cast<ReconstructedParticle*>(colJet->getElementAt(i));
                if ( recP->getEnergy() < _EjCut) return;
            }

            // Get MCParticle
            for(int i = 0; i < NMCP; i++)
            {
                MCParticle * a_MCP = dynamic_cast<MCParticle*>(MCPart->getElementAt(i));
                if(a_MCP->getParents().size() == 0) continue;
                
                MCPID = a_MCP->getPDG();
                int *p = std::find(Leptons, Leptons + 6, MCPID);
                if(p == Leptons + 6) continue;

                MCPVertex[0] = a_MCP->getVertex()[0];
                MCPVertex[1] = a_MCP->getVertex()[1];
                MCPVertex[2] = a_MCP->getVertex()[2];
                if(MCPVertex[0] > 0.0001 || MCPVertex[1] > 0.0001 || MCPVertex[2] > 0.0001) continue;

                MCPEn = a_MCP->getEnergy();
                MCPP[0] = a_MCP->getMomentum()[0];
                MCPP[1] = a_MCP->getMomentum()[1];
                MCPP[2] = a_MCP->getMomentum()[2];
                TLorentzVector Vlepton;
                Vlepton.SetPxPyPzE( MCPP[0],
                                    MCPP[1],
                                    MCPP[2],
                                    MCPEn);
                Vlepton_table.push_back( Vlepton );
            }
            Vlepton_table_size = Vlepton_table.size();

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
        
        // exclude lepton jet
        if ( _DeltaRjlMax >=0 )
        {
            for (int i = 0; i < NJetsNum; i++)
            {
                Vj1.SetPxPyPzE( vjPx[i],
                                vjPy[i],
                                vjPz[i],
                                vjE[i]);

                for (int j = 0; j < Vlepton_table_size; j++)
                {
                    float DeltaRjl = Vj1.DeltaR(Vlepton_table.at(j));
                    if ( DeltaRjl < _DeltaRjlMax ) return;
                }
            }
        }

        if ( NJetsNum == 4 ) {
            // reconstruct singlet scaler
            for (int i = 0; i < NJetsNum; i++)
            {
                for (int j = i + 1; j < NJetsNum; j++)
                {
                    // mass jj
                    Mjj[i][j] = getMassjj(vjE.at(i),
                                        vjPx.at(i),
                                        vjPy.at(i),
                                        vjPz.at(i),
                                        vjE.at(j),
                                        vjPx.at(j),
                                        vjPy.at(j),
                                        vjPz.at(j));
                    // calculate deltaR
                    Vj1.SetPxPyPzE( vjPx[ i ],
                                    vjPy[ i ],
                                    vjPz[ i ],
                                    vjE[ i ]);
                    Vj2.SetPxPyPzE( vjPx[ j ],
                                    vjPy[ j ],
                                    vjPz[ j ],
                                    vjE[ j ]);
                    DeltaRjj[i][j] = Vj1.DeltaR(Vj2);
                }
            }

            _deltaM = getDeltaM(0, 1, 2, 3);
            //_dMsM = getdMsM(0, 1, 2, 3);
            jIndex[0] = 0;
            jIndex[1] = 1;
            jIndex[2] = 2;
            jIndex[3] = 3;
            if ( getDeltaM(0, 2, 1, 3) < _deltaM )
            //if ( getdMsM(0, 2, 1, 3) < _dMsM )
            {
                _deltaM = getDeltaM(0, 2, 1, 3);
                //_dMsM = getdMsM(0, 2, 1, 3);
                jIndex[0] = 0;
                jIndex[1] = 2;
                jIndex[2] = 1;
                jIndex[3] = 3;
            }
            if ( getDeltaM(0, 3, 1, 2) < _deltaM )
            //if ( getdMsM(0, 3, 1, 2) < _dMsM )
            {
                _deltaM = getDeltaM(0, 3, 1, 2);
                //_dMsM = getdMsM(0, 3, 1, 2);
                jIndex[0] = 0;
                jIndex[1] = 3;
                jIndex[2] = 1;
                jIndex[3] = 2;
            }
            if ( _deltaMCut >=0 && _deltaM >= _deltaMCut ) return;

            // check Rm
            _Rm = getRm(jIndex[0], jIndex[1], jIndex[2], jIndex[3]);
            if ( _Rm >= _RmCut ) return;

            // check DeltaR_jj
            if ( _DeltaRMax >=0 )
            {
                if ( DeltaRjj[ jIndex[0] ][ jIndex[1] ] >= _DeltaRMax ) return;
                if ( DeltaRjj[ jIndex[2] ][ jIndex[3] ] >= _DeltaRMax ) return;
            }

            // fill the first singlet
            for (int i = 0; i < 4; i += 2)
            {
                j1I = jIndex[i];
                j2I = jIndex[i + 1];

                _h1InvMass = Mjj[ j1I ][ j2I ];
                _DeltaR = DeltaRjj[ j1I ][ j2I ];

                _j1Tag = vjTag[ j1I ];
                _j2Tag = vjTag[ j2I ];

                _outputTree->Fill();
            }
        }

        if ( NJetsNum == 2 )
        {
            _h2InvMass = getMassjj( vjE.at(0),
                                    vjPx.at(0),
                                    vjPy.at(0),
                                    vjPz.at(0),
                                    vjE.at(1),
                                    vjPx.at(1),
                                    vjPy.at(1),
                                    vjPz.at(1));
            Vj1.SetPxPyPzE( vjPx.at(0),
                            vjPy.at(0),
                            vjPz.at(0),
                            vjE.at(0) );
            Vj2.SetPxPyPzE( vjPx.at(1),
                            vjPy.at(1),
                            vjPz.at(1),
                            vjE.at(1));
            _DeltaR = Vj1.DeltaR(Vj2);
            _j1Tag = vjTag.at(0);
            _j2Tag = vjTag.at(1);

            _outputTree2->Fill();
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