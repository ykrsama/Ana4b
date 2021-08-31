#include "LSSAna.hh"

LSSAna LSSAna_instance;

LSSAna::LSSAna() : Processor("LSSAna"),
                   _output(0) {
    _description = "LSS Ana";

    _treeFileName = "LSSAna.root";
    registerProcessorParameter("TreeOutputFile",
        "The name of the file to which the ROOT tree will be written",
        _treeFileName,
        _treeFileName);

    _treeName_event = "event";
    registerProcessorParameter("EventThreeName",
        "event level",
        _treeName_event,
        _treeName_event);

    _treeName_rech1 = "rech1";
    registerProcessorParameter("rech1TreeName",
        "reconstructed h1",
        _treeName_rech1,
        _treeName_rech1);

    _colName = "RefinedJets";
    registerProcessorParameter( "CollectionName",
        "The name of Jet Collection.",
        _colName,
        _colName);

    _overwrite = 1;
    registerProcessorParameter( "OverwriteFile",
        "If zero an already existing file will not be overwriteen.",
        _overwrite,
        _overwrite);
}

LSSAna::~LSSAna() {
    delete tree_file;
}

void LSSAna::init() {
    printParameters();

    tree_file = new TFile(_treeFileName.c_str(), (_overwrite ? "RECREATE" : "UPDATE") );

    if ( !tree_file->IsOpen() ) {
        delete tree_file;
        tree_file = new TFile(_treeFileName.c_str(), "NEW");
    }

    _outputTree_event = new TTree( _treeName_event.c_str(), _treeName_event.c_str());
    _outputTree_rech1 = new TTree( _treeName_rech1.c_str(), _treeName_rech1.c_str());

    // Register Branch
    _outputTree_event->Branch("EventNum", &feventNum, "EventNum/I");
    _outputTree_event->Branch("Evisible", &fEvisible);
    _outputTree_event->Branch("Emiss", &fEmiss);
    _outputTree_event->Branch("Mll", &fmll);
    _outputTree_event->Branch("Ell", &fEll);
    _outputTree_event->Branch("deltaMllMZ", &fdelta_mll_mz);
    _outputTree_event->Branch("Mrecoil", &fmrecoil);
    _outputTree_event->Branch("deltaMrecoilMH", &fdelta_mrecoil_mh);
    // jets
    _outputTree_event->Branch("Jet0pT", &fjet0_pT);
    _outputTree_event->Branch("Jet1pT", &fjet1_pT);
    _outputTree_event->Branch("Jet2pT", &fjet2_pT);
    _outputTree_event->Branch("Jet3pT", &fjet3_pT);
    _outputTree_event->Branch("Jet0E", &fjet0_E);
    _outputTree_event->Branch("Jet1E", &fjet1_E);
    _outputTree_event->Branch("Jet2E", &fjet2_E);
    _outputTree_event->Branch("Jet3E", &fjet3_E);

    _outputTree_rech1 = new TTree(_treeName_rech1.c_str(), _treeName_rech1.c_str());
    _outputTree_rech1->Branch("h1InvMass", &fh1InvMass);

    // init constants
    Leptons[0]=11; Leptons[1]=-11; Leptons[2]=13; Leptons[3]=-13; Leptons[4]=15; Leptons[5]=-15;
    fmz = 91.1876; // PDG2020
    fmh = 125.25;
    S = 240 * 240;

}

void LSSAna::processEvent( LCEvent *evtP ) {
    if (evtP) {
        // declear containers
        IntVec iyth;
        DoubleVec vJpT;
        DoubleVec vJE;
        std::vector<ReconstructedParticle*> vAllJets; //sorted by DeltaRjl
        TVector3 P_sum, VTX, EndP;
        // init vars
        fmll = 0;
        fdelta_mll_mz = fmz;
        ffound_Zlepton = false;
        fEvisible = 0;
        P_sum.SetXYZ(0,0,0);
        // read lcio data
        try {
            LCCollection* colJet = evtP->getCollection(_colName);
            LCCollection* MCPart = evtP->getCollection("MCParticle");
            //LCCollection* ClusterCol = evtP->getCollection("ArborCharged");
            NMCP = MCPart->getNumberOfElements();
            NJetsNum = colJet->getNumberOfElements();
            //fCluster_num = ClusterCol->getNumberOfElements();
            feventNum = evtP->getEventNumber();
            
            //=================================================================
            // Get MCParticle ( i.e. Z leptons)
            //-----------------------------------------------------------------
            for (int i = 0; i < NMCP; i++) {
                MCParticle *MCP_i = dynamic_cast<MCParticle*>(MCPart->getElementAt(i));
                int MCPDG_i = MCP_i->getPDG();
                VTX = MCP_i->getVertex();
                EndP = MCP_i->getEndpoint();
                // sum visible E
                if ( VTX.Mag() < 1 && EndP.Mag() > 1 && fabs(MCPDG_i) != 12 && fabs(MCPDG_i) != 14 && fabs(MCPDG_i) != 16 ) {
                    fEvisible += MCP_i->getEnergy();
                    P_sum += MCP_i->getMomentum();
                }
                // find Z leptons
                int *p = std::find(Leptons, Leptons + 6, MCPDG_i);
                if(p == Leptons + 6) continue;

                for (int j = i + 1; j < NMCP; j++) {
                    MCParticle *MCP_j = dynamic_cast<MCParticle*>(MCPart->getElementAt(j));
                    int MCPID_j = MCP_j->getPDG();
                    if (MCPID_j != -MCPDG_i ) continue; // same flavor differnt charge

                    double mll_temp = getInvMass(MCP_i, MCP_j);
                    if ( fabs( mll_temp - fmz ) < fdelta_mll_mz ) {
                        fmll = mll_temp;
                        fdelta_mll_mz = fabs( fmll - fmz );
                        Zleptons[0] = MCP_i;
                        Zleptons[1] = MCP_j;
                        ffound_Zlepton = true;
                    }
                }
            }

            if ( ! ffound_Zlepton ) {
                std::cout << "[INFO] Zleptons not found in event " << feventNum << ", skip" << std::endl;
                return;
            }
            
            for (int i = 0; i < 2; i++){
                Vl[i].SetPxPyPzE(   Zleptons[i]->getMomentum()[0],
                                    Zleptons[i]->getMomentum()[1],
                                    Zleptons[i]->getMomentum()[2],
                                    Zleptons[i]->getEnergy());
            }
            //============================================================================
            // Get Jets
            //----------------------------------------------------------------------------
            PIDHandler pidH (colJet);
            // get algorithm IDs
            alcfiplus = pidH.getAlgorithmID("lcfiplus");
            ayth = pidH.getAlgorithmID("yth");
            // get lcfiplus parameter indiices
            ibtag = pidH.getParameterIndex(alcfiplus, "BTag");
            ictag = pidH.getParameterIndex(alcfiplus, "CTag");
            // get yth parameter indicies
            for (int i = 0; i < 10; i++ ) {
                std::stringstream ss;
                ss << "y" << i << (i + 1);
                iyth.push_back( pidH.getParameterIndex(ayth, ss.str() ) );
            }
            // get jets PID
            for (int i = 0; i < NJetsNum; i++) {
                ReconstructedParticle* recP = dynamic_cast<ReconstructedParticle*>(colJet->getElementAt(i));
                vAllJets.push_back(recP);
            }
        } catch (lcio::DataNotAvailableException err) {}
        if ( vAllJets.size() < 6 ) { // sanity check
            std::cout << "[INFO] All Jets in event " << feventNum << " less than 6, skip" << std::endl;
            return;
        } 
        //================================================================
        // MCParticle
        //----------------------------------------------------------------
        fMass_visible = sqrt(fEvisible * fEvisible - P_sum.Mag2() );
        fEll = Vl[0].Energy() + Vl[1].Energy();
        fmrecoil = sqrt( S - 2 * sqrt( S ) * fEll + fmll * fmll );
        fdelta_mrecoil_mh = fabs( fmrecoil - fmh );
        //================================================================
        // jets
        //----------------------------------------------------------------
        vAllJets = sortLessDeltaRjl(vAllJets);
        leptonJets[0] = vAllJets.at(0);
        leptonJets[1] = vAllJets.at(1);
        realJets[0] = vAllJets.at(2);
        realJets[1] = vAllJets.at(3);
        realJets[2] = vAllJets.at(4);
        realJets[3] = vAllJets.at(5);
        std::vector<ReconstructedParticle*> vRealJets(realJets, realJets + 4);
        vRealJets = sortGreaterPT(vRealJets);
        std::vector<TLorentzVector> Vj;
        for (int i = 0; i < 4; i++) {
            Vj.push_back(getTLorentzVector(vRealJets.at(i)));
        }
        fjet0_pT = Vj.at(0).Pt();
        fjet1_pT = Vj.at(1).Pt();
        fjet2_pT = Vj.at(2).Pt();
        fjet3_pT = Vj.at(3).Pt();
        fjet0_E = Vj.at(0).E();
        fjet1_E = Vj.at(1).E();
        fjet2_E = Vj.at(2).E();
        fjet3_E = Vj.at(3).E();


        _outputTree_event->Fill();
        //_outputTree_rech1->Fill();

    }


}

void LSSAna::end() {
    if (_outputTree_event)
    {
        tree_file = _outputTree_event->GetCurrentFile(); // just in case we switched to a new file.
        tree_file->Write();
    }

    delete tree_file;
}

double LSSAna::getInvMass(MCParticle* part1, MCParticle* part2) {
    TVector3 P12 = part1->getMomentum();
    P12 += part2->getMomentum();
    double E12 = part1->getEnergy() + part2->getEnergy();
    double M12 = sqrt( E12 * E12 - P12.Mag2() );

    return M12;
}

double LSSAna::getInvMass(ReconstructedParticle* part1, ReconstructedParticle* part2) {
    TVector3 P12 = part1->getMomentum();
    P12 += part2->getMomentum();
    double E12 = part1->getEnergy() + part2->getEnergy();
    double M12 = sqrt( E12 * E12 - P12.Mag2() );

    return M12;
}

TLorentzVector LSSAna::getTLorentzVector(MCParticle* part) {
    TLorentzVector Vpart;
    Vpart.SetPxPyPzE(   part->getMomentum()[0],
                        part->getMomentum()[1],
                        part->getMomentum()[2],
                        part->getEnergy());
    return Vpart;
}

TLorentzVector LSSAna::getTLorentzVector(ReconstructedParticle* part) {
    TLorentzVector Vpart;
    Vpart.SetPxPyPzE(   part->getMomentum()[0],
                        part->getMomentum()[1],
                        part->getMomentum()[2],
                        part->getEnergy());
    return Vpart;
}

bool LSSAna::lessDeltaRjl(ReconstructedParticle *part1, ReconstructedParticle *part2) {
    TLorentzVector Vjet1;
    TLorentzVector Vjet2;
    Vjet1.SetPxPyPzE(   part1->getMomentum()[0],
                        part1->getMomentum()[1],
                        part1->getMomentum()[2],
                        part1->getEnergy());
    Vjet2.SetPxPyPzE(   part2->getMomentum()[0],
                        part2->getMomentum()[1],
                        part2->getMomentum()[3],
                        part2->getEnergy());
    double DeltaRlj1 = std::min( Vjet1.DeltaR( Vl[0] ), Vjet1.DeltaR( Vl[1] ) );
    double DeltaRlj2 = std::min( Vjet2.DeltaR( Vl[0] ), Vjet2.DeltaR( Vl[1] ) );

    return ( DeltaRlj1 < DeltaRlj2 );
}

bool LSSAna::greaterPT(ReconstructedParticle *part1, ReconstructedParticle *part2) {
    TLorentzVector Vpart1;
    TLorentzVector Vpart2;
    Vpart1.SetPxPyPzE(  part1->getMomentum()[0],
                        part1->getMomentum()[1],
                        part1->getMomentum()[2],
                        part1->getEnergy());
    Vpart2.SetPxPyPzE(   part2->getMomentum()[0],
                        part2->getMomentum()[1],
                        part2->getMomentum()[3],
                        part2->getEnergy());
    return (Vpart1.Pt() > Vpart2.Pt());
}

std::vector<ReconstructedParticle*> LSSAna::sortLessDeltaRjl(std::vector<ReconstructedParticle*> &ivec) {
    const int vsize = ivec.size();
    ReconstructedParticle* temp;
    for(int i=1; i < vsize; i++) {
        for(int j = 0; j < vsize - i; j++) {
            if ( ! lessDeltaRjl( ivec[j], ivec[j+1] ) ) {
                temp = ivec[j];
                ivec[j] = ivec[j+1];
                ivec[j+1] = temp;
            }
        }
      }
    return ivec;
}

std::vector<ReconstructedParticle*> LSSAna::sortGreaterPT(std::vector<ReconstructedParticle*> &ivec) {
    const int vsize = ivec.size();
    ReconstructedParticle* temp;
    for(int i=1; i < vsize; i++) {
        for(int j = 0; j < vsize - i; j++) {
            if ( ! greaterPT( ivec[j], ivec[j+1] ) ) {
                temp = ivec[j];
                ivec[j] = ivec[j+1];
                ivec[j+1] = temp;
            }
        }
      }
    return ivec;
}
