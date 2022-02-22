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
    registerProcessorParameter("EventTreeName",
        "event level",
        _treeName_event,
        _treeName_event);

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

    //_outputTree_truth = new TTree( _treeName_truth.c_str(), _treeName_truth.c_str());
    _outputTree_event = new TTree( _treeName_event.c_str(), _treeName_event.c_str());

    // Register Branch
    _outputTree_event->Branch("EventNum", &feventNum, "EventNum/I");
    _outputTree_event->Branch("Evisible", &fEvisible);
    _outputTree_event->Branch("Emiss", &fEmiss);
    _outputTree_event->Branch("Mvisible", &fMass_visible);
    _outputTree_event->Branch("Mll", &fmll);
    _outputTree_event->Branch("Ell", &fEll);
    _outputTree_event->Branch("deltaMllMZ", &fdelta_mll_mz);
    _outputTree_event->Branch("Mrecoil", &fmrecoil);
    _outputTree_event->Branch("deltaMrecoilMH", &fdelta_mrecoil_mh);
    _outputTree_event->Branch("bNum", &fMCbNum);
    _outputTree_event->Branch("cNum", &fMCcNum);
    _outputTree_event->Branch("qNum", &fMCqNum);
    _outputTree_event->Branch("Eb", &fEb);
    // jets
    _outputTree_event->Branch("NJetsNum", &NJetsNum);
    _outputTree_event->Branch("PTJet", &fPTjet);
    _outputTree_event->Branch("EJet", &fEjet);
    _outputTree_event->Branch("DeltaRBJ", &fDeltaRBJ);
    _outputTree_event->Branch("EBJ", &fEBJ);
    _outputTree_event->Branch("Jet0E", &fjet0_E);
    _outputTree_event->Branch("Jet1E", &fjet1_E);
    _outputTree_event->Branch("Jet2E", &fjet2_E);
    _outputTree_event->Branch("Jet3E", &fjet3_E);
    _outputTree_event->Branch("Jet0pT", &fjet0_pT);
    _outputTree_event->Branch("Jet1pT", &fjet1_pT);
    _outputTree_event->Branch("Jet2pT", &fjet2_pT);
    _outputTree_event->Branch("Jet3pT", &fjet3_pT);
    // h1
    _outputTree_event->Branch("Mj00j01", &Mj00j01);
    _outputTree_event->Branch("Mj10j11", &Mj10j11);
    _outputTree_event->Branch("deltaMjj", &fdeltaMjj);
    _outputTree_event->Branch("DeltaRj00j01", &fDeltaRj00j01);
    _outputTree_event->Branch("DeltaRj10j11", &fDeltaRj10j11);
    _outputTree_event->Branch("JetTags", &fJetTags);
    _outputTree_event->Branch("lcfiplus0", &flcfiplus0);
    _outputTree_event->Branch("lcfiplus1", &flcfiplus1);
    _outputTree_event->Branch("lcfiplus2", &flcfiplus2);
    _outputTree_event->Branch("lcfiplus3", &flcfiplus3);
    //_outputTree_event->Branch("h1InvMass", &fh1InvMass);

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
        fMCbNum = 0;
        fMCcNum = 0;
        fMCqNum = 0;
        fbTagNum = 0;
        P_sum.SetXYZ(0,0,0);
        fJetTags.clear();
        vlcfiplus.clear();
        fEb.clear();
        fPTjet.clear();
        fEjet.clear();
        fDeltaRBJ.clear();
        fEBJ.clear();

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
                // get Parent info
                MCParticle *MCP_parent = nullptr;
                int MCPDG_parent = 0;
                if ( MCP_i->getParents().size() > 0 ) {
                    MCP_parent = MCP_i->getParents()[0];
                    MCPDG_parent = MCP_parent->getPDG();
                }
                // sum visible E
                if ( VTX.Mag() < 1 && fabs(MCPDG_i) != 12 && fabs(MCPDG_i) != 14 && fabs(MCPDG_i) != 16 ) {
                    if ( EndP.Mag() > 1) {
                        fEvisible += MCP_i->getEnergy();
                        P_sum += MCP_i->getMomentum();
                    }
                    if ( NJetsNum == 2 ) {
                        // for ee->ZH, H->bb ,find b by check its parents.
                        if( MCP_parent != nullptr &&
                            MCP_parent->getParents().size() > 0 && 
                            MCP_parent->getParents()[0]->getParents().size() > 0 &&
                            MCP_parent->getParents()[0]->getParents()[0]->getParents().size() == 0 ) {
                            if (fabs(MCPDG_i) >= 1 && fabs(MCPDG_i) <= 8) {
                            //if (fabs(MCPDG_i) == 5) {
                                MCbs[fMCbNum] = MCP_i;
                                fMCbNum++;
                            }
                            if (fabs(MCPDG_i) == 4) fMCcNum++;
                            if (fabs(MCPDG_i) >= 1 && fabs(MCPDG_i) <= 8) fMCqNum++;
                        }
                    } else if ( MCPDG_parent == 9999992 ) {
                        if (fabs(MCPDG_i) >= 1 && fabs(MCPDG_i) <= 8) {
                        //if (fabs(MCPDG_i) == 5) {
                            MCbs[fMCbNum] = MCP_i;
                            fMCbNum++;
                        }
                        if (fabs(MCPDG_i) == 4) fMCcNum++;
                        if (fabs(MCPDG_i) >= 1 && fabs(MCPDG_i) <= 8) fMCqNum++;
                    }
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
            if ( fMCbNum < 2 ) return;
            // if ( ! ffound_Zlepton ) {
            //     std::cout << "[INFO] Zleptons not found in event " << feventNum << ", skip" << std::endl;
            //     return;
            // }
            
            for (int i = 0; i < 2; i++){
                Vl[i].SetPxPyPzE(   Zleptons[i]->getMomentum()[0],
                                    Zleptons[i]->getMomentum()[1],
                                    Zleptons[i]->getMomentum()[2],
                                    Zleptons[i]->getEnergy());
            }
            //================================================================
            // Get Jets
            //-----------------------------------------------------------------
            PIDHandler pidH (colJet);
            // get algorithm IDs
            alcfiplus = pidH.getAlgorithmID("lcfiplus");
            ayth = pidH.getAlgorithmID("yth");
            // get lcfiplus parameter indiices
            ibtag = pidH.getParameterIndex(alcfiplus, "BTag");
            ictag = pidH.getParameterIndex(alcfiplus, "CTag");
            icategory = pidH.getParameterIndex(alcfiplus, "Category");
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
            // if ( vAllJets.size() < 6 ) { // sanity check
            //     std::cout << "[INFO] All Jets in event " << feventNum << " less than 6, skip" << std::endl;
            //     return;
            // }
            //================================================================
            // Process MCParticle
            //----------------------------------------------------------------
            fMass_visible = sqrt(fEvisible * fEvisible - P_sum.Mag2() );
            fEmiss = sqrt( S ) - fEvisible;
            fEll = Vl[0].Energy() + Vl[1].Energy();
            fmrecoil = sqrt( S - 2 * sqrt( S ) * fEll + fmll * fmll );
            fdelta_mrecoil_mh = fabs( fmrecoil - fmh );
            //================================================================
            // Process Jets
            //-----------------------------------------------------------------
            std::vector<ReconstructedParticle*> vRealJets;
            if ( NJetsNum >=4 ) {
                int leptonIndex = NJetsNum - 4;
                if ( NJetsNum > 4 && ffound_Zlepton )
                    vAllJets = sortLessDeltaRjl(vAllJets);
                for ( int j = 0; j < leptonIndex; j++)
                    leptonJets[j] = vAllJets.at(j);
                for ( int j = 0; j < 4; j++ )
                    vRealJets.push_back( vAllJets.at(j + leptonIndex) );
            } else {
                for ( int j = 0; j < NJetsNum; j++ )
                    vRealJets.push_back( vAllJets.at(j) );
            } 
            // sorted by greater E
            std::vector<TLorentzVector> Vj;
            std::sort(vRealJets.begin(), vRealJets.end(), greaterE);
            for (int i = 0; i < NJetsNum; i++) {
                Vj.push_back(getTLorentzVector(vRealJets.at(i)));
                fPTjet.push_back(Vj.at(i).Pt());
                fEjet.push_back(Vj.at(i).E());
            }
            if ( NJetsNum >= 4 ) {
                fjet0_pT = Vj.at(0).Pt();
                fjet1_pT = Vj.at(1).Pt();
                fjet2_pT = Vj.at(2).Pt();
                fjet3_pT = Vj.at(3).Pt();
                fjet0_E = Vj.at(0).E();
                fjet1_E = Vj.at(1).E();
                fjet2_E = Vj.at(2).E();
                fjet3_E = Vj.at(3).E();
            }
            // For Energy resolution
            //---------------------------------------------------------
            if ( fMCbNum == 4 || fMCbNum == 2) {
                std::vector<MCParticle*> vMCbs(MCbs, MCbs + fMCbNum);
                std::vector<TLorentzVector> VecB; // TLorentzVector of vMCbs
                std::sort(vMCbs.begin(), vMCbs.end(), greaterE);
                for (int i = 0; i < fMCbNum; i++)
                    VecB.push_back(getTLorentzVector(vMCbs.at(i)));
                //VecB = arrangeMCbByEBJ(VecB, Vj);
                for (int i = 0; i < fMCbNum; i++) {
                    fEb.push_back(VecB.at(i).E());
                    fDeltaRBJ.push_back(Vj.at(i).DeltaR(VecB.at(i)));
                    fEBJ.push_back(VecB.at(i).E() - Vj.at(i).E());
                }
            }
            if (NJetsNum >= 4 ) {
                // sorted by DeltaR
                vRealJets = arrangeJets( vRealJets );
                Vj.clear();
                for (int i = 0; i < 4; i++) {
                    Vj.push_back(getTLorentzVector(vRealJets.at(i)));
                    // get PID
                    const ParticleID &pid_lcfiplus = pidH.getParticleID(vRealJets.at(i), alcfiplus);
                    fbTag = pid_lcfiplus.getParameters()[ibtag];
                    fcTag = pid_lcfiplus.getParameters()[ictag];
                    fCategory = pid_lcfiplus.getParameters()[icategory];
                    flTag = 1 - fbTag - fcTag;
                    // find most likely jet tag
                    ftempTagParam = fbTag;
                    ftempTag = "b";
                    if ( fcTag > ftempTagParam )
                    {
                        ftempTagParam = fcTag;
                        ftempTag = "c";
                    }
                    if ( flTag > ftempTagParam )
                    {
                        ftempTagParam = flTag;
                        ftempTag = "l";
                    }
                    if (ftempTag == "b") fbTagNum++;
                    fJetTags.push_back(ftempTag);

                    // lcfiplus
                    lcfiplus_temp.clear();
                    lcfiplus_temp.push_back(fbTag);
                    lcfiplus_temp.push_back(fcTag);
                    lcfiplus_temp.push_back(fCategory);
                    vlcfiplus.push_back(lcfiplus_temp);
                }
                // sort vlcfiplus by bTag
                std::sort(vlcfiplus.begin(), vlcfiplus.end(), greaterBTag);
                flcfiplus0 = vlcfiplus.at(0);
                flcfiplus1 = vlcfiplus.at(1);
                flcfiplus2 = vlcfiplus.at(2);
                flcfiplus3 = vlcfiplus.at(3);

                Mj00j01 = getMassjj( vRealJets.at(0), vRealJets.at(1) );
                Mj10j11 = getMassjj( vRealJets.at(2), vRealJets.at(3) );
                fdeltaMjj = fabs(Mj00j01 - Mj10j11);
                fDeltaRj00j01 = Vj.at(0).DeltaR(Vj.at(1));
                fDeltaRj10j11 = Vj.at(2).DeltaR(Vj.at(3));
            }
        } catch (lcio::DataNotAvailableException err) {}
        //_outputTree_truth->Fill();
        _outputTree_event->Fill();
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

bool LSSAna::lessDeltaRjl(MCParticle *part1, ReconstructedParticle *part2) {
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

// bool LSSAna::greaterPT(ReconstructedParticle *part1, ReconstructedParticle *part2) {
//     TLorentzVector Vpart1 = getTLorentzVector(part1);
//     TLorentzVector Vpart2 = getTLorentzVector(part2);
//     return (Vpart1.Pt() > Vpart2.Pt());
// }

// bool LSSAna::greaterE(ReconstructedParticle *part1, ReconstructedParticle *part2) {
//     TLorentzVector Vpart1 = getTLorentzVector(part1);
//     TLorentzVector Vpart2 = getTLorentzVector(part2);
//     return (Vpart1.E() > Vpart2.E());
// }

// bool LSSAna::greaterE(MCParticle *part1, MCParticle *part2) {
//     TLorentzVector Vpart1 = getTLorentzVector(part1);
//     TLorentzVector Vpart2 = getTLorentzVector(part2);
//     return (Vpart1.E() > Vpart2.E());
// }

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

MCParticle* LSSAna::getCloserMCP(ReconstructedParticle *part1, std::vector<MCParticle*> &ivec) {
    TLorentzVector Vpart1 = getTLorentzVector(part1);
    const int vsize = ivec.size();
    MCParticle* temp;
    temp = ivec[0];
    for (int i = 1; i < vsize; i++) {
        TLorentzVector Vtemp = getTLorentzVector(temp);
        TLorentzVector Viveci = getTLorentzVector(ivec[i]);
        if ( Viveci.DeltaR(Vpart1) < Vtemp.DeltaR(Vpart1) ) {
            temp = ivec[i];
        }
    }
    return temp;
}

std::vector<TLorentzVector> LSSAna::arrangeMCbByEBJ(std::vector<TLorentzVector> &vecMCb, std::vector<TLorentzVector> &vecJet) {
    const int vsize = vecMCb.size();
    std::vector<TLorentzVector> ivec;
    float EBJ;
    float sd2;
    float sd2_temp;
    // for permutation
    std::vector<int> index;
    for (int i = 0; i < vsize; i++) index.push_back(i);
    // calculate standard deviation of EBJ
    sd2 = -1;
    do {
        // get standard deviation
        sd2_temp = 0;
        for ( int i = 0; i < vsize; i++) {
            EBJ = vecMCb.at( index.at(i) ).E() - vecJet.at( index.at(i) ).E();
            sd2_temp += EBJ * EBJ;
        }
        // add ivec
        if ( sd2_temp < sd2 || sd2 < 0 ) {
            sd2 = sd2_temp;
            ivec.clear();
            for ( int i = 0; i < vsize; i++ )
                ivec.push_back( vecMCb.at( index.at(i) ) );
        }
    } while ( std::next_permutation( index.begin(), index.end() ) );
    
    return ivec;
}

double LSSAna::getMassjj( ReconstructedParticle* part1, ReconstructedParticle* part2 ) {
    TVector3 Pjj = part1->getMomentum();
    Pjj += part2->getMomentum();
    double Ejj = part1->getEnergy() + part2->getEnergy();
    double Mjj = sqrt( Ejj * Ejj - Pjj.Mag2() );
    return Mjj;
}

std::vector<ReconstructedParticle*> LSSAna::arrangeJets(std::vector<ReconstructedParticle*> &ivec) {
    // arrange jets for reconstruction of h1. (ivec[0] ivec[1]) is a pair, (ivec[2] ivec[3]) is another pair.
    double Mjj[4][4];
    double deltaM, deltaM_temp;
    int jIndex[4];
    std::vector<ReconstructedParticle*> resVec;
    TLorentzVector Vj1, Vj2;
    for (int i = 0; i < 4; i++) {
        for (int j = i; j < 4; j++) {
            Mjj[i][j] = getMassjj(ivec.at(i), ivec.at(j));
        }
    }
    // 0 1, 2 3
    deltaM = fabs( Mjj[0][1] - Mjj[2][3] );
    jIndex[0] = 0;
    jIndex[1] = 1;
    jIndex[2] = 2;
    jIndex[3] = 3;
    // 0 2, 1 3
    deltaM_temp = fabs( Mjj[0][2] - Mjj[1][3] );
    if ( deltaM_temp < deltaM ) {
        deltaM = deltaM_temp;
        jIndex[0] = 0;
        jIndex[1] = 2;
        jIndex[2] = 1;
        jIndex[3] = 3;
    }
    // 0 3, 1 2
    deltaM_temp = fabs( Mjj[0][3] - Mjj[1][2] );
    if ( deltaM_temp < deltaM ) {
        deltaM = deltaM_temp;
        jIndex[0] = 0;
        jIndex[1] = 3;
        jIndex[2] = 1;
        jIndex[3] = 2;
    }

    // arrange
    for (int i = 0; i < 4; i++) {
        resVec.push_back( ivec.at( jIndex[i] ) );
    }
    return resVec;
}
