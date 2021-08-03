#include "LSSAna.hh"

LSSAna LSSAna_instance;

LSSAna::LSSAna() : Processor("LSSAna"),
                   foutput(0),
                   fmll(0),
                   fdelta_mll_mz(91.1876),
                   ffound_Zlepton(false) {
    _description = "LSS Ana";

    ftree_file_name = "LSSAna.root";
    registerProcessorParameter("TreeOutputFile",
        "The name of the file to which the ROOT tree will be written",
        ftree_file_name,
        ftree_file_name);

    fevent_tree_name = "event";
    registerProcessorParameter("ThreeName",
        "event level",
        fevent_tree_name,
        fevent_tree_name);

    frech1_tree_name = "rech1";
    registerProcessorParameter("rech1",
        "reconstructed h1",
        frech1_tree_name,
        frech1_tree_name);

    fcol_name = "RefinedJets";
    registerProcessorParameter( "CollectionName",
        "The name of Jet Collection.",
        fcol_name,
        fcol_name);

    foverwrite = 1;
    registerProcessorParameter( "OverwriteFile",
        "If zero an already existing file will not be overwriteen.",
        foverwrite,
        foverwrite);
    
    // init constants
    Leptons[0]=11; Leptons[1]=-11; Leptons[2]=13; Leptons[3]=-13; Leptons[4]=15; Leptons[5]=-15;
    fmz = 91.1876; // PDG2020
    fmh = 125.25;
}

LSSAna::~LSSAna() {

}

void LSSAna::init() {
    printParameters();

    ptree_file = new TFile(ftree_file_name.c_str(), (foverwrite ? "RECREATE" : "UPDATE") );

    if ( !ptree_file->IsOpen() ) {
        delete ptree_file;
        ptree_file = new TFile(ftree_file_name.c_str(), "NEW");
    }

    pevent_tree = new TTree( fevent_tree_name.c_str(), fevent_tree_name.c_str());
    prech1_tree = new TTree( frech1_tree_name.c_str(), frech1_tree_name.c_str());

    // Register Branch
    pevent_tree->Branch("EventNum", &feventNum, "EventNum/I");
    pevent_tree->Branch("Evisible", &fEvisible);
    pevent_tree->Branch("Emiss", &fEmiss);
    pevent_tree->Branch("Mll", &fmll);
    pevent_tree->Branch("deltaMllMZ", &fdelta_mll_mz);
    pevent_tree->Branch("Mrecoil", &fmrecoil);
    pevent_tree->Branch("deltaMrecoilMH", &fdelta_mrecoil_mh);
    // jets
    pevent_tree->Branch("Jet0pT", &fjet0_pT);
    pevent_tree->Branch("Jet1pT", &fjet1_pT);
    pevent_tree->Branch("Jet2pT", &fjet2_pT);
    pevent_tree->Branch("Jet3pT", &fjet3_pT);
    pevent_tree->Branch("Jet0E", &fjet0_E);
    pevent_tree->Branch("Jet1E", &fjet1_E);
    pevent_tree->Branch("Jet2E", &fjet2_E);
    pevent_tree->Branch("Jet3E", &fjet3_E);

}

void LSSAna::processEvent( LCEvent *evtP ) {
    if (evtP) {
        // declear containers
        IntVec iyth;
        DoubleVec vJpT;
        DoubleVec vJE;
        std::vector<ReconstructedParticle*> vAllJets; //sorted by DeltaRjl
        // read lcio data
        try {
            LCCollection* colJet = evtP->getCollection(fcol_name);
            LCCollection* MCPart = evtP->getCollection("MCParticle");
            NMCP = MCPart->getNumberOfElements();
            NJetsNum = colJet->getNumberOfElements();
            feventNum = evtP->getEventNumber();

            // Get MCParticle ( focus on Z leptons)
            for (int i = 0; i < NMCP; i++) {
                MCParticle *MCP_i = dynamic_cast<MCParticle*>(MCPart->getElementAt(i));
                int MCPID_i = MCP_i->getPDG();
                int *p = std::find(Leptons, Leptons + 6, MCPID_i);
                if(p == Leptons + 6) continue;

                for (int j = i + 1; j < NMCP; j++) {
                    MCParticle *MCP_j = dynamic_cast<MCParticle*>(MCPart->getElementAt(j));
                    int MCPID_j = MCP_j->getPDG();
                    if (MCPID_j != -MCPID_i ) continue; // same flavor differnt charge

                    double mll_temp = getInvMass(MCP_i, MCP_j);
                    if ( abs( mll_temp - fmz ) < fdelta_mll_mz ) {
                        fmll = mll_temp;
                        fdelta_mll_mz = abs( fmll - fmz );
                        Zleptons[0] = MCP_i;
                        Zleptons[1] = MCP_j;
                        ffound_Zlepton = true;
                    }
                }
            }
            if (!ffound_Zlepton) return;

            for (int i = 0; i < 2; i++){
                Vl[i].SetPxPyPzE(   Zleptons[i]->getMomentum()[0],
                                    Zleptons[i]->getMomentum()[1],
                                    Zleptons[i]->getMomentum()[2],
                                    Zleptons[i]->getEnergy());
            }


            // Get Jets
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
        // jets
        std::sort(vAllJets.begin(), vAllJets.end(), lessDeltaRjl );
        leptonJets[0] = vAllJets.at(0);
        leptonJets[1] = vAllJets.at(1);
        realJets[0] = vAllJets.at(2);
        realJets[1] = vAllJets.at(3);
        realJets[2] = vAllJets.at(4);
        realJets[3] = vAllJets.at(5);
        std::vector<ReconstructedParticle*> vRealJets(realJets, realJets + 4);
        std::sort(vRealJets.begin(), vRealJets.end(), greaterPT);

        
    }


}

void LSSAna::end() {

}

double LSSAna::getInvMass(MCParticle* part1, MCParticle* part2) {
    double E12 = part1->getEnergy() + part2->getEnergy();
    double P12x = part1->getMomentum()[0] + part2->getMomentum()[0];
    double P12y = part1->getMomentum()[1] + part2->getMomentum()[1];
    double P12z = part1->getMomentum()[2] + part2->getMomentum()[2];
    double P12sqr = P12x * P12x + P12y * P12y + P12z * P12z;
    double M12 = sqrt( E12 * E12 - P12sqr );

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