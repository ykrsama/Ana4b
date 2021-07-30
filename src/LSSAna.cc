#include "LSSAna.hh"

LSSAna LSSAna_instance;

LSSAna::LSSAna() : Processor("LSSAna"),
                   foutput(0) {
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

}

LSSAna::~LSSAna() {

}

void LSSAna::processEvent( LCEvent *evtP ) {

}

void LSSAna::end() {

}