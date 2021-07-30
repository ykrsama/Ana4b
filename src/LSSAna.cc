#include "LSSAna.hh"

LSSAna LSSAna_instance;

LSSAna::LSSAna() 
: Processor("LSSAna"), 
_output(0) {
    _description = "LSS Ana";

    _treeFileName = "LSSAna.root";

}