/*!  globals.h
 *   globals.h -- includes libraries, defines verbosity level of output stream.
 *	
 */
 
#ifndef _GLOBALS_H_      
#define _GLOBALS_H_   


//cpp
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>
#include <map>
#include <memory>
#include <ostream>
#include <set>

//ROOT 
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include "TSystem.h"
#include "TROOT.h"
#include "TEnv.h"
#include "TH2D.h"
#include <TH2.h>
#include "TMath.h"
#include "TBenchmark.h"
#include "TF1.h"
#include <TH1.h>
#include "TGraph.h"
#include "TMultiGraph.h"
#include <Math/MinimizerOptions.h>
#include "TCanvas.h"
#include "TDirectory.h"
#include <TStyle.h>

const int VERBOSITY_LEVEL = 4;

#endif

