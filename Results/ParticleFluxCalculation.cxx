#include "TreeReader.h"
#include <TCanvas.h>
#include <TH1F.h>
#include <TROOT.h>
#include <string>

// A derived class for specific analysis
class ParticleFluxCalculation : public TreeReader {
private:
  //************************************************************************************//
  //  Create the neceesary variables to read the Tree information
  //************************************************************************************//
  Int_t evt;
  Char_t ParticleName[20];
  Int_t ParentID;
  Int_t ParticleID;
  Int_t StepNumber;
  Double_t posParticle[3];
  Double_t KineticEnergy;
  Char_t InteractionType[20];
  Char_t TargetIsotope[20];
  Char_t CreatorProcessName[30];
  //************************************************************************************//

public:
  ParticleFluxCalculation(const char *filename, const char *treename = "tree")
      : TreeReader(filename, treename) {

    //************************************************************************************//
    // Initialize variables
    //************************************************************************************//
    evt = 0;
    ParticleName[20] = {0};
    ParentID = 0;
    ParticleID = 0;
    StepNumber = 0;
    posParticle[0] = 0.0;
    posParticle[1] = 0.0;
    posParticle[2] = 0.0;
    // KineticEnergy = 0.0;
    std::fill(std::begin(InteractionType), std::end(InteractionType), 0.0);
    std::fill(std::begin(TargetIsotope), std::end(TargetIsotope), 0.0);
    std::fill(std::begin(CreatorProcessName), std::end(CreatorProcessName),
              0.0);
    //************************************************************************************//
  }

  void SetBranchAddresses() override {
    //************************************************************************************//
    // Set branch addresses for your the variables
    //************************************************************************************//
    GetTree()->SetBranchAddress("fEvent", &evt);
    GetTree()->SetBranchAddress("fParticleName", &ParticleName);
    GetTree()->SetBranchAddress("fParentID", &ParentID);
    GetTree()->SetBranchAddress("fParticleID", &ParticleID);
    GetTree()->SetBranchAddress("fStepNumber", &StepNumber);
    GetTree()->SetBranchAddress("fX", &posParticle[0]);
    GetTree()->SetBranchAddress("fY", &posParticle[1]);
    GetTree()->SetBranchAddress("fZ", &posParticle[2]);
    GetTree()->SetBranchAddress("fKinEnergy", &KineticEnergy);
    GetTree()->SetBranchAddress("fInteractionType", &InteractionType);
    GetTree()->SetBranchAddress("targetIsotope", &TargetIsotope);
    GetTree()->SetBranchAddress("fCreatorProcessName", &CreatorProcessName);
    //************************************************************************************//
  }

  void ProcessEvent(Long64_t entry) override {
    //************************************************************************************//
    // Call base method to load the entry
    //************************************************************************************//
    TreeReader::ProcessEvent(entry);
    //************************************************************************************//
  }

  std::vector<std::string>
  GetUniqueElements(const std::vector<std::string> &vec) {
    std::set<std::string> uniqueSet(vec.begin(), vec.end());
    return std::vector<std::string>(uniqueSet.begin(), uniqueSet.end());
  }

  void Analyze(std::string histtitle, std::string histoutname,
               Int_t Savecanvas) {

    //************************************************************************************//
    // Call base Analyze method
    TreeReader::Analyze();
    //************************************************************************************//

    std::cout << "Starting custom analysis..." << std::endl;
    std::cout << std::endl;

    //************************************************************************************//
    // Create vectors to hold the variables in the tree
    std::vector<Int_t> v_evt;
    std::vector<std::string> v_fParticleName;
    std::vector<Int_t> v_fParentID;
    std::vector<Int_t> v_fParticleID;
    std::vector<Int_t> v_fStepNumber;
    std::vector<Double_t> v_posParticleX;
    std::vector<Double_t> v_posParticleY;
    std::vector<Double_t> v_posParticleZ;
    std::vector<Double_t> v_KineticEnergy;
    std::vector<std::string> v_interactionType;
    std::vector<std::string> v_targetIsotope;
    std::vector<std::string> v_fCreatorProcessName;
    //************************************************************************************//
    //************************************************************************************//
    // Loop through events and fill vectors
    for (Long64_t i = 0; i < GetEntries(); i++) {
      ProcessEvent(i);
      // std::cout<<"Processing entry " << i << std::endl;
      v_evt.push_back(evt);
      v_fParticleName.push_back(ParticleName);
      v_fParentID.push_back(ParentID);
      v_fParticleID.push_back(ParticleID);
      v_fStepNumber.push_back(StepNumber);
      v_posParticleX.push_back(posParticle[0]);
      v_posParticleY.push_back(posParticle[1]);
      v_posParticleZ.push_back(posParticle[2]);
      v_KineticEnergy.push_back(KineticEnergy);
      v_interactionType.push_back(InteractionType);
      v_targetIsotope.push_back(TargetIsotope);
      v_fCreatorProcessName.push_back(CreatorProcessName);
      //   std::cout<< i << " " << InteractionType << std::endl;
      //   std::cin.get();
    }
    //************************************************************************************//
    // Get the unique particle types
    //************************************************************************************//

    // Check the particle name vector to see which particles are present
    std::vector<std::string> uniqueParticles =
        GetUniqueElements(v_fParticleName);
    std::cout << "Total Number of Particles: " << v_fParticleName.size()
              << std::endl;

    // Print the unique particle types
    std::cout << "Particle types: " << std::endl;
    for (const auto &element : uniqueParticles) {
      std::cout << element << std::endl;
    }
    std::cout << std::endl;

    //************************************************************************************//
    // Create Histograms for particle flux and Fill
    //************************************************************************************//
    std::map<std::string, TH2F *> particleHistMap;
    std::map<std::string, TH1F *> particlePosMap;

    // std::string histtitle = "XY Plane +Z";
    for (const auto &element : uniqueParticles) {
      TString H_Title =
          TString::Format("%s Map for %s", element.c_str(), histtitle.c_str());
      TH2F *hist = new TH2F(element.c_str(), H_Title.Data(), 510, -510., 510.,
                            510, -510., 510.);

      TH1F *posZHist = new TH1F(
          TString::Format("Pos_%s", element.c_str()),
          TString::Format("Z Position Distribution for %s", element.c_str()),
          510, -510., 510.);

      particleHistMap[element] = hist;
      particlePosMap[element] = posZHist;
    } // end of creating histograms

    // Fill the histograms
    for (long unsigned int j = 0; j < v_evt.size(); j++) {
      std::string particle = v_fParticleName[j];
      Double_t x = v_posParticleX[j];
      Double_t y = v_posParticleY[j];
      Double_t z = v_posParticleZ[j];
      if (particleHistMap.find(particle) != particleHistMap.end() && z > 499.) {
        particleHistMap[particle]->Fill(x, y);
        // particlePosMap[particle]->Fill(v_posParticleZ[j]);
      }
      if (particlePosMap.find(particle) != particlePosMap.end() && z > 499.) {
        particlePosMap[particle]->Fill(v_posParticleZ[j]);
      }
    } // end of filling histograms

    //************************************************************************************//
    // Plot the histograms for particle positions along the direction of motion
    //************************************************************************************//
    for (const auto &pair : particlePosMap) {
      const std::string &particle = pair.first;
      TH1F *posZHist = pair.second;
      TCanvas *c2 =
          new TCanvas(TString::Format("Pos_%s", particle.c_str()),
                      TString::Format("Pos_%s", particle.c_str()), 1368, 1126);

      gStyle->SetOptStat(0);
      posZHist->GetXaxis()->SetTitle("Z (mm)");
      posZHist->GetXaxis()->SetTickLength(0.02);
      posZHist->GetXaxis()->SetLabelSize(0.03);
      posZHist->GetYaxis()->SetTitle("Counts");
      posZHist->GetYaxis()->SetTickLength(0.02);
      posZHist->GetYaxis()->SetLabelSize(0.03);
      posZHist->SetLineWidth(2);
      // Center the axis titles
      posZHist->GetXaxis()->CenterTitle(true);
      posZHist->GetYaxis()->CenterTitle(true);
      posZHist->Draw("HIST"); // draw with error bars

      // save the canvas
      TString Hist_out_name = TString::Format(
          "ParticlePosZ_%s_%s.png", particle.c_str(), histoutname.c_str());
      if (Savecanvas)
        c2->SaveAs(Hist_out_name.Data());

      // Cleanup
      delete c2;
      delete posZHist;
    }
    //************************************************************************************//
    // Plot the particle flux histograms and calculate the flux values and write
    // to output file
    //************************************************************************************//
    // Normalize the histogram
    // Neutron source strength is 2.8E10 n/s
    // Number of primaries created in the simulation is 10E7
    double scale = 2800.0;          // scale factor to get per second rate
    double area_cm2 = 100. * 100.0; // area in cm^2
    std::ofstream outFile;
    TString outputFileName =
        TString::Format("ParticleFlux_Results_%s.txt", histoutname.c_str());
    outFile.open(outputFileName.Data());
    if (!outFile.is_open()) {
      std::cerr << "Error opening output file: " << outputFileName.Data()
                << std::endl;
      exit(0);
    }

    std::cout << "Output file created: " << outputFileName.Data() << std::endl;
    outFile << "Particle Flux Results for " << histoutname.c_str() << ":\n";
    outFile << "----------------------------------------\n";

    for (const auto &pair : particleHistMap) {
      const std::string &particle = pair.first;
      TH2F *hist = pair.second;
      TCanvas *c = new TCanvas(particle.c_str(), particle.c_str(), 1368, 1126);

      Int_t totalEntries = hist->GetEntries();
      TString H_Title =
          TString::Format("%s Map for %s - \n # %s: %d", particle.c_str(),
                          histtitle.c_str(), particle.c_str(), totalEntries);

      hist->SetTitle(H_Title.Data());
      gStyle->SetOptStat(0);
      hist->GetXaxis()->SetTitle("X (mm)");
      hist->GetXaxis()->SetTickLength(0.02);
      hist->GetXaxis()->SetLabelSize(0.03);
      hist->GetYaxis()->SetTitle("Y (mm)");
      hist->GetYaxis()->SetTickLength(0.02);
      hist->GetYaxis()->SetLabelSize(0.03);
      hist->SetLineWidth(2);
      // Center the axis titles
      hist->GetXaxis()->CenterTitle(true);
      hist->GetYaxis()->CenterTitle(true);
      hist->Draw("COLZ"); // draw with color map
      // save the canvas
      TString Hist_out_name = TString::Format(
          "ParticleFlux_%s_%s.png", particle.c_str(), histoutname.c_str());
      if (Savecanvas)
        c->SaveAs(Hist_out_name.Data());

      // calculate the flux and save to output file
      double flux =
          totalEntries * scale / area_cm2; // flux in particles per second
      outFile << "Particle: " << particle << ", Total Counts: " << totalEntries
              << ", Flux (particles/s): " << flux << "\n";
      // Cleanup
      delete c;
      delete hist;
    }
    outFile.close();
    std::cout << "Custom analysis completed." << std::endl;
  } // End of Anlyze function
}; // enf od the class

// // Read Branches for checking
// void readBranchesExample() {
//   TreeReader
//   reader("NeutronSource_run0_thermalScattering_1E7_Primaries_Si_dep_"
//                     "DEDX_Calculation.root",
//                     "SiliconEdep_Z_1");
//   if (!reader.IsOpen())
//     return;

//   // Print branch information
//   reader.ListBranches();

//   // Read specific branch values for first few events
//   for (Long64_t i = 0; i < 5; i++) {
//     reader.GetTree()->GetEntry(i);

//     // Read different types of branches
//     Int_t eventId;
//     Int_t ParticleID;
//     Int_t ParentID;
//     // Char_t ParticleName[20];

//     if (reader.GetBranchValue("fEvent", eventId) &&
//         // reader.GetBranchValue("fParticleName", ParticleName) &&
//         reader.GetBranchValue("fParentID", ParentID) &&
//         reader.GetBranchValue("fParticleID", ParticleID)) {

//       std::cout << "Event " << i << ": ID=" << eventId
//                 // << ", ParticleName=" << ParticleName
//                 << ", ParentID=" << ParentID << ", ParticleID=" << ParticleID
//                 << std::endl;
//     }
//   }
// }

// Main function
int main() {
  //************************************************************************************//
  // Initialize configuration parameters
  //************************************************************************************//
  const char *config_file = "config_file_flux.cfg"; // name of the confog file
  TreeReader Config;               // create an object for TreeReader class
  Config.init_config(config_file); // initialize the configuration
  std::string RootFileName = Config.RootFileName; // get the root file name
  const char *TreeObjectName =
      Config.TreeObjectName;                // get the tree object name
  std::string HistTitle = Config.HistTitle; // get the histogram title
  std::string HistOutputName =
      Config.HistOutputName;            // get the histogram output name
  Int_t SaveCanvas = Config.SaveCanvas; // get the save canvas option
  std::cout << "Root File: " << RootFileName << std::endl;
  std::cout << "Tree Object Name: " << TreeObjectName << std::endl;
  std::cout << "HistTitle: " << HistTitle << std::endl;
  std::cout << "HistOutputName: " << HistOutputName << std::endl;
  //************************************************************************************//

  //************************************************************************************//
  // Create analysis object and run analysis
  //************************************************************************************//
  ParticleFluxCalculation analysis(RootFileName.c_str(), TreeObjectName);
  if (analysis.IsOpen()) {
    analysis.SetBranchAddresses();
    analysis.Analyze(HistTitle, HistOutputName,
                     SaveCanvas); // run analysis and pass the inputs
  }
  //************************************************************************************//

  //   // read branches for checking
  //   readBranchesExample();

  return 0;
}
