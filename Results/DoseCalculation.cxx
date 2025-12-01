#include "TreeReader.h"
#include <TCanvas.h>
#include <TH1F.h>
#include <TROOT.h>
#include <iostream>

// A derived class for specific analysis
class DoseCalculation : public TreeReader {
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
  Char_t InteractionType[20];
  Char_t TargetIsotope[20];
  Double_t edepStep;
  Double_t StopingTable;
  Double_t StopingFull;
  Double_t MeandEdX;
  Double_t StopingPower;
  Char_t CreatorProcessName[30];
  Char_t PVatVertexname[30];
  //************************************************************************************//

public:
  DoseCalculation(const char *filename, const char *treename = "tree")
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
    std::fill(std::begin(InteractionType), std::end(InteractionType), 0.0);
    std::fill(std::begin(TargetIsotope), std::end(TargetIsotope), 0.0);
    edepStep = 0;
    StopingTable = 0;
    StopingFull = 0;
    MeandEdX = 0;
    StopingPower = 0;
    std::fill(std::begin(CreatorProcessName), std::end(CreatorProcessName),
              0.0); // Fills all elements with 0
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
    GetTree()->SetBranchAddress("fInteractionType", &InteractionType);
    GetTree()->SetBranchAddress("targetIsotope", &TargetIsotope);
    GetTree()->SetBranchAddress("Edep", &edepStep);
    GetTree()->SetBranchAddress("StopTable", &StopingTable);
    GetTree()->SetBranchAddress("StopFull", &StopingFull);
    GetTree()->SetBranchAddress("MeandEdx", &MeandEdX);
    GetTree()->SetBranchAddress("StopPower", &StopingPower);
    GetTree()->SetBranchAddress("fCreatorProcessName", &CreatorProcessName);
    GetTree()->SetBranchAddress("fPVatVertexname", &PVatVertexname);
    //************************************************************************************//
  }

  void ProcessEvent(Long64_t entry) override {
    //************************************************************************************//
    // Call base method to load the entry
    //************************************************************************************//
    TreeReader::ProcessEvent(entry);
    //************************************************************************************//
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
    std::vector<std::string> v_interactionType;
    std::vector<std::string> v_targetIsotope;
    std::vector<Double_t> v_edepStep;
    std::vector<Double_t> v_StopTable;
    std::vector<Double_t> v_StopFull;
    std::vector<Double_t> v_MeandEdx;
    std::vector<Double_t> v_StopPower;
    std::vector<std::string> v_fCreatorProcessName;
    std::vector<std::string> v_fPVatVertexname;
    //************************************************************************************//

    //************************************************************************************//
    // Loop through events and fill vectors
    for (Long64_t i = 0; i < GetEntries(); i++) {
      ProcessEvent(i);
      v_evt.push_back(evt);
      v_fParticleName.push_back(ParticleName);
      v_fParentID.push_back(ParentID);
      v_fParticleID.push_back(ParticleID);
      v_fStepNumber.push_back(StepNumber);
      v_posParticleX.push_back(posParticle[0]);
      v_posParticleY.push_back(posParticle[1]);
      v_posParticleZ.push_back(posParticle[2]);
      v_interactionType.push_back(InteractionType);
      v_targetIsotope.push_back(TargetIsotope);
      v_edepStep.push_back(edepStep);
      v_StopTable.push_back(StopingTable);
      v_StopFull.push_back(StopingFull);
      v_MeandEdx.push_back(MeandEdX);
      v_StopPower.push_back(StopingPower);
      v_fCreatorProcessName.push_back(CreatorProcessName);
      v_fPVatVertexname.push_back(PVatVertexname);
      //   std::cout<< i << " " << InteractionType << std::endl;
      //   std::cin.get();
    }
    //************************************************************************************//
    // Create the Histograms
    //************************************************************************************//
    TString H_Title_1 =
        TString::Format("LET Spectrum for %s (MeV / cm) ", histtitle.c_str());
    TH1F *Hist_LET = new TH1F("Hist_LET", H_Title_1.Data(), 400, 0., 4.);
    TString H_Title_2 =
        TString::Format("TID Spectrum for %s (rad) ", histtitle.c_str());
    TH1F *Hist_TID = new TH1F("Hist_TID", H_Title_2.Data(), 400, 0., 4.);
    TString H_Title_3 = TString::Format(
        "Stopping Power From Table Restricted for %s (MeV / cm / g)",
        histtitle.c_str());
    auto Hist_StopTable =
        new TH1F("Hist_StopTable", H_Title_3.Data(), 400, 0., 4.);

    TString H_Title_4 =
        TString::Format("Stopping Power From Table Total for %s (MeV / cm / g)",
                        histtitle.c_str());
    TH1F *Hist_StopFull =
        new TH1F("Hist_StopFull", H_Title_4.Data(), 400, 0., 4.);

    TString H_Title_5 = TString::Format(
        "dE/dX From Simulations for %s (MeV / cm ) ", histtitle.c_str());
    TH1F *DEDX_Hist = new TH1F("DEDX_Hist", H_Title_5.Data(), 400, 0., 4.);

    TString H_Title_6 = TString::Format(
        "Stopping Power From Simulations for %s (MeV / cm / g) ",
        histtitle.c_str());
    TH1F *Hist_StopPower =
        new TH1F("Hist_StopPower", H_Title_6.Data(), 400, 0., 4.);

    //************************************************************************************//
    // Calculate the energy deposition per event and divide by thickness to get
    // LET and TID and fill the histograms
    //************************************************************************************//
    Double_t edep_per_event = 0;
    edep_per_event = v_edepStep[0]; // initialize with the first entry
    Double_t StopTable_per_event = 0;
    StopTable_per_event = v_StopTable[0]; // initialize with the first entry
    Double_t StopFull_per_event = 0;
    StopFull_per_event = v_StopFull[0]; // initialize with the first entry
    Double_t DEDX_per_event = 0;
    DEDX_per_event = v_MeandEdx[0]; // initialize with the first entry
    Double_t StopPower_per_event = 0;
    StopPower_per_event = v_StopPower[0];     // initialize with the first entry
    Double_t Joule_conversion = 1.602E-13;    // From MeV to Joules
    Double_t area = 100 * 0.99 * 100. * 0.99; // in cm2, Area of the slab
    Double_t density = 2329.4;                // in mg/cm3, Silicon density
    Double_t mass =
        area * 0.2 * density * 1E-6; // in kg, 2 mm (0.2 cm) thickness,  Silicon
                                     // mass, 1E-6 for mg to g conversion

    for (long unsigned int j = 1; j < v_evt.size(); j++) {
      if (v_evt[j] == v_evt[j - 1]) {
        edep_per_event = edep_per_event + v_edepStep[j];
        StopTable_per_event = StopTable_per_event + v_StopTable[j];
        StopFull_per_event = StopFull_per_event + v_StopFull[j];
        StopPower_per_event = StopPower_per_event + v_StopPower[j];
        DEDX_per_event = DEDX_per_event + v_MeandEdx[j];
      } else {

        //   std::cout << v_evt[j] << "    " << edep_per_event << std::endl;
        Hist_LET->Fill(edep_per_event / 0.2); // in MeV/cm
        Hist_TID->Fill(edep_per_event, edep_per_event *
                                           (Joule_conversion / (mass)) *
                                           100.); // in rad
        edep_per_event = v_edepStep[j];

        Hist_StopTable->Fill(StopTable_per_event); // in MeV*cm2/g
        StopTable_per_event = v_StopTable[j];

        Hist_StopFull->Fill(StopFull_per_event); // in MeV*cm2/g
        StopFull_per_event = v_StopFull[j];

        Hist_StopPower->Fill(StopPower_per_event); // in MeV*cm2/g
        StopPower_per_event = v_StopPower[j];

        DEDX_Hist->Fill(DEDX_per_event); // in MeV/cm
        DEDX_per_event = v_MeandEdx[j];
      }
    } // End of the loop over all entries
    //************************************************************************************//

    //************************************************************************************//
    // Normalize the histogram
    // Neutron source strength is 2.8E10 n/s
    // Number of primaries created in the simulation is 10E7
    //************************************************************************************//
    double scale = 2800.0; // scale factor to get per second rate
    Hist_LET->Scale(scale);
    Hist_TID->Scale(scale * 1E6); // It is also scaled to urad
    Hist_StopTable->Scale(scale);
    Hist_StopFull->Scale(scale);
    Hist_StopPower->Scale(scale);
    DEDX_Hist->Scale(scale);
    //************************************************************************************//

    //************************************************************************************/
    // Draw the histogram for LET
    //************************************************************************************/
    TCanvas *c1 = new TCanvas("c1", "c1", 1368, 1126);
    gStyle->SetOptStat(0);
    Hist_LET->GetXaxis()->SetTitle("MeV/cm");
    Hist_LET->GetXaxis()->SetTickLength(0.02);
    Hist_LET->GetXaxis()->SetLabelSize(0.03);
    Hist_LET->GetYaxis()->SetTitle("# particles per second (Normalized)");
    //   Hist_LET->GetYaxis()->SetMoreLogLabels();
    Hist_LET->GetYaxis()->SetLabelSize(0.03);
    //   Hist_LET->GetYaxis()->SetMaxDigits(4);
    Hist_LET->SetLineWidth(2);
    gPad->SetLogy(1);
    //   Hist_LET->GetYaxis()->SetMoreLogLabels(false);
    // Center the axis titles
    Hist_LET->GetXaxis()->CenterTitle(true);
    Hist_LET->GetYaxis()->CenterTitle(true);
    Hist_LET->Draw("HISTE"); // draw with error bars

    // save the canvas
    TString Hist_out_name_1 =
        TString::Format("LET_Spectrum_%s.png", histoutname.c_str());
    if (Savecanvas)
      c1->SaveAs(Hist_out_name_1.Data());
    //************************************************************************************/
    // Draw the histograms for TID
    //************************************************************************************/
    TCanvas *c2 = new TCanvas("c2", "c2", 1368, 1126);
    gStyle->SetOptStat(0);
    Hist_TID->GetXaxis()->SetTitle("Energy (MeV)");
    Hist_TID->GetXaxis()->SetTickLength(0.02);
    Hist_TID->GetXaxis()->SetLabelSize(0.03);
    Hist_TID->GetYaxis()->SetTitle("TID (#murad)");
    //   Hist_TID->GetYaxis()->SetMoreLogLabels();
    Hist_TID->GetYaxis()->SetLabelSize(0.03);
    Hist_TID->GetYaxis()->SetMaxDigits(4);
    Hist_TID->SetLineWidth(2);
    // Center the axis titles
    Hist_TID->GetXaxis()->CenterTitle(true);
    Hist_TID->GetYaxis()->CenterTitle(true);
    Hist_TID->Draw("HISTE"); // draw with error bars

    // save the canvas
    TString Hist_out_name_2 =
        TString::Format("TID_Spectrum_%s.png", histoutname.c_str());
    if (Savecanvas)
      c2->SaveAs(Hist_out_name_2.Data());
    //************************************************************************************/
    // Draw the histograms for Stopping Power Restricted from Table by using
    // GetDEDX() method
    //************************************************************************************/
    TCanvas *c3 = new TCanvas("c3", "c3", 1368, 1126);
    gStyle->SetOptStat(0);
    Hist_StopTable->GetXaxis()->SetTitle("MeV * cm2 / g");
    Hist_StopTable->GetXaxis()->SetTickLength(0.02);
    Hist_StopTable->GetXaxis()->SetLabelSize(0.03);
    Hist_StopTable->GetYaxis()->SetTitle("# particles per second (Normalized)");
    //   Hist_StopTable->GetYaxis()->SetMoreLogLabels();
    Hist_StopTable->GetYaxis()->SetLabelSize(0.03);
    Hist_StopTable->GetYaxis()->SetMaxDigits(4);
    Hist_StopTable->SetLineWidth(2);
    // Center the axis titles
    Hist_StopTable->GetXaxis()->CenterTitle(true);
    Hist_StopTable->GetYaxis()->CenterTitle(true);
    Hist_StopTable->Draw("HISTE"); // draw with error bars

    // save the canvas
    TString Hist_out_name_3 =
        TString::Format("StoppingTable_Spectrum_%s.png", histoutname.c_str());
    if (Savecanvas)
      c3->SaveAs(Hist_out_name_3.Data());

    //************************************************************************************/
    // Draw the histograms for Stopping Power Unrestricted from Table by using
    // ComputeTotalDEDX() method
    //************************************************************************************/
    TCanvas *c4 = new TCanvas("c4", "c4", 1368, 1126);
    gStyle->SetOptStat(0);
    Hist_StopFull->GetXaxis()->SetTitle("MeV * cm2 / g");
    Hist_StopFull->GetXaxis()->SetTickLength(0.02);
    Hist_StopFull->GetXaxis()->SetLabelSize(0.03);
    Hist_StopFull->GetYaxis()->SetTitle("# particles per second (Normalized)");
    //   Hist_StopFull->GetYaxis()->SetMoreLogLabels();
    Hist_StopFull->GetYaxis()->SetLabelSize(0.03);
    Hist_StopFull->GetYaxis()->SetMaxDigits(4);
    Hist_StopFull->SetLineWidth(2);
    // Center the axis titles
    Hist_StopFull->GetXaxis()->CenterTitle(true);
    Hist_StopFull->GetYaxis()->CenterTitle(true);
    Hist_StopFull->Draw("HISTE"); // draw with error bars

    // save the canvas
    TString Hist_out_name_4 =
        TString::Format("StoppingFull_Spectrum_%s.png", histoutname.c_str());
    if (Savecanvas)
      c4->SaveAs(Hist_out_name_4.Data());
    //************************************************************************************/
    // Draw the histograms for Stopping Power from Simulations, dE/dX divided by
    // the material density
    //************************************************************************************/
    TCanvas *c5 = new TCanvas("c5", "c5", 1368, 1126);
    gStyle->SetOptStat(0);
    Hist_StopPower->GetXaxis()->SetTitle("MeV * cm2 / g");
    Hist_StopPower->GetXaxis()->SetTickLength(0.02);
    Hist_StopPower->GetXaxis()->SetLabelSize(0.03);
    Hist_StopPower->GetYaxis()->SetTitle("# particles per second (Normalized)");
    //   Hist_StopPower->GetYaxis()->SetMoreLogLabels();
    Hist_StopPower->GetYaxis()->SetLabelSize(0.03);
    Hist_StopPower->GetYaxis()->SetMaxDigits(4);
    Hist_StopPower->SetLineWidth(2);
    // Center the axis titles
    Hist_StopPower->GetXaxis()->CenterTitle(true);
    Hist_StopPower->GetYaxis()->CenterTitle(true);
    Hist_StopPower->Draw("HISTE"); // draw with error bars

    // save the canvas
    TString Hist_out_name_5 =
        TString::Format("StoppingPower_Spectrum_%s.png", histoutname.c_str());
    if (Savecanvas)
      c5->SaveAs(Hist_out_name_5.Data());

    //************************************************************************************/
    // Draw the histograms for dE/dX from Simulations, edep divided by the step
    // lenght
    //************************************************************************************//
    TCanvas *c6 = new TCanvas("c6", "c6", 1368, 1126);
    gStyle->SetOptStat(0);
    DEDX_Hist->GetXaxis()->SetTitle("MeV / cm");
    DEDX_Hist->GetXaxis()->SetTickLength(0.02);
    DEDX_Hist->GetXaxis()->SetLabelSize(0.03);
    DEDX_Hist->GetYaxis()->SetTitle("# particles per second (Normalized)");
    //   DEDX_Hist->GetYaxis()->SetMoreLogLabels();
    DEDX_Hist->GetYaxis()->SetLabelSize(0.03);
    DEDX_Hist->GetYaxis()->SetMaxDigits(4);
    DEDX_Hist->SetLineWidth(2);
    // Center the axis titles
    DEDX_Hist->GetXaxis()->CenterTitle(true);
    DEDX_Hist->GetYaxis()->CenterTitle(true);
    DEDX_Hist->Draw("HISTE"); // draw with error bars

    // save the canvas
    TString Hist_out_name_6 = TString::Format(
        "DEDX_with_steplenght_Spectrum_%s.png", histoutname.c_str());
    if (Savecanvas)
      c6->SaveAs(Hist_out_name_6.Data());
    //************************************************************************************/
    // Cleanup
    //************************************************************************************//
    delete Hist_LET;
    delete Hist_TID;
    delete Hist_StopTable;
    delete Hist_StopFull;
    delete DEDX_Hist;
    delete Hist_StopPower;
    delete c1;
    delete c2;
    delete c3;
    delete c4;
    delete c5;
    delete c6;
    //************************************************************************************/

    //************************************************************************************//
    // Calculate TID for each particle type
    // 0.2 for 2 mm thickness in cm
    // 100 for Gy to rad conversion
    // Conversions
    // 1 MeV = 1.602×10E-13 Joules
    //  1 mg = 10E-6 kg
    // J/kg = Gy
    // 1 Gy = 100 rad
    //************************************************************************************//
    Double_t TID_Total = 0.;
    Double_t TID_electron = 0.;
    Double_t TID_gamma = 0.;
    Double_t TID_proton = 0.;
    Double_t TID_nucleus = 0.;
    //************************************************************************************//
    // Total TID
    //************************************************************************************//
    for (long unsigned int j = 0; j < v_evt.size(); j++) {
      TID_Total =
          TID_Total + v_edepStep[j] * (Joule_conversion / (mass)) * 100.;
    }
    //************************************************************************************//
    // TID e+ end e+
    //************************************************************************************//
    // // Finding electron creator process
    // for (long unsigned int j = 0; j < v_evt.size(); j++) {
    //   if (v_fParticleName[j] == "e-" &&
    //       v_fCreatorProcessName[j] != "RadioactiveDecay" &&
    //       v_fCreatorProcessName[j] != "compt" &&
    //       v_fCreatorProcessName[j] != "phot" &&
    //       v_fCreatorProcessName[j] != "conv") {
    //     std::cout << v_fCreatorProcessName[j] << std::endl;
    //   }
    // }

    // For electrons, is it not compton and photo electric then it is comming
    // from the environment. Electrons depositing energy which are created due
    // to for instance compton scattering, are leaving energy deposition in the
    // slab. This electrons are secondary particles and this energy deposition
    // must  be counted for gammas. When gamma scattred by compton process, teh
    // energy deposition os not counted, it is zero. The energy deposition due
    // to this compton scattering is calculated thorugh the secondary electrons.
    // We would like to calculate the energy deposition of electron coming to
    // the slab from outside of the slab. Secondary electrons due to
    // comptonscattering are created inside the slab. So their deposition will
    // be counted as gamma deposition.

    for (long unsigned int j = 0; j < v_evt.size(); j++) {
      if (v_fParticleName[j] == "e-" && v_fCreatorProcessName[j] != "compt" &&
          v_fCreatorProcessName[j] != "phot") {
        TID_electron = TID_electron + v_edepStep[j] *
                                          (Joule_conversion / (mass)) *
                                          100.; // in rad
      }
      if (v_fParticleName[j] == "e+") {
        TID_electron = TID_electron + v_edepStep[j] *
                                          (Joule_conversion / (mass)) *
                                          100.; // in rad
      }
    }
    //************************************************************************************//
    // TID GAMMA
    //************************************************************************************//
    // // Finding gamma creator process
    // for (long unsigned int j = 0; j < v_evt.size(); j++) {
    //   if (v_fParticleName[j] == "gamma" &&
    //       v_fCreatorProcessName[j] != "nCapture" &&
    //       v_fCreatorProcessName[j] != "neutronInelastic" &&
    //       v_fCreatorProcessName[j] != "RadioactiveDecay") {
    //     std::cout << v_fCreatorProcessName[j] << std::endl;
    //   }
    // }

    for (long unsigned int j = 0; j < v_evt.size(); j++) {
      if (v_fParticleName[j] == "gamma") {
        TID_gamma =
            TID_gamma + v_edepStep[j] * (Joule_conversion / (mass)) * 100.;
      }
      if (v_fParticleName[j] == "e-" && (v_fCreatorProcessName[j] == "compt" ||
                                         v_fCreatorProcessName[j] == "phot")) {
        TID_gamma = TID_gamma + v_edepStep[j] * (Joule_conversion / (mass)) *
                                    100.; // in rad
      }
    }
    //************************************************************************************//
    // TID Proton and Nucleus
    //************************************************************************************//
    // // Finding proton creator process
    // for (long unsigned int j = 0; j < v_evt.size(); j++) {
    //   if (v_fParticleName[j] == "proton" &&
    //       v_fCreatorProcessName[j] != "neutronInelastic") {
    //     std::cout << v_fCreatorProcessName[j] << std::endl;
    //   }
    // }

    for (long unsigned int j = 0; j < v_evt.size(); j++) {
      if (v_fParticleName[j] == "proton") {
        TID_proton =
            TID_proton + v_edepStep[j] * (Joule_conversion / (mass)) * 100.;
      }
      if (v_fParticleName[j] == "Si28" || v_fParticleName[j] == "Si29" ||
          v_fParticleName[j] == "Si30" || v_fParticleName[j] == "Si31" ||
          v_fParticleName[j] == "P31") {
        TID_nucleus =
            TID_nucleus + v_edepStep[j] * (Joule_conversion / (mass)) * 100.;
      }
    }
    //************************************************************************************//
    // Print the results
    //************************************************************************************//
    std::cout << "Total TID in the " << histtitle.c_str()
              << " slab: " << TID_Total * scale * 1E6 << " urad" << std::endl;
    std::cout << "Total TID e-&e+ in the " << histtitle.c_str()
              << " slab: " << TID_electron * scale * 1E6 << " urad"
              << std::endl;
    std::cout << "Total TID gamma in the " << histtitle.c_str()
              << " slab: " << TID_gamma * scale * 1E6 << " urad" << std::endl;
    std::cout << "Total TID proton in the " << histtitle.c_str()
              << " slab: " << TID_proton * scale * 1E6 << " urad" << std::endl;
    std::cout << "Total TID nucleus in the " << histtitle.c_str()
              << " slab: " << TID_nucleus * scale * 1E6 << " urad" << std::endl;
    //************************************************************************************//

    //************************************************************************************//
    // Save the results to a text file
    //************************************************************************************//
    std::ofstream outFile;
    TString outputFileName =
        TString::Format("TID_Results_%s.txt", histoutname.c_str());
    outFile.open(outputFileName.Data());
    outFile << "TID Results for " << histtitle.c_str() << " slab\n";
    outFile << "----------------------------------------\n";
    outFile << "Total TID (urad): " << TID_Total * scale * 1E6 << "\n";
    outFile << "TID e- & e+ (urad): " << TID_electron * scale * 1E6 << "\n";
    outFile << "TID gamma (urad): " << TID_gamma * scale * 1E6 << "\n";
    outFile << "TID proton (urad): " << TID_proton * scale * 1E6 << "\n";
    outFile << "TID nucleus (urad): " << TID_nucleus * scale * 1E6 << "\n";
    outFile.close();
    //************************************************************************************//
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
  const char *config_file = "config_file_dose.cfg"; // name of the confog file
  TreeReader Config;               // create an object for TreeReader class
  Config.init_config(config_file); // initialize the configuration
  std::string RootFileName = Config.RootFileName; // get the root file name
  const char *TreeObjectName =
      Config.TreeObjectName;                // get the tree object name
  std::string HistTitle = Config.HistTitle; // get the histogram title
  std::string HistOutputName =
      Config.HistOutputName;            // get the histogram output name
  Int_t SaveCanvas = Config.SaveCanvas; // get the save canvas option
  // std::cout << "Root File: " << RootFileName << std::endl;
  // std::cout << "Tree Object Name: " << TreeObjectName << std::endl;
  // std::cout << "HistTitle: " << HistTitle << std::endl;
  //************************************************************************************//

  //************************************************************************************//
  // Create analysis object and run analysis
  //************************************************************************************//
  DoseCalculation analysis(RootFileName.c_str(), TreeObjectName);
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
