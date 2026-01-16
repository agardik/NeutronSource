#include "TreeReader.h"
#include <TCanvas.h>
#include <TH1F.h>
#include <TH3F.h>
#include <TFile.h>
#include <TROOT.h>
#include <TStyle.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <TPolyLine3D.h>
#include <TLegend.h>
#include <TGraph.h>


// A derived class for specific analysis
class DoseCalculation : public TreeReader {
private:
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

public:
  DoseCalculation(const char *filename, const char *treename = "tree")
      : TreeReader(filename, treename) {

    evt = 0;
    std::memset(ParticleName, 0, sizeof(ParticleName));
    ParentID = 0;
    ParticleID = 0;
    StepNumber = 0;
    posParticle[0] = posParticle[1] = posParticle[2] = 0.0;
    std::memset(InteractionType, 0, sizeof(InteractionType));
    std::memset(TargetIsotope, 0, sizeof(TargetIsotope));
    edepStep = 0;
    StopingTable = 0;
    StopingFull = 0;
    MeandEdX = 0;
    StopingPower = 0;
    std::memset(CreatorProcessName, 0, sizeof(CreatorProcessName));
  }

  void SetBranchAddresses() override {
    GetTree()->SetBranchAddress("fEvent", &evt);
    GetTree()->SetBranchAddress("fParticleName", ParticleName);
    GetTree()->SetBranchAddress("fParentID", &ParentID);
    GetTree()->SetBranchAddress("fParticleID", &ParticleID);
    GetTree()->SetBranchAddress("fStepNumber", &StepNumber);
    GetTree()->SetBranchAddress("fX", &posParticle[0]);
    GetTree()->SetBranchAddress("fY", &posParticle[1]);
    GetTree()->SetBranchAddress("fZ", &posParticle[2]);
    GetTree()->SetBranchAddress("fInteractionType", InteractionType);
    GetTree()->SetBranchAddress("targetIsotope", TargetIsotope);
    GetTree()->SetBranchAddress("Edep", &edepStep);
    GetTree()->SetBranchAddress("StopTable", &StopingTable);
    GetTree()->SetBranchAddress("StopFull", &StopingFull);
    GetTree()->SetBranchAddress("MeandEdx", &MeandEdX);
    GetTree()->SetBranchAddress("StopPower", &StopingPower);
    GetTree()->SetBranchAddress("fCreatorProcessName", CreatorProcessName);
  }

  void Analyze(std::string histtitle, std::string histoutname, Int_t Savecanvas) {

    // IMPORTANT: prevent ROOT ownership bugs
    TH1::AddDirectory(kFALSE);

    TString rootOutName =
        TString::Format("DoseResults_%s.root", histoutname.c_str());
    TFile *rootFile = new TFile(rootOutName, "RECREATE");

    TreeReader::Analyze();

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

    for (Long64_t i = 0; i < GetEntries(); i++) {
      ProcessEvent(i);
      v_evt.push_back(evt);
      v_fParticleName.emplace_back(ParticleName);
      v_fParentID.push_back(ParentID);
      v_fParticleID.push_back(ParticleID);
      v_fStepNumber.push_back(StepNumber);
      v_posParticleX.push_back(posParticle[0]);
      v_posParticleY.push_back(posParticle[1]);
      v_posParticleZ.push_back(posParticle[2]);
      v_interactionType.emplace_back(InteractionType);
      v_targetIsotope.emplace_back(TargetIsotope);
      v_edepStep.push_back(edepStep);
      v_StopTable.push_back(StopingTable);
      v_StopFull.push_back(StopingFull);
      v_MeandEdx.push_back(MeandEdX);
      v_StopPower.push_back(StopingPower);
      v_fCreatorProcessName.emplace_back(CreatorProcessName);
    }

    // Guard against empty tree
    if (v_evt.empty()) {
      std::cerr << "ERROR: Tree contains no entries!" << std::endl;
      rootFile->Close();
      delete rootFile;
      return;
    }
    // Create histograms
    TH1F *Hist_LET = new TH1F("alphadep", "Alpha Energy Deposition", 150, 0., 1.8);
    TH1F *Hist_TID = new TH1F("Li7dep", "Li7 Energy Deposition", 150, 0., 1.8);
    TH1F *Hist_StopTable =
        new TH1F("all_dep", "All Energy Deposition", 150, 0., 1.8);
    TH1F *Hist_StopFull =
        new TH1F("Hist_StopFull", "Alpha interaction X", 100, 0., 4.);
    TH1F *DEDX_Hist =
        new TH1F("DEDX_Hist", "dE/dx Alpha", 100, 0., 2.);
    TH1F *Hist_StopPower =
        new TH1F("Hist_StopPower", "Stopping Power", 100, 0., 4.);
    TH1F *Hist_AlphaAngle =
    new TH1F("AlphaAngle",
             "Alpha angular distribution;cos(#theta);Counts",
             50, -1.0, 1.0);

    
    Double_t edep_per_event = v_edepStep[0];
    Double_t StopPower_per_event = v_StopPower[0];
    bool have_first_step = false;
    bool have_first_alpha = true;
    double x0=0, y0=0, z0=0;
    // ---- Alpha track selection parameters ----
    const double ALPHA_EDEP_MIN = 1.24;   // MeV (adjust)
    const double ALPHA_EDEP_MAX = 1.4;   // MeV (adjust)
    const int MAX_ALPHA_TRACKS = 30;
    // ---- Storage for selected alpha tracks ----
    std::vector<std::vector<double>> alpha_x;
    std::vector<std::vector<double>> alpha_y;
    std::vector<std::vector<double>> alpha_z;

    alpha_x.reserve(MAX_ALPHA_TRACKS);
    alpha_y.reserve(MAX_ALPHA_TRACKS);
    alpha_z.reserve(MAX_ALPHA_TRACKS);

    // Temporary storage for current event
    std::vector<double> tmp_x, tmp_y, tmp_z;

    // Storage for all selected alpha tracks for track length edep calculation
    std::vector<double> alpha_track_lengths;
    std::vector<double> alpha_total_edep;





    int count1 = 0;
    // Loop over all entries
    for (size_t j = 1; j < v_evt.size(); j++) {
      if (v_evt[j] == v_evt[j - 1]) {
        // ---- Store alpha step positions ----
        if (v_fParticleName[j] == "alpha") {
          tmp_x.push_back(v_posParticleX[j]);
          tmp_y.push_back(v_posParticleY[j]);
          tmp_z.push_back(v_posParticleZ[j]);
        }

        // ----- Alpha angular distribution -----
        if (v_fParticleName[j] == "alpha" && have_first_alpha) {

          if (!have_first_step) {
            // store first step
            x0 = v_posParticleX[j];
            y0 = v_posParticleY[j];
            z0 = v_posParticleZ[j];
            have_first_step = true;
          }
          else {
            // second step → compute direction
            double dx = v_posParticleX[j] - x0;
            double dy = v_posParticleY[j] - y0;
            double dz = v_posParticleZ[j] - z0;
            have_first_step=false;
            double norm = std::sqrt(dx*dx + dy*dy + dz*dz);
            if (norm > 0) {
              double theta = std::acos(dx / norm) * 180.0 / M_PI;
              Hist_AlphaAngle->Fill(cos(theta*M_PI/180.0));
              // only once per alpha
              have_first_alpha = false;
            }

          }
        }
        edep_per_event += v_edepStep[j];
        StopPower_per_event += v_StopPower[j];
        if (v_fParticleName[j] == "alpha") {
          DEDX_Hist->Fill(v_MeandEdx[j]);
          Hist_StopFull->Fill(v_posParticleX[j]);
        }
      } else {
        have_first_alpha = true;
        if (v_fParticleName[j - 1] == "alpha")
          Hist_LET->Fill(edep_per_event);
              // Apply energy band cut
          if (edep_per_event >= ALPHA_EDEP_MIN &&
              edep_per_event <= ALPHA_EDEP_MAX &&
              (int)alpha_x.size() < MAX_ALPHA_TRACKS) {

            alpha_x.push_back(tmp_x);
            alpha_y.push_back(tmp_y);
            alpha_z.push_back(tmp_z);
          }
        if (v_fParticleName[j - 1] == "Li7")
          Hist_TID->Fill(edep_per_event);

        Hist_StopTable->Fill(edep_per_event);
        Hist_StopPower->Fill(StopPower_per_event);

        // --- Compute track length for previous alpha event ---
        if (v_fParticleName[j - 1] == "alpha" && tmp_x.size() > 1) {
          double trackLength = 0.0;
          for (size_t k = 1; k < tmp_x.size(); k++) {
            double dx = tmp_x[k] - tmp_x[k - 1];
            double dy = tmp_y[k] - tmp_y[k - 1];
            double dz = tmp_z[k] - tmp_z[k - 1];
            trackLength += std::sqrt(dx*dx + dy*dy + dz*dz);
          }
          alpha_track_lengths.push_back(trackLength);
          alpha_total_edep.push_back(edep_per_event);
        }

        if (v_fParticleName[j - 1] == "alpha" ||
            v_fParticleName[j - 1] == "Li7")
          count1++;

        edep_per_event = v_edepStep[j];
        StopPower_per_event = v_StopPower[j];
        tmp_x.clear();
        tmp_y.clear();
        tmp_z.clear();
      }
    }

    rootFile->cd();
    Hist_LET->Write();
    Hist_TID->Write();
    Hist_StopTable->Write();
    Hist_StopFull->Write();
    Hist_StopPower->Write();
    DEDX_Hist->Write();
    Hist_AlphaAngle->Write();
    // ---------- Plot selected alpha tracks 3d----------
    // ---------- Plot selected alpha tracks with axes ----------
    if (!alpha_x.empty()) {

      TCanvas *cTracks = new TCanvas("cAlphaTracks",
                                    "Alpha Tracks (Energy Selected)",
                                    800, 600);
      cTracks->cd();

      // ---- Determine plot bounds ----
      double xmin =  1e9, xmax = -1e9;
      double ymin =  1e9, ymax = -1e9;
      double zmin =  1e9, zmax = -1e9;

      for (size_t i = 0; i < alpha_x.size(); i++) {
        for (size_t p = 0; p < alpha_x[i].size(); p++) {
          xmin = std::min(xmin, alpha_x[i][p]);
          xmax = std::max(xmax, alpha_x[i][p]);
          ymin = std::min(ymin, alpha_y[i][p]);
          ymax = std::max(ymax, alpha_y[i][p]);
          zmin = std::min(zmin, alpha_z[i][p]);
          zmax = std::max(zmax, alpha_z[i][p]);
        }
      }

      // Add margins
      double dx = 0.1 * (xmax - xmin);
      double dy = 0.1 * (ymax - ymin);
      double dz = 0.1 * (zmax - zmin);

      // ---- 3D frame with axes ----
      TH3F *frame = new TH3F("frame",
                            "Alpha tracks;X;Y;Z",
                            10, xmin - dx, xmax + dx,
                            10, ymin - dy, ymax + dy,
                            10, zmin - dz, zmax + dz);

      frame->SetStats(0);
      frame->Draw("BOX");

      // ---- Draw tracks ----
      for (size_t i = 0; i < alpha_x.size(); i++) {
        int n = alpha_x[i].size();
        if (n < 2) continue;

        TPolyLine3D *pl = new TPolyLine3D(n);
        for (int p = 0; p < n; p++) {
          pl->SetPoint(p,
                      alpha_x[i][p],
                      alpha_y[i][p],
                      alpha_z[i][p]);
        }

        pl->SetLineWidth(2);
        pl->SetLineColor(1 + i % 9);
        pl->Draw("same");
      }

      cTracks->Write();
    }


    // ---------- XZ projection of selected alpha tracks ----------
    if (!alpha_x.empty()) {

      TCanvas *cTracksXZ = new TCanvas("cAlphaTracksXZ",
                                      "Alpha Tracks XZ Projection",
                                      800, 600);
      cTracksXZ->cd();

      TLegend *legXZ = new TLegend(0.65, 0.65, 0.88, 0.88);
      legXZ->SetBorderSize(0);

      bool first = true;

      for (size_t i = 0; i < alpha_x.size(); i++) {
        int n = alpha_x[i].size();
        if (n < 2) continue;

        TGraph *gr = new TGraph(n);
        for (int p = 0; p < n; p++) {
          gr->SetPoint(p,
                      alpha_x[i][p],  // X
                      alpha_z[i][p]); // Z
        }

        gr->SetLineWidth(2);
        gr->SetLineColor(1 + i % 9);
        gr->SetMarkerStyle(20);
        gr->SetMarkerSize(0.6);
        gr->SetMarkerColor(1 + i % 9);

        if (first) {
          gr->SetTitle("Alpha tracks XZ projection;X;Z");
          gr->Draw("AL");

          // ---- SET AXIS RANGES HERE ----
          gr->GetXaxis()->SetLimits(1.2, 3.1);  // X range
          gr->GetYaxis()->SetRangeUser(-16, 16); // Z range

          first = false;
        } else {
          gr->Draw("L same");
        }


        legXZ->AddEntry(gr,
                        Form("Alpha track %zu", i + 1),
                        "l");
      }

      legXZ->Draw();
      cTracksXZ->Write();
    }
    // -----------------------------------------------
    // ---------- Track Length vs Total Edep ----------
    if (!alpha_track_lengths.empty()) {
      TCanvas *cLengthEdep = new TCanvas("cLengthEdep",
                                        "Alpha Track Length vs Total Edep",
                                        800, 600);
      cLengthEdep->cd();

      // --- Determine axis ranges ---
      double minLength = *std::min_element(alpha_track_lengths.begin(), alpha_track_lengths.end());
      double maxLength = *std::max_element(alpha_track_lengths.begin(), alpha_track_lengths.end());
      double minEdep  = *std::min_element(alpha_total_edep.begin(), alpha_total_edep.end());
      double maxEdep  = *std::max_element(alpha_total_edep.begin(), alpha_total_edep.end());

      // Add 5-10% margin
      double xMargin = 0.05 * (maxLength - minLength);
      double yMargin = 0.05 * (maxEdep - minEdep);

      // --- Create graph ---
      int nPoints = alpha_track_lengths.size();
      TGraph *grLengthEdep = new TGraph(nPoints);
      grLengthEdep->SetTitle("Alpha Track Length vs Total Edep;Track length (mm);Total Edep (MeV)");
      grLengthEdep->SetMarkerStyle(20);
      grLengthEdep->SetMarkerSize(0.8);
      grLengthEdep->SetMarkerColor(kBlue);

      for (int i = 0; i < nPoints; i++) {
        grLengthEdep->SetPoint(i, alpha_track_lengths[i], alpha_total_edep[i]);
      }

      // --- Draw with explicit axis ranges ---
      grLengthEdep->GetXaxis()->SetLimits(minLength - xMargin, maxLength + xMargin);
      grLengthEdep->SetMinimum(minEdep - yMargin);
      grLengthEdep->SetMaximum(maxEdep + yMargin);

      grLengthEdep->Draw("AP"); // A=axes, P=points
      cLengthEdep->Write();
    }

    // -----------------------------------------------


    rootFile->Close();
    delete rootFile;

    std::cout << "Number of alpha and Li7 events: " << count1 << std::endl;

    // IMPORTANT: DO NOT delete histograms (ROOT + canvases still use them)
  }
};

int main() {
  const char *config_file = "config_file_dose.cfg";
  TreeReader Config;
  Config.init_config(config_file);

  DoseCalculation analysis(Config.RootFileName.c_str(),
                           Config.TreeObjectName);

  if (analysis.IsOpen()) {
    analysis.SetBranchAddresses();
    analysis.Analyze(Config.HistTitle,
                     Config.HistOutputName,
                     Config.SaveCanvas);
  }

  return 0;
}
