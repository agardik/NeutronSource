#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <iostream>
#include <map>
#include <string>
#include <fstream>
#include <cstring>

void particle_contributions()
{
    // -----------------------------
    // Open file
    // -----------------------------
    TFile *f = TFile::Open("testoutput2.root");   // new file

    if (!f || f->IsZombie()) {
        std::cout << "Error opening file!" << std::endl;
        return;
    }

    TTree *t = (TTree*)f->Get("Hits");

    if (!t) {
        std::cout << "Tree not found!" << std::endl;
        return;
    }

    // -----------------------------
    // Variables
    // (branches are null-terminated strings, so use zero-initialised buffers)
    // -----------------------------
    Double_t Edep = 0;
    Char_t fParticleName[100];
    Char_t fCreatorProcessName[100];
    Char_t fVolumeName[100];
    memset(fParticleName, 0, sizeof(fParticleName));
    memset(fCreatorProcessName, 0, sizeof(fCreatorProcessName));
    memset(fVolumeName, 0, sizeof(fVolumeName));

    t->SetBranchAddress("Edep", &Edep);
    t->SetBranchAddress("fParticleName", fParticleName);
    t->SetBranchAddress("fCreatorProcessName", fCreatorProcessName);
    t->SetBranchAddress("fVolumeName", fVolumeName);

    // Optional: only count hits in one volume (empty string = all hits)
    const std::string volumeFilter = "";   // e.g. "Electronics"

    // -----------------------------
    // Geometry (electronics slab), matching the XY map macro
    // -----------------------------
    double rho = 2.33e-3;   // g/mm^3 (Si)
    double thickness = 2.0; // mm

    int nx = 9, ny = 6;
    double xmin = 0,    xmax = 200;
    double ymin = -120, ymax = 0;

    double binArea = ((xmax - xmin)/nx) * ((ymax - ymin)/ny);
    double volume = binArea * thickness;
    double mass_kg = rho * volume * 1e-3;

    double MeV_to_J = 1.602e-13;
    double Gy_to_urad = 1e8;

    double doseFactor = (MeV_to_J / mass_kg) * Gy_to_urad;

    // -----------------------------
    // Accumulators
    // -----------------------------
    std::map<std::string, double> edep_map;
    std::map<std::string, double> dose_map;

    double gamma_total_edep = 0.0;
    double gamma_total_dose = 0.0;

    // -----------------------------
    // Loop over entries
    // -----------------------------
    Long64_t nentries = t->GetEntries();

    for (Long64_t i = 0; i < nentries; i++)
    {
        t->GetEntry(i);

        std::string particle = fParticleName;
        std::string creator  = fCreatorProcessName;

        if (!volumeFilter.empty() && std::string(fVolumeName) != volumeFilter)
            continue;

        // ---- Per particle ----
        edep_map[particle] += Edep;
        dose_map[particle] += Edep * doseFactor;

        // ---- Gamma-like contribution ----
        bool isGammaLike = false;

        // Direct gamma
        if (particle == "gamma")
            isGammaLike = true;

        // Gamma-induced electrons/positrons
        if ((particle == "e-" || particle == "e+") &&
            (creator == "compt" || creator == "phot" || creator == "conv"))
            isGammaLike = true;

        if (isGammaLike)
        {
            gamma_total_edep += Edep;
            gamma_total_dose += Edep * doseFactor;
        }
    }

    // Add combined gamma contribution
    edep_map["gamma_total"] = gamma_total_edep;
    dose_map["gamma_total"] = gamma_total_dose;

    // -----------------------------
    // Histograms
    // -----------------------------
    int n = edep_map.size();

    TH1D *hE = new TH1D("hE",
        "Total Energy Deposition per Particle;Particle Type;E_{dep} [MeV]",
        n, 0, n);

    TH1D *hD = new TH1D("hD",
        "Total Dose Contribution per Particle;Particle Type;Dose [#murad]",
        n, 0, n);

    int bin = 1;

    for (auto const& kv : edep_map)
    {
        const std::string &particle = kv.first;

        hE->GetXaxis()->SetBinLabel(bin, particle.c_str());
        hD->GetXaxis()->SetBinLabel(bin, particle.c_str());

        hE->SetBinContent(bin, kv.second);
        hD->SetBinContent(bin, dose_map[particle]);

        bin++;
    }

    // -----------------------------
    // Plot Energy
    // -----------------------------
    gStyle->SetOptStat(0);

    TCanvas *c1 = new TCanvas("c1","Edep per particle",900,700);
    c1->SetLogy();

    hE->LabelsOption("v");
    hE->Draw("HIST");

    c1->SaveAs("particle_edep.png");

    // -----------------------------
    // Plot Dose
    // -----------------------------
    TCanvas *c2 = new TCanvas("c2","Dose per particle",900,700);
    c2->SetLogy();

    hD->LabelsOption("v");
    hD->Draw("HIST");

    c2->SaveAs("particle_dose.png");

    // -----------------------------
    // CSV output
    // -----------------------------
    std::ofstream out("particle_contributions.csv");

    out << "Particle,Edep[MeV],Dose[urad]\n";

    for (auto const& kv : edep_map)
    {
        out << kv.first << ","
            << kv.second << ","
            << dose_map[kv.first] << "\n";
    }

    out.close();

    std::cout << "Saved: particle_contributions.csv" << std::endl;
}