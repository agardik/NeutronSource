#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>

void radial_energy_density_with_dose()
{
    // ---- INPUT ----
    std::vector<std::string> files = {
        "Centered_Sphere_1.root",
        "Centered_Sphere_2.root",
        "Centered_Sphere_3.root",
        "Centered_Sphere_4.root",
        "Centered_Sphere_5.root"
    };

    std::vector<double> R = {300, 450, 600, 750, 900}; // mm

    // Silicon properties
    double thickness = 10.0; // mm
    double rho = 2.33e-3;    // g/mm^3
    double MeV_to_J = 1.602e-13;
    double g_to_kg = 1e-3;

    std::vector<double> energy_density;
    std::vector<double> gamma_density;
    std::vector<double> dose;
    std::vector<double> gamma_dose;

    for (size_t i = 0; i < files.size(); i++)
    {
        TFile *f = TFile::Open(files[i].c_str());

        if (!f || f->IsZombie()) {
            std::cout << "Error opening file: " << files[i] << std::endl;
            continue;
        }

        TTree *t = (TTree*)f->Get("Sphere");
        if (!t) {
            std::cout << "Tree 'Sphere' not found in " << files[i] << std::endl;
            f->Close();
            continue;
        }

        // Variables
        Double_t Edep;
        Char_t fParticleName[50];
        Char_t fCreatorProcessName[50]; 

        // Set branches
        t->SetBranchAddress("Edep", &Edep);
        t->SetBranchAddress("fParticleName", &fParticleName);
        t->SetBranchAddress("fCreatorProcessName", &fCreatorProcessName);

        double totalE = 0.0;
        double gammaE = 0.0;

        Long64_t nentries = t->GetEntries();

        for (Long64_t j = 0; j < nentries; j++)
        {
            t->GetEntry(j);

            std::string particle = fParticleName;
            std::string creator  = fCreatorProcessName;

            totalE += Edep;

            // Gamma + gamma-induced electrons
            if (particle == "gamma" ||
               (particle == "e-" && (creator == "compt" || creator == "phot")))
            {
                gammaE += Edep;
            }
        }

        // ---- Geometry ----
        double area = 4.0 * M_PI * R[i] * R[i]; // mm^2
        double volume = area * thickness;       // mm^3
        double mass = rho * volume * g_to_kg;   // kg

        // ---- Energy density ----
        energy_density.push_back(totalE / area);
        gamma_density.push_back(gammaE / area);

        // ---- Dose ----
        double totalDose_Gy = (totalE * MeV_to_J) / mass;
        double gammaDose_Gy = (gammaE * MeV_to_J) / mass;

        double totalDose_urad = totalDose_Gy * 1e8; // total dose in urad (assuming 1 Gy = 100 rad)
        double gammaDose_urad = gammaDose_Gy * 1e8; // gamma dose in urad (assuming 1 Gy = 100 rad)

        dose.push_back(totalDose_urad);
        gamma_dose.push_back(gammaDose_urad);

        std::cout << "R = " << R[i]
                  << " mm | Dose = " << totalDose_urad
                  << " urad | Gamma Dose = " << gammaDose_urad << std::endl;

        f->Close();
    }

    // ---- Plot ----
    gStyle->SetOptStat(0);

    TGraph *g_dose = new TGraph(R.size(), &R[0], &dose[0]);
    TGraph *g_gamma_dose = new TGraph(R.size(), &R[0], &gamma_dose[0]);

    g_dose->SetTitle("Dose vs Radius;Radius [mm];Dose [#murad]");
    g_dose->SetMarkerStyle(20);
    g_dose->SetLineWidth(2);

    g_gamma_dose->SetMarkerStyle(21);
    g_gamma_dose->SetLineStyle(2);

    TCanvas *c1 = new TCanvas("c1", "Dose Study", 900, 700);

    g_dose->Draw("APL");
    g_gamma_dose->Draw("PL SAME");

    gPad->SetGrid();
    c1->BuildLegend();

    c1->SaveAs("Dose_vs_Radius.png");

    // ---- CSV ----
    std::ofstream out("dose_results.csv");

    out << "Radius[mm],E/A[MeV/mm2],Gamma_E/A,Dose[urad],Gamma_Dose[urad]\n";

    for (size_t i = 0; i < R.size(); i++)
    {
        out << R[i] << ","
            << energy_density[i] << ","
            << gamma_density[i] << ","
            << dose[i] << ","
            << gamma_dose[i] << "\n";
    }

    out.close();

    std::cout << "Saved: dose_results.csv" << std::endl;
}