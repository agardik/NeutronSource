#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

void radial_energy_density()
{
    // ---- USER INPUT ----
    std::vector<std::string> files = {
        "Centered_Sphere_1.root",
        "Centered_Sphere_2.root",
        "Centered_Sphere_3.root",
        "Centered_Sphere_4.root",
        "Centered_Sphere_5.root"
    };

    // Radii in mm (EDIT THESE!)
    std::vector<double> R = {300, 450, 600, 750, 900};

    // ---------------------

    std::vector<double> energy_density;
    std::vector<double> gamma_density;

    for (int i = 0; i < files.size(); i++)
    {
        TFile *f = TFile::Open(files[i].c_str());
        TTree *t = (TTree*)f->Get("Sphere");

        Double_t Edep;
        Char_t fParticleName[50];

        t->SetBranchAddress("Edep", &Edep);
        t->SetBranchAddress("fParticleName", &fParticleName);

        double totalE = 0.0;
        double gammaE = 0.0;

        Long64_t nentries = t->GetEntries();

        for (Long64_t j = 0; j < nentries; j++)
        {
            t->GetEntry(j);

            totalE += Edep;

            // Check for gamma
            if (std::string(fParticleName) == "gamma")
                gammaE += Edep;
        }

        // Sphere area (mm^2)
        double area = 4.0 * M_PI * R[i] * R[i];

        // Energy per unit area
        energy_density.push_back(totalE / area);
        gamma_density.push_back(gammaE / area);

        std::cout << "File: " << files[i]
                  << "  TotalE = " << totalE
                  << "  GammaE = " << gammaE << std::endl;

        f->Close();
    }

    // ---- Create graphs ----
    TGraph *g_total = new TGraph(R.size(), &R[0], &energy_density[0]);
    TGraph *g_gamma = new TGraph(R.size(), &R[0], &gamma_density[0]);

    g_total->SetTitle("Energy Deposition per Unit Area vs Radius;Radius [mm];E/A [MeV/mm^{2}]");
    g_total->SetMarkerStyle(20);

    g_gamma->SetMarkerStyle(21);

    // ---- Draw ----
    TCanvas *c1 = new TCanvas("c1", "Radial Study", 900, 700);

    g_total->Draw("APL");
    g_gamma->Draw("PL SAME");

    gPad->SetGrid();

    c1->BuildLegend();

    c1->SaveAs("Radial_Energy_Density.png");

    // ---- Save CSV ----
    std::ofstream out("radial_energy_density.csv");

    out << "Radius[mm],Total_EnergyDensity[MeV/mm^2],Gamma_EnergyDensity[MeV/mm^2]\n";

    for (int i = 0; i < R.size(); i++)
    {
        out << R[i] << ","
            << energy_density[i] << ","
            << gamma_density[i] << "\n";
    }

    out.close();

    std::cout << "CSV saved: radial_energy_density.csv" << std::endl;
}

 radial_energy_density();