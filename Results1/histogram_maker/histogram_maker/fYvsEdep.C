#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <iostream>
#include <fstream>

void plotFX_Edep_and_exportCSV() {

    // Open ROOT file
    TFile *f = new TFile("Electronics.root");
    if (!f || f->IsZombie()) {
        std::cout << "Error opening file!" << std::endl;
        return;
    }

    // Get tree
    TTree *t = (TTree*)f->Get("Electronics");
    if (!t) {
        std::cout << "Error: TTree not found!" << std::endl;
        return;
    }

    // Style
    gStyle->SetOptStat(0);

    // Canvas
    TCanvas *c1 = new TCanvas("c1", "Edep vs X", 900, 700);

    // Histogram (Y axis with energy weighting)
    TH1D *h1 = new TH1D("h1",
        "Energy Deposition vs X;X [mm];Deposited Energy [MeV]",
        100, -30, 230);

    // Fill histogram (Edep as weight)
    t->Draw("fX>>h1", "Edep", "goff");

    // Draw
    h1->Draw("hist");
    gPad->SetGrid();

    // Save plot
    c1->SaveAs("FX_Edep.png");

    // ---- Export histogram to CSV ----
    std::ofstream out("FX_Edep_histogram.csv");

    // Header
    out << "BinCenter_X[mm],Edep[MeV]\n";

    // Loop over bins
    for (int i = 1; i <= h1->GetNbinsX(); i++) {
        double y = h1->GetBinCenter(i);
        double edep = h1->GetBinContent(i);

        out << y << "," << edep << "\n";
    }

    out.close();

    std::cout << "CSV file 'FX_Edep_histogram.csv' created." << std::endl;

    f->Close();
}

plotFX_Edep_and_exportCSV();