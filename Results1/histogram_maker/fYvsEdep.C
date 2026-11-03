#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <iostream>
#include <fstream>

void plotFY_Edep_and_exportCSV() {

    // Open ROOT file
    TFile *f = new TFile("Electronics_Strip.root");
    if (!f || f->IsZombie()) {
        std::cout << "Error opening file!" << std::endl;
        return;
    }

    // Get tree
    TTree *t = (TTree*)f->Get("Hits");
    if (!t) {
        std::cout << "Error: TTree not found!" << std::endl;
        return;
    }

    // Style
    gStyle->SetOptStat(0);

    // Canvas
    TCanvas *c1 = new TCanvas("c1", "Edep vs Y", 900, 700);

    // Histogram (Y axis with energy weighting)
    TH1D *h1 = new TH1D("h1",
        "Energy Deposition vs Y;Y [mm];Deposited Energy [MeV]",
        100, -120, 0);

    // Fill histogram (Edep as weight)
    t->Draw("fY>>h1", "Edep", "goff");

    // Draw
    h1->Draw("hist");
    gPad->SetGrid();

    // Save plot
    c1->SaveAs("FY_Edep.png");

    // ---- Export histogram to CSV ----
    std::ofstream out("FY_Edep_histogram.csv");

    // Header
    out << "BinCenter_Y[mm],Edep[MeV]\n";

    // Loop over bins
    for (int i = 1; i <= h1->GetNbinsX(); i++) {
        double y = h1->GetBinCenter(i);
        double edep = h1->GetBinContent(i);

        out << y << "," << edep << "\n";
    }

    out.close();

    std::cout << "CSV file 'FY_Edep_histogram.csv' created." << std::endl;

    f->Close();
}

plotFY_Edep_and_exportCSV();