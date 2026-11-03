#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <map>
#include <iostream>
#include <cstring>

void plot_total_edep()
{
    // Open ROOT file
    TFile *file = TFile::Open(
        "Electronics_Top_Strip_8_test_kill.root"
    );

    // Get tree
    TTree *tree = (TTree*)file->Get("Hits");

    // Variables
    Int_t fEvent;
    Double_t Edep;

    // Particle name
    Char_t fParticleName[64];

    // Threshold to ignore zero depositions
    double threshold = 1e-9;

    // Set branches
    tree->SetBranchAddress("fEvent", &fEvent);
    tree->SetBranchAddress("Edep", &Edep);
    tree->SetBranchAddress("fParticleName", fParticleName);

    // Map to accumulate total Edep per event
    std::map<int,double> eventEdep;

    // Running total
    double totalEdep = 0.0;

    Long64_t nentries = tree->GetEntries();

    for(Long64_t i = 0; i < nentries; i++)
    {
        tree->GetEntry(i);

        // Keep only gammas
        if(strcmp(fParticleName, "gamma") != 0)
            continue;

        if(Edep > threshold)
        {
            eventEdep[fEvent] += Edep;
            totalEdep += Edep;
        }
    }

    // Create histogram
    TH1D *hEdep = new TH1D(
        "hEdep",
        "Total Gamma Energy Deposition per Event;Total Edep;Counts",
        500,0,2
    );

    // Fill histogram
    for(auto const& ev : eventEdep)
    {
        if(ev.second > threshold)
            hEdep->Fill(ev.second);
    }

    // Draw histogram
    TCanvas *c1 = new TCanvas(
        "c1",
        "Total Gamma Edep",
        800,
        600
    );

    hEdep->Draw();

    c1->SaveAs("total_gamma_edep_histogram.png");

    // Print summary
    std::cout << "===================================="
              << std::endl;

    std::cout << "Gamma Energy Deposition Summary"
              << std::endl;

    std::cout << "===================================="
              << std::endl;

    std::cout << "Total gamma events processed: "
              << eventEdep.size()
              << std::endl;

    std::cout << "Total Gamma Edep in tree = "
              << totalEdep
              << std::endl;

    std::cout << "Histogram entries = "
              << hEdep->GetEntries()
              << std::endl;

    std::cout << "Histogram integral = "
              << hEdep->Integral()
              << std::endl;

    std::cout << "Histogram integral (width) = "
              << hEdep->Integral("width")
              << std::endl;

    std::cout << "===================================="
              << std::endl;
}