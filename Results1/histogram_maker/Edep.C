#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <map>
#include <iostream>

void plot_total_edep()
{
    // Open ROOT file
    TFile *file = TFile::Open("Electronics8.root");

    // Get tree
    TTree *tree = (TTree*)file->Get("Electronics");

    // Variables
    Int_t fEvent;
    Double_t Edep;

    // Threshold to ignore zero depositions
    double threshold = 1e-9;

    // Set branches
    tree->SetBranchAddress("fEvent", &fEvent);
    tree->SetBranchAddress("Edep", &Edep);

    // Map to accumulate total Edep per event
    std::map<int,double> eventEdep;


    // Running total
    double totalEdep = 0.0;

    Long64_t nentries = tree->GetEntries();

    for(Long64_t i = 0; i < nentries; i++)
    {
        tree->GetEntry(i);

        if(Edep > threshold)
        {
            eventEdep[fEvent] += Edep;
            totalEdep += Edep;
        }
    }

    // Create histogram
    TH1D *hEdep = new TH1D("hEdep",
                           "Total Energy Deposition per Event;Total Edep;Counts",
                           500,0,2);

    // Fill histogram
    for(auto const& ev : eventEdep)
    {
        if(ev.second > threshold)
            hEdep->Fill(ev.second);
    }

    // Draw histogram
    TCanvas *c1 = new TCanvas("c1","Total Edep",800,600);
    hEdep->Draw();

    c1->SaveAs("total_edep_histogram.png");

    std::cout << "Total events processed: " << eventEdep.size() << std::endl;

    std::cout << "Total Edep in tree = " << totalEdep << std::endl;
}