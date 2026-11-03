#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>

#include <set>
#include <map>
#include <string>
#include <iostream>
#include <fstream>   // <-- for CSV

void countParticles() {

    // Open file
    TFile *f = TFile::Open("Centered_Sphere_1.root");
    if (!f || f->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    // Get tree
    TTree *tree = (TTree*)f->Get("Sphere");
    if (!tree) {
        std::cerr << "Tree 'Sphere' not found!" << std::endl;
        return;
    }

    // Variables
    Int_t fEvent;
    Int_t fParticleID;
    Char_t fParticleName[100];

    tree->SetBranchAddress("fEvent", &fEvent);
    tree->SetBranchAddress("fParticleID", &fParticleID);
    tree->SetBranchAddress("fParticleName", fParticleName);

    // Unique particles per event
    std::set<std::pair<int,int>> uniqueParticles;

    // Count per particle type
    std::map<std::string, int> particleCounts;

    Long64_t nEntries = tree->GetEntries();

    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);

        std::pair<int,int> key = std::make_pair(fEvent, fParticleID);

        if (uniqueParticles.find(key) != uniqueParticles.end())
            continue;

        uniqueParticles.insert(key);

        std::string pname = std::string(fParticleName);
        particleCounts[pname]++;
    }

    // =========================
    // WRITE CSV FILE
    // =========================
    std::ofstream outFile("particle_counts.csv");

    if (!outFile.is_open()) {
        std::cerr << "Error opening CSV file!" << std::endl;
        return;
    }

    // Header
    outFile << "Particle,Count\n";

    // Data
    for (auto &p : particleCounts) {
        outFile << p.first << "," << p.second << "\n";
    }

    outFile.close();

    std::cout << "CSV file written: particle_counts.csv" << std::endl;

    // =========================
    // OPTIONAL: Histogram
    // =========================
    int nTypes = particleCounts.size();

    TH1F *h = new TH1F("hParticles",
                       "Unique Particle Counts;Particle Type;Count",
                       nTypes, 0, nTypes);

    int bin = 1;
    for (auto &p : particleCounts) {
        h->SetBinContent(bin, p.second);
        h->GetXaxis()->SetBinLabel(bin, p.first.c_str());
        bin++;
    }

    TCanvas *c = new TCanvas("c", "Particle Histogram", 900, 600);
    h->LabelsOption("v");
    h->SetStats(0);
    h->Draw("hist");

    c->SaveAs("particle_histogram.png");
}

 countParticles()