#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>

#include <set>
#include <map>
#include <string>
#include <iostream>
#include <fstream>

void foreignParticles() {

    // Open file
    TFile *f = TFile::Open("Electronics_Full.root");
    if (!f || f->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    // Get Electronics tree
    TTree *tree = (TTree*)f->Get("Electronics");
    if (!tree) {
        std::cerr << "Tree 'Electronics' not found!" << std::endl;
        return;
    }

    // Branches
    Int_t fEvent;
    Int_t fParticleID;
    Int_t fParentID;
    Char_t fParticleName[100];
    Char_t fPVatVertexname[100];

    tree->SetBranchAddress("fEvent", &fEvent);
    tree->SetBranchAddress("fParticleID", &fParticleID);
    tree->SetBranchAddress("fParentID", &fParentID);
    tree->SetBranchAddress("fParticleName", fParticleName);
    tree->SetBranchAddress("fPVatVertexname", fPVatVertexname);

    // Unique particles per event (considering parent ID to avoid counting multiple steps of the same particle)
    std::set<std::tuple<int,int,int>> uniqueParticles;

    // Count foreign particles by type 
    std::map<std::string,int> foreignCounts;

    Long64_t nEntries = tree->GetEntries();

    for (Long64_t i=0; i<nEntries; i++) {
        tree->GetEntry(i);

        std::tuple<int,int,int> key = std::make_tuple(fEvent, fParticleID, fParentID);

        if (uniqueParticles.find(key) != uniqueParticles.end())
            continue; // skip repeated steps

        uniqueParticles.insert(key);

        std::string pvname(fPVatVertexname);
        std::string pname(fParticleName);

        if (pvname != "LV_Electronics") {
            foreignCounts[pname]++;
        }
    }

    // Print foreign particles
    std::cout << "\n=== Foreign particles entering Electronics ===" << std::endl;
    for (auto &p : foreignCounts) {
        std::cout << p.first << " : " << p.second << std::endl;
    }

    // Save CSV
    std::ofstream outFile("foreign_particles.csv");
    outFile << "Particle,Count\n";
    for (auto &p : foreignCounts) {
        outFile << p.first << "," << p.second << "\n";
    }
    outFile.close();
    std::cout << "CSV saved: foreign_particles.csv" << std::endl;

    // Histogram
    int nTypes = foreignCounts.size();
    TH1F *h = new TH1F("hForeign","Foreign Particles;Particle Type;Count", nTypes, 0, nTypes);

    int bin = 1;
    for (auto &p : foreignCounts) {
        h->SetBinContent(bin, p.second);
        h->GetXaxis()->SetBinLabel(bin, p.first.c_str());
        bin++;
    }

    TCanvas *c = new TCanvas("cForeign","Foreign Particle Histogram", 900,600);
    h->LabelsOption("v");
    h->SetStats(0);
    h->Draw("hist");
    c->SaveAs("foreign_particles_hist.png");
}

foreignParticles()