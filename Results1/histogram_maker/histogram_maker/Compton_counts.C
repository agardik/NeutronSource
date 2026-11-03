#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>

#include <set>
#include <tuple>
#include <string>
#include <iostream>

void countComptonElectrons() {

    TFile *f = TFile::Open("Centered_Sphere_1.root");
    if (!f || f->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    TTree *tree = (TTree*)f->Get("Sphere");
    if (!tree) {
        std::cerr << "Tree 'Sphere' not found!" << std::endl;
        return;
    }

    // Branches
    Int_t fEvent;
    Int_t fParticleID;
    Int_t fParentID;
    Char_t fParticleName[100];
    Char_t fCreatorProcess[100];

    tree->SetBranchAddress("fEvent", &fEvent);
    tree->SetBranchAddress("fParticleID", &fParticleID);
    tree->SetBranchAddress("fParentID", &fParentID);
    tree->SetBranchAddress("fParticleName", fParticleName);
    tree->SetBranchAddress("fCreatorProcessName", fCreatorProcess);

    // Unique particle tracking
    std::set<std::tuple<int,int,int>> uniqueParticles;

    int comptonElectrons = 0;

    Long64_t nEntries = tree->GetEntries();

    for (Long64_t i=0; i<nEntries; i++) {
        tree->GetEntry(i);

        std::tuple<int,int,int> key =
            std::make_tuple(fEvent, fParticleID, fParentID);

        if (uniqueParticles.find(key) != uniqueParticles.end())
            continue;

        uniqueParticles.insert(key);

        std::string pname(fParticleName);
        std::string creator(fCreatorProcess);

        // 🔥 THE IMPORTANT CONDITION
        if (pname == "e-" && creator == "compt") {
            comptonElectrons++;
        }
    }

    std::cout << "\n=== Compton electrons ===" << std::endl;
    std::cout << "Total: " << comptonElectrons << std::endl;

    // Histogram (optional: trivial single bin)
    TH1F *h = new TH1F("h","Compton Electrons",1,0,1);
    h->SetBinContent(1, comptonElectrons);

    TCanvas *c = new TCanvas("c","Compton electrons",600,400);
    h->Draw("hist");
    c->SaveAs("compton_electrons.png");
}

countComptonElectrons();