#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TCanvas.h>

#include <set>
#include <tuple>
#include <string>
#include <iostream>
#include <sstream>

void gammaFlux2D() {

    const int NX = 9;
    const int NY = 6;

    // Create 2D histogram (bins correspond to grid positions)
    TH2F *h2 = new TH2F("h2","Gamma Flux Map;X index;Y index",
                       NX, 0, NX,
                       NY, 0, NY);

    // Loop over all files
    for (int ix = 0; ix < NX; ix++) {
        for (int iy = 0; iy < NY; iy++) {

            // Build filename
            std::stringstream fname;
            fname << "Electronics_x" << ix << "_y" << iy << ".root";

            TFile *f = TFile::Open(fname.str().c_str());
            if (!f || f->IsZombie()) {
                std::cerr << "Cannot open " << fname.str() << std::endl;
                continue;
            }

            TTree *tree = (TTree*)f->Get("Electronics");
            if (!tree) {
                std::cerr << "Tree missing in " << fname.str() << std::endl;
                f->Close();
                continue;
            }

            // Branches
            Int_t fEvent, fParticleID, fParentID;
            Char_t fParticleName[100];
            Char_t fPVatVertexname[100];

            tree->SetBranchAddress("fEvent", &fEvent);
            tree->SetBranchAddress("fParticleID", &fParticleID);
            tree->SetBranchAddress("fParentID", &fParentID);
            tree->SetBranchAddress("fParticleName", fParticleName);
            tree->SetBranchAddress("fPVatVertexname", fPVatVertexname);

            std::set<std::tuple<int,int,int>> uniqueParticles;

            int gammaCount = 0;

            Long64_t nEntries = tree->GetEntries();

            for (Long64_t i=0; i<nEntries; i++) {
                tree->GetEntry(i);

                std::tuple<int,int,int> key =
                    std::make_tuple(fEvent, fParticleID, fParentID);

                if (uniqueParticles.find(key) != uniqueParticles.end())
                    continue;

                uniqueParticles.insert(key);

                std::string pname(fParticleName);
                std::string pvname(fPVatVertexname);

                // Select ONLY gammas entering electronics
                if (pname == "gamma" && pvname != "LV_Electronics") {
                    gammaCount++;
                }
            }

            // Fill histogram
            h2->SetBinContent(ix+1, iy+1, gammaCount);

            f->Close();
        }
    }

    // Draw result
    TCanvas *c = new TCanvas("c2D","Gamma Flux 2D",900,700);
    h2->SetStats(0);
    h2->Draw("COLZ TEXT");

    c->SaveAs("gamma_flux_2D.png");
}

 gammaFlux2D();