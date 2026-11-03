#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TLegend.h>
#include <iostream>
#include <string>
#include <map>

void sphere_energy_spectra()
{
    // Open file
    TFile *f = TFile::Open("Centered_Sphere_1.root");
    if (!f || f->IsZombie()) {
        std::cout << "Error opening file!" << std::endl;
        return;
    }

    TTree *t = (TTree*)f->Get("Particles_Exit_World");
    if (!t) {
        std::cout << "Tree not found!" << std::endl;
        return;
    }

    // Variables
    Double_t fKinEnergy = 0;
    Char_t fParticleName[50];
    Char_t fPVatVertexname[100];

    t->SetBranchAddress("fKinEnergy", &fKinEnergy);
    t->SetBranchAddress("fParticleName", &fParticleName);
    t->SetBranchAddress("fPVatVertexname", &fPVatVertexname);

    // Map: volume name → histogram
    std::map<std::string, TH1D*> hMap;

    Long64_t nentries = t->GetEntries();

    for (Long64_t i = 0; i < nentries; i++)
    {
        t->GetEntry(i);

        std::string particle = fParticleName;

        if (particle != "gamma") continue;

        std::string volume = fPVatVertexname;

        // Create histogram if it doesn't exist
        if (hMap.find(volume) == hMap.end())
        {
            std::string hname = "h_" + volume;

            hMap[volume] = new TH1D(hname.c_str(),
                (volume + ";E_{kin} [MeV];Counts").c_str(),
                400, 0, 10);
        }

        hMap[volume]->Fill(fKinEnergy);
    }

    // Style
    gStyle->SetOptStat(0);

    // Canvas
    TCanvas *c1 = new TCanvas("c1", "Gamma Spectra per Volume", 1000, 700);
    c1->SetLogy();

    // Legend
    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);

    int color = 1;
    bool first = true;

    for (auto &it : hMap)
    {
        TH1D *h = it.second;

        h->SetLineColor(color++);
        h->SetLineWidth(2);

        if (first)
        {
            h->Draw();
            first = false;
        }
        else
        {
            h->Draw("SAME");
        }

        leg->AddEntry(h, it.first.c_str(), "l");

        // Avoid ROOT ugly colors overflow
        if (color == 10) color = 2;
    }

    leg->Draw();

    c1->SaveAs("gamma_spectra_by_volume.png");

    // Print summary
    for (auto &it : hMap)
    {
        std::cout << "Volume: " << it.first
                  << " | Entries: " << it.second->GetEntries()
                  << std::endl;
    }

    // =========================
// SAVE ALL HISTOGRAMS TO ONE CSV
// =========================
std::ofstream out("gamma_spectra_by_volume.csv");

// Header
out << "bin_center";

for (auto &it : hMap)
{
    out << "," << it.first;
}
out << "\n";

// Assume all histograms have same binning
TH1D* h_ref = hMap.begin()->second;
int nbins = h_ref->GetNbinsX();

for (int i = 1; i <= nbins; i++)
{
    double x = h_ref->GetBinCenter(i);
    out << x;

    for (auto &it : hMap)
    {
        TH1D *h = it.second;
        out << "," << h->GetBinContent(i);
    }

    out << "\n";
}

out.close();

std::cout << "Saved: gamma_spectra_by_volume.csv" << std::endl;
}