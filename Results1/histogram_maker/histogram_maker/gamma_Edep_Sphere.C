#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <iostream>
#include <string>

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

    // Branch variables
    Double_t fKinEnergy = 0;
    Double_t Edep = 0;
    Char_t fParticleName[50];

    t->SetBranchAddress("fKinEnergy", &fKinEnergy);
    t->SetBranchAddress("Edep", &Edep);
    t->SetBranchAddress("fParticleName", &fParticleName);

    // Histograms (log scale friendly)
    TH1D *h_gamma = new TH1D("h_gamma",
        "Gamma Kinetic Energy;E_{kin} [MeV];Counts",
        400, 0, 10);

    TH1D *h_neutron = new TH1D("h_neutron",
        "Neutron Kinetic Energy;E_{kin} [MeV];Counts",
        200, 0, 10);

    TH1D *h_gamma_edep = new TH1D("h_gamma_edep",
        "Gamma Energy Deposition;E_{dep} [MeV];Counts",
        200, 0, 10);

    TH1D *h_neutron_edep = new TH1D("h_neutron_edep",
        "Neutron Energy Deposition;E_{dep} [MeV];Counts",
        200, 0, 10);

    Long64_t nentries = t->GetEntries();

    for (Long64_t i = 0; i < nentries; i++)
    {
        t->GetEntry(i);

        std::string particle = fParticleName;

        if (particle == "gamma")
        {
            h_gamma->Fill(fKinEnergy);
            h_gamma_edep->Fill(Edep);
        }
        else if (particle == "neutron")
        {
            h_neutron->Fill(fKinEnergy);
            h_neutron_edep->Fill(Edep);
        }
    }

    // Style
    gStyle->SetOptStat(0);

    // Canvas
    TCanvas *c1 = new TCanvas("c1", "Kinetic Energy Spectra", 1000, 700);
    //c1->SetLogx();
    c1->SetLogy();

    h_gamma->SetLineColor(kRed);
    h_neutron->SetLineColor(kBlue);

    h_gamma->Draw();
    h_neutron->Draw("SAME");

    c1->BuildLegend();
    c1->SaveAs("kinetic_energy_gamma_neutron.png");

    // Second canvas for Edep
    TCanvas *c2 = new TCanvas("c2", "Energy Deposition", 1000, 700);
    c2->SetLogy();

    h_gamma_edep->SetLineColor(kRed);
    h_neutron_edep->SetLineColor(kBlue);

    h_gamma_edep->Draw();
    h_neutron_edep->Draw("SAME");

    c2->BuildLegend();
    c2->SaveAs("edep_gamma_neutron.png");

    // Print summary
    std::cout << "Gamma entries: " << h_gamma->GetEntries() << std::endl;
    std::cout << "Neutron entries: " << h_neutron->GetEntries() << std::endl;
}