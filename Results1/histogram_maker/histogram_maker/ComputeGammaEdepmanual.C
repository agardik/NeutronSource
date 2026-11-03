{
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <iostream>
#include <cstring>
#include <string>

// =====================================================
// Approximate Silicon Mass Attenuation Coefficient
// Returns mu/rho in cm^2/g
// Energy input in MeV
// =====================================================

auto MassAttenuationSi = [](double E)
{
{
    if (E <= 0)
        return 0.0;

    if (E < 0.07)
    {
        return TMath::Power(
            10.0,
            (-2.8 * TMath::Log10(E) - 3.8)
        );
    }
    else if (E < 10)
    {
        return TMath::Power(
            10.0,
            (-0.42 * TMath::Log10(E) - 1.19)
        );
    }
    else 
    {
        return 0.025;
    }
};
};

// -------------------------------------------------
// Open input ROOT file
// -------------------------------------------------

TFile *fin = new TFile(
    "Electronics_Top_Strip_8_test_kill.root",
    "READ"
);

if (!fin || fin->IsZombie())
{
    std::cout << "ERROR: Cannot open input ROOT file" << std::endl;
    return;
}

// -------------------------------------------------
// Get Hits tree
// -------------------------------------------------

TTree *Hits = (TTree*)fin->Get("Hits");

if (!Hits)
{
    std::cout << "ERROR: Hits tree not found" << std::endl;
    return;
}

// -------------------------------------------------
// Input branches
// -------------------------------------------------

Int_t fEvent;
Double_t fKinEnergy;

Char_t fParticleName[64];
Char_t fVolumeName[128];

Hits->SetBranchAddress("fEvent", &fEvent);
Hits->SetBranchAddress("fKinEnergy", &fKinEnergy);
Hits->SetBranchAddress("fParticleName", fParticleName);
Hits->SetBranchAddress("fVolumeName", fVolumeName);

// -------------------------------------------------
// Create output ROOT file
// -------------------------------------------------

TFile *fout = new TFile(
    "GammaEdep_Silicon.root",
    "RECREATE"
);

TTree *outTree = new TTree(
    "GammaEdep",
    "Computed gamma energy deposition in silicon"
);

// -------------------------------------------------
// Output variables
// -------------------------------------------------

Double_t E_gamma;
Double_t mu_over_rho;
Double_t mu;
Double_t thickness_cm;
Double_t interactionProb;
Double_t EdepComputed;

// -------------------------------------------------
// Silicon properties
// -------------------------------------------------

const double rhoSi = 2.33; // g/cm^3

// 500 micron silicon detector
thickness_cm = 4; // cm

// -------------------------------------------------
// Output branches
// -------------------------------------------------

outTree->Branch("Event", &fEvent, "Event/I");
outTree->Branch("Egamma_MeV", &E_gamma, "Egamma_MeV/D");
outTree->Branch("MuOverRho", &mu_over_rho, "MuOverRho/D");
outTree->Branch("Mu", &mu, "Mu/D");
outTree->Branch("Thickness_cm", &thickness_cm, "Thickness_cm/D");
outTree->Branch(
    "InteractionProbability",
    &interactionProb,
    "InteractionProbability/D"
);
outTree->Branch(
    "ComputedEdep_MeV",
    &EdepComputed,
    "ComputedEdep_MeV/D"
);

// -------------------------------------------------
// Total energy deposition accumulator
// -------------------------------------------------

double TotalEdep = 0.0;

// -------------------------------------------------
// Event loop
// -------------------------------------------------

Long64_t nentries = Hits->GetEntries();

std::cout << "Processing "
          << nentries
          << " entries..."
          << std::endl;

for (Long64_t i = 0; i < nentries; i++)
{
    Hits->GetEntry(i);

    // Keep only gammas
    if (strcmp(fParticleName, "gamma") != 0)
        continue;

    // Keep only silicon volumes
    //std::string vol = fVolumeName;

    //if (vol.find("Silicon") == std::string::npos)
        //continue;

    // Incoming gamma energy
    E_gamma = fKinEnergy;

    if (E_gamma <= 0)
        continue;

    // Compute attenuation coefficient
    mu_over_rho = MassAttenuationSi(E_gamma);

    // Linear attenuation coefficient
    mu = mu_over_rho * rhoSi;

    // Interaction probability
    interactionProb =
        1.0 - TMath::Exp(-mu * thickness_cm);

    // Expected deposited energy
    EdepComputed =
        E_gamma * interactionProb;

    // Accumulate total deposited energy
    TotalEdep += EdepComputed;

    // Fill output tree
    outTree->Fill();
}

// -------------------------------------------------
// Write output file
// -------------------------------------------------

fout->cd();
outTree->Write();

fout->Close();
fin->Close();

// -------------------------------------------------
// Print results
// -------------------------------------------------

std::cout << std::endl;
std::cout << "======================================="
          << std::endl;

std::cout << "Gamma Energy Deposition Summary"
          << std::endl;

std::cout << "======================================="
          << std::endl;

std::cout << "Total Gamma Edep = "
          << TotalEdep
          << " MeV"
          << std::endl;

std::cout << "Output file = GammaEdep_Silicon.root"
          << std::endl;

std::cout << "======================================="
          << std::endl;
}