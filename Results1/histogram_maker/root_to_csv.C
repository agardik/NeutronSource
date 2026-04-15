#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <fstream>

void exportElectronicsToCSV() {

    // Open ROOT file
    TFile *f = new TFile("Electronics_Full.root");
    if (!f || f->IsZombie()) {
        std::cout << "Error opening file!" << std::endl;
        return;
    }

    // Get tree
    TTree *t = (TTree*)f->Get("Electronics");
    if (!t) {
        std::cout << "Error: TTree 'Electronics' not found!" << std::endl;
        return;
    }

    // Variables
    Int_t fEvent, fParentID, fParticleID, fStepNumber;
    Double_t fX, fY, fZ, Edep, MeandEdx, StopPower;
    Char_t fParticleName[50];
    Char_t fInteractionType[50];
    Char_t fCreatorProcessName[50];

    // Set branch addresses
    t->SetBranchAddress("fEvent", &fEvent);
    t->SetBranchAddress("fParentID", &fParentID);
    t->SetBranchAddress("fParticleID", &fParticleID);
    t->SetBranchAddress("fStepNumber", &fStepNumber);
    t->SetBranchAddress("fX", &fX);
    t->SetBranchAddress("fY", &fY);
    t->SetBranchAddress("fZ", &fZ);
    t->SetBranchAddress("Edep", &Edep);
    t->SetBranchAddress("MeandEdx", &MeandEdx);
    t->SetBranchAddress("StopPower", &StopPower);
    t->SetBranchAddress("fParticleName", &fParticleName);
    t->SetBranchAddress("fInteractionType", &fInteractionType);
    t->SetBranchAddress("fCreatorProcessName", &fCreatorProcessName);

    // Open CSV file
    std::ofstream out("Electronics.csv");

    // Header with units
    out << "Event,Particle,ParentID,ParticleID,Step,"
        << "X[mm],Y[mm],Z[mm],"
        << "Edep[MeV],Mean_dEdx[MeV/mm],StopPower[MeV/mm],"
        << "Interaction,CreatorProcess\n";

    // Loop over entries
    Long64_t nentries = t->GetEntries();

    for (Long64_t i = 0; i < nentries; i++) {
        t->GetEntry(i);

        out << fEvent << ","
            << fParticleName << ","
            << fParentID << ","
            << fParticleID << ","
            << fStepNumber << ","
            << fX << ","
            << fY << ","
            << fZ << ","
            << Edep << ","
            << MeandEdx << ","
            << StopPower << ","
            << fInteractionType << ","
            << fCreatorProcessName
            << "\n";
    }

    out.close();

    std::cout << "CSV file 'Electronics.csv' created with "
              << nentries << " entries." << std::endl;

    f->Close();
}

 exportElectronicsToCSV();