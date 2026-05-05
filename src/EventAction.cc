//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file EventAction.cc
/// \brief Implementation of the EventAction class

#include "EventAction.hh"
#include "HistoManager.hh"
#include "Run.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event *anEvent) {

  //  Print the Run status
  G4int eventID = anEvent->GetEventID();
  G4Run *run = static_cast<G4Run *>(
      G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  G4int nOfEvents = run->GetNumberOfEventToBeProcessed();
  G4double evperCent = 10.; // status increment in percent

  if (fmod(eventID, double(nOfEvents * evperCent * 0.001)) == 0) {
    time_t my_time = time(NULL);
    tm *ltm = localtime(&my_time);
    G4double status = (100 * (eventID / double(nOfEvents)));
    std::cout << "=> Run " << eventID << " starts (" << status << "%, "
              << ltm->tm_hour << ":" << ltm->tm_min << ":" << ltm->tm_sec << ")"
              << std::endl;
  }

  fTotalEnergyDeposit = 0.;
  fTotalEnergyFlow = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::AddEdep(G4double Edep) { fTotalEnergyDeposit += Edep; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::AddEflow(G4double Eflow) { fTotalEnergyFlow += Eflow; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4SDManager.hh"
#include "SensitiveDetector.hh"

void EventAction::EndOfEventAction(const G4Event* event)
{
    Run* run = static_cast<Run*>(
        G4RunManager::GetRunManager()->GetNonConstCurrentRun());

    run->AddEdep(fTotalEnergyDeposit);
    run->AddEflow(fTotalEnergyFlow);

    auto analysisManager = G4AnalysisManager::Instance();

    analysisManager->FillH1(1, fTotalEnergyDeposit);
    analysisManager->FillH1(3, fTotalEnergyFlow);

    // ==========================
    // GET HITS COLLECTIONS
    // IDs are cached after the first event — avoids repeated SDManager
    // lookups. These are member variables (defined in EventAction.hh),
    // NOT static locals, to avoid the MT data race on first initialization.
    // ==========================
    auto* hce = event->GetHCofThisEvent();
    if (!hce) return;

    auto* sdMan = G4SDManager::GetSDMpointer();

    // Use G4HCtable::GetCollectionID(G4String) with "SDname/HCname" format.
    // This returns -1 silently if the collection was never registered
    // (i.e. the SD volume was not found in the GDML), unlike
    // G4SDManager::GetCollectionID() which prints a warning every event.
    auto* hcTable = sdMan->GetHCtable();

    if (fGasHCID < 0)
        fGasHCID = hcTable->GetCollectionID("GasSD/GasHitsCollection");
    if (fElectronicsHCID < 0)
        fElectronicsHCID = hcTable->GetCollectionID("ElectronicsSD/ElectronicsHitsCollection");
    if (fSphereHCID < 0)
        fSphereHCID = hcTable->GetCollectionID("SphereSD/SphereHitsCollection");

    // ==========================
    // LOOP OVER ALL THREE COLLECTIONS
    // All hits go into ntuple 9. Use hit->GetSDName() in ROOT to filter:
    //   tree->Draw("Edep", "SDName==\"GasSD\"")
    // ==========================
    for (G4int hcid : {fGasHCID, fElectronicsHCID, fSphereHCID})
    {
        if (hcid < 0) continue;  // collection not registered (SD not found)

        auto* hc = static_cast<MyHitsCollection*>(hce->GetHC(hcid));
        if (!hc) continue;       // volume not hit this event

        for (size_t i = 0; i < hc->GetSize(); i++)
        {
            auto* hit = (*hc)[i];

            analysisManager->FillNtupleIColumn(9, 0,  hit->GetEventID());
            analysisManager->FillNtupleSColumn(9, 1,  hit->GetParticleName());
            analysisManager->FillNtupleIColumn(9, 2,  hit->GetParentID());
            analysisManager->FillNtupleIColumn(9, 3,  hit->GetTrackID());
            analysisManager->FillNtupleIColumn(9, 4,  hit->GetStepNumber());

            auto pos = hit->GetPosition();
            analysisManager->FillNtupleDColumn(9, 5,  pos.x());
            analysisManager->FillNtupleDColumn(9, 6,  pos.y());
            analysisManager->FillNtupleDColumn(9, 7,  pos.z());

            analysisManager->FillNtupleDColumn(9, 8,  hit->GetEdep());
            analysisManager->FillNtupleDColumn(9, 9,  hit->GetKineticEnergy());

            analysisManager->FillNtupleSColumn(9, 10, hit->GetInteractionType());
            analysisManager->FillNtupleSColumn(9, 11, hit->GetTargetIsotope());
            analysisManager->FillNtupleSColumn(9, 12, hit->GetCreatorProcess());
            analysisManager->FillNtupleSColumn(9, 13, hit->GetVertexVolume());

            analysisManager->FillNtupleDColumn(9, 14, hit->GetStopTable());
            analysisManager->FillNtupleDColumn(9, 15, hit->GetStopFull());
            analysisManager->FillNtupleDColumn(9, 16, hit->GetMeanDEDX());
            analysisManager->FillNtupleDColumn(9, 17, hit->GetStopPower());
            analysisManager->FillNtupleSColumn(9, 18, hit->GetVolumeName());
            analysisManager->FillNtupleSColumn(9, 19, hit->GetSDName());

            analysisManager->AddNtupleRow(9);
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......