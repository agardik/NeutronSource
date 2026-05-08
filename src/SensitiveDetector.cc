#include "SensitiveDetector.hh"

#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4VProcess.hh"
#include "G4HadronicProcess.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4EmCalculator.hh"

// Allocator definition (MANDATORY)
G4ThreadLocal G4Allocator<MyHit>* MyHitAllocator = nullptr;

// ==========================
// Constructor
// ==========================
SensitiveDetector::SensitiveDetector(const G4String& name,
                                     const G4String& hitsCollectionName,
                                     G4bool excludeNeutrons)
    : G4VSensitiveDetector(name),
      fExcludeNeutrons(excludeNeutrons)
{
    collectionName.insert(hitsCollectionName);
}

// ==========================
// Initialize (called once per event)
// ==========================
void SensitiveDetector::Initialize(G4HCofThisEvent* hce)
{
    fHitsCollection = new MyHitsCollection(SensitiveDetectorName, collectionName[0]);

    if (fHCID < 0)
        fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);

    hce->AddHitsCollection(fHCID, fHitsCollection);

    // Cache event ID here (once per event) instead of calling
    // GetRunManager()->GetCurrentEvent() inside ProcessHits every step
    fCurrentEventID =
        G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
}

// ==========================
// ProcessHits (called every step inside the sensitive volume)
// ==========================
G4bool SensitiveDetector::ProcessHits(G4Step* step,
                                      G4TouchableHistory*)
{
    auto track = step->GetTrack();
    auto pre   = step->GetPreStepPoint();
    auto post  = step->GetPostStepPoint();

    if (!pre || !post) return false;

    // ======================
    // VOLUME FILTER
    // Both pre and post step must be inside this volume.
    // This replicates the ProcessVolume() condition from SteppingAction:
    //   if (thePostPVname != volumeName || thePrePVname != volumeName) return;
    // Note: boundary steps (entering/exiting) have pre XOR post inside —
    // they are correctly excluded here.
    // ======================
    G4VPhysicalVolume* prePV  = pre->GetPhysicalVolume();
    G4VPhysicalVolume* postPV = post->GetPhysicalVolume();

    //if (!prePV || !postPV) return false;

    const G4String& preName  = prePV->GetName();
    const G4String& postName = postPV->GetName();

    // The SD fires only for its own volume, but boundary steps can have
    // postPV pointing to the next volume. Reject those.
    //if (preName != postName) return false;

    // ======================
    // NEUTRON FILTER
    // Mirrors NeutronFlag == 1 logic in ProcessVolume().
    // Pass excludeNeutrons=true when registering the gas volume SD,
    // false for electronics/sphere.
    // ======================
    const G4String& particleName =
        track->GetDefinition()->GetParticleName();

    if (fExcludeNeutrons && particleName == "neutron") return false;

    // ======================
    // STEP LENGTH (cm)
    // ======================
    G4double stepLength = step->GetStepLength() / cm;
    if (stepLength <= 0.) return false;

    // ======================
    // ENERGY (MeV)
    // ======================
    G4double preE  = pre->GetKineticEnergy();
    G4double postE = post->GetKineticEnergy();
    G4double edep  = (preE - postE) / MeV;

    if (edep <= 0.) return false;

    // ======================
    // MATERIAL
    // ======================
    G4Material* material = post->GetMaterial();
    if (!material) return false;

    // Density in g/cm3 — needed for mass stopping power
    G4double density = material->GetDensity() / (g / cm3);

    // ======================
    // POSITION (mm) — store in mm, divide at fill time
    // ======================
    G4ThreeVector pos = post->GetPosition() / mm;

    // ======================
    // BASIC INFO
    // ======================
    G4int parentID  = track->GetParentID();
    G4int trackID   = track->GetTrackID();
    G4int stepNum   = track->GetCurrentStepNumber();

    // ======================
    // CREATOR PROCESS
    // ======================
    G4String creatorProcessName = "NoCreator";
    if (track->GetCreatorProcess())
        creatorProcessName = track->GetCreatorProcess()->GetProcessName();

    // ======================
    // VERTEX VOLUME
    // ======================
    const G4LogicalVolume* PVatVertex = track->GetLogicalVolumeAtVertex();
    G4String vertexName = PVatVertex ? PVatVertex->GetName() : "None";

    // ======================
    // INTERACTION TYPE & TARGET ISOTOPE
    // ======================
    G4String interactionType = "unknown";
    G4String targetIsotope   = "unknown";

    // Use const — no need for const_cast
    const G4VProcess* postProcess = post->GetProcessDefinedStep();

    if (postProcess)
    {
        G4HadronicProcess* hproc =
        dynamic_cast<G4HadronicProcess*>(const_cast<G4VProcess*>(postProcess));

        if (hproc && hproc->GetTargetIsotope())
        {
            interactionType = hproc->GetProcessName();
            targetIsotope   = hproc->GetTargetIsotope()->GetName();
        }
        else
        {
            interactionType = postProcess->GetProcessName();
            targetIsotope   = postPV ? postPV->GetName() : "None";
        }
    }

    // ======================
    // dE/dx — use G4EmCalculator (member, not local — constructed once)
    // Units: MeV/cm, then divide by density for mass stopping power MeV·cm²/g
    // This matches exactly what SteppingAction fills:
    //   stopTable / (CLHEP::MeV * CLHEP::cm2 / CLHEP::g)
    // ======================
    //G4double preE = pre->GetKineticEnergy();  // internal units

    G4double dEdxTable = 0.;
    G4double dEdxFull  = 0.;

    if (track->GetDefinition()->GetPDGCharge() != 0.)
    {
        // GetDEDX and ComputeTotalDEDX return Geant4 internal units (MeV/mm)
        // Divide by (MeV/cm) to get MeV/cm as a plain double
        dEdxTable = fEmCalc.GetDEDX(
            preE, track->GetDefinition(), material) / (MeV / cm);

        dEdxFull = fEmCalc.ComputeTotalDEDX(
            preE, track->GetDefinition(), material) / (MeV / cm);
    }

    // Mass stopping power: (MeV/cm) / (g/cm3) = MeV·cm²/g
    G4double stopTable = dEdxTable / density;   // MeV·cm²/g
    G4double stopFull  = dEdxFull  / density;   // MeV·cm²/g

    // Simulated stopping power from this step
    G4double meanDEDX  = edep / stepLength;     // MeV/cm
    G4double stopPower = meanDEDX / density;    // MeV·cm²/g

    // ======================
    // VOLUME NAME
    // ======================
    G4String volumeName = postPV ? postPV->GetName() : "OutOfWorld";

    // ======================
    // CREATE AND FILL HIT
    // ======================
    MyHit* hit = new MyHit();

    hit->SetEventID(fCurrentEventID);
    hit->SetParticleName(particleName);
    hit->SetParentID(parentID);
    hit->SetTrackID(trackID);
    hit->SetStepNumber(stepNum);

    hit->SetPosition(pos);               // mm

    hit->SetInteractionType(interactionType);
    hit->SetTargetIsotope(targetIsotope);
    hit->SetCreatorProcess(creatorProcessName);
    hit->SetVertexVolume(vertexName);

    hit->SetEdep(edep);                  // MeV
    hit->SetKineticEnergy(preE / MeV);   // MeV

    hit->SetStopTable(stopTable);        // MeV·cm²/g
    hit->SetStopFull(stopFull);          // MeV·cm²/g
    hit->SetMeanDEDX(meanDEDX);          // MeV/cm
    hit->SetStopPower(stopPower);        // MeV·cm²/g
    hit->SetVolumeName(volumeName);
    hit->SetSDName(SensitiveDetectorName);

    fHitsCollection->insert(hit);

    // ==================================
    // KILL TRACK WHEN ENTERING SD VOLUME
    // ==================================
    //step->GetTrack()->SetTrackStatus(fStopAndKill);

    return true;
}

// ==========================
// End of Event
// ==========================
void SensitiveDetector::EndOfEvent(G4HCofThisEvent*)
{
    // Data is read and written to ntuple in EventAction::EndOfEventAction
}