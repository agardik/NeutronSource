#ifndef SensitiveDetector_hh
#define SensitiveDetector_hh

#include "G4VSensitiveDetector.hh"
#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4Allocator.hh"
#include "G4EmCalculator.hh"
#include "globals.hh"

// ==========================
// HIT CLASS (FULL DATA)
// ==========================
class MyHit : public G4VHit
{
public:
    MyHit() = default;
    virtual ~MyHit() = default;

    MyHit(const MyHit&) = default;
    MyHit& operator=(const MyHit&) = default;

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // -------- SETTERS --------
    void SetEventID(G4int v)                   { fEventID = v; }
    void SetParticleName(const G4String& v)    { fParticleName = v; }
    void SetParentID(G4int v)                  { fParentID = v; }
    void SetTrackID(G4int v)                   { fTrackID = v; }
    void SetStepNumber(G4int v)                { fStepNumber = v; }

    void SetPosition(const G4ThreeVector& v)   { fPos = v; }

    void SetInteractionType(const G4String& v) { fInteractionType = v; }
    void SetTargetIsotope(const G4String& v)   { fTargetIsotope = v; }
    void SetCreatorProcess(const G4String& v)  { fCreatorProcess = v; }
    void SetVertexVolume(const G4String& v)    { fVertexVolume = v; }

    void SetEdep(G4double v)                   { fEdep = v; }
    void SetKineticEnergy(G4double v)          { fKineticEnergy = v; }
    void SetStopTable(G4double v)              { fStopTable = v; }
    void SetStopFull(G4double v)               { fStopFull = v; }
    void SetMeanDEDX(G4double v)               { fMeanDEDX = v; }
    void SetStopPower(G4double v)              { fStopPower = v; }
    void SetVolumeName(const G4String& v)      { fVolumeName = v; }
    void SetSDName(const G4String& v)          { fSDName = v; }

    // -------- GETTERS --------

    // --- Identification ---
    G4int           GetEventID()         const { return fEventID; }
    const G4String& GetParticleName()    const { return fParticleName; }
    G4int           GetParentID()        const { return fParentID; }
    G4int           GetTrackID()         const { return fTrackID; }
    G4int           GetStepNumber()      const { return fStepNumber; }

    // --- Position ---
    G4ThreeVector   GetPosition()        const { return fPos; }

    // --- Physics / process ---
    const G4String& GetInteractionType() const { return fInteractionType; }
    const G4String& GetTargetIsotope()   const { return fTargetIsotope; }
    const G4String& GetCreatorProcess()  const { return fCreatorProcess; }
    const G4String& GetVertexVolume()    const { return fVertexVolume; }
    const G4String& GetVolumeName()      const { return fVolumeName; }
    const G4String& GetSDName()          const { return fSDName; }

    // --- Energy ---
    G4double        GetEdep()            const { return fEdep; }
    G4double        GetKineticEnergy()   const { return fKineticEnergy; }

    // --- Stopping power ---
    G4double        GetStopTable()       const { return fStopTable; }
    G4double        GetStopFull()        const { return fStopFull; }
    G4double        GetMeanDEDX()        const { return fMeanDEDX; }
    G4double        GetStopPower()       const { return fStopPower; }

private:

    // ===== IDENTIFICATION =====
    G4int    fEventID      = 0;
    G4String fParticleName = "";
    G4String fVolumeName   = "";
    G4String fSDName       = "";
    G4int    fParentID     = 0;
    G4int    fTrackID      = 0;
    G4int    fStepNumber   = 0;

    // ===== POSITION =====
    G4ThreeVector fPos;

    // ===== PHYSICS =====
    G4String fInteractionType = "";
    G4String fTargetIsotope   = "";
    G4String fCreatorProcess  = "";
    G4String fVertexVolume    = "";

    // ===== ENERGY =====
    G4double fEdep          = 0.;
    G4double fKineticEnergy = 0.;

    // ===== STOPPING POWER =====
    G4double fStopTable = 0.;
    G4double fStopFull  = 0.;
    G4double fMeanDEDX  = 0.;
    G4double fStopPower = 0.;
};

// ==========================
// ALLOCATOR (mandatory)
// ==========================
typedef G4THitsCollection<MyHit> MyHitsCollection;

extern G4ThreadLocal G4Allocator<MyHit>* MyHitAllocator;

inline void* MyHit::operator new(size_t)
{
    if (!MyHitAllocator)
        MyHitAllocator = new G4Allocator<MyHit>;
    return (void*)MyHitAllocator->MallocSingle();
}

inline void MyHit::operator delete(void* hit)
{
    MyHitAllocator->FreeSingle((MyHit*)hit);
}

// ==========================
// SENSITIVE DETECTOR
// ==========================
class SensitiveDetector : public G4VSensitiveDetector
{
public:
    // excludeNeutrons=true  → gas volume (mirrors NeutronFlag=1 in SteppingAction)
    // excludeNeutrons=false → electronics, sphere, etc.
    SensitiveDetector(const G4String& name,
                      const G4String& hitsCollectionName,
                      G4bool excludeNeutrons = false);

    virtual ~SensitiveDetector() = default;

    virtual void   Initialize(G4HCofThisEvent* hce) override;
    virtual G4bool ProcessHits(G4Step* step,
                               G4TouchableHistory* history) override;
    virtual void   EndOfEvent(G4HCofThisEvent* hce) override;

private:
    MyHitsCollection* fHitsCollection  = nullptr;
    G4int             fHCID            = -1;
    G4int             fCurrentEventID  = 0;    // cached once in Initialize()
    G4bool            fExcludeNeutrons = false;
    G4EmCalculator    fEmCalc;                 // constructed once, not per step
};

#endif
