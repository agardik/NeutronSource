// XSectionBiasingOperator.hh
#include "G4VBiasingOperator.hh"
#include "G4BOptnChangeCrossSection.hh"

class XSectionBiasingOperator : public G4VBiasingOperator {
public:
  XSectionBiasingOperator(const G4String& name = "B4CCaptureBiasing");
  virtual ~XSectionBiasingOperator() = default;

  void SetBiasFactor(G4double f) { fBiasFactor = f; }

private:
  virtual G4VBiasingOperation*
  ProposeOccurenceBiasingOperation(const G4Track* track,
                                    const G4BiasingProcessInterface* callingProcess) override;

  // unused for pure occurrence biasing:
  virtual G4VBiasingOperation*
  ProposeFinalStateBiasingOperation(const G4Track*, const G4BiasingProcessInterface*) override
  { return nullptr; }
  virtual G4VBiasingOperation*
  ProposeNonPhysicsBiasingOperation(const G4Track*, const G4BiasingProcessInterface*) override
  { return nullptr; }

  G4BOptnChangeCrossSection* fXSChangeOperation;
  G4double fBiasFactor = 10.0; // default enhancement factor
};