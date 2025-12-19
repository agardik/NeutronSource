#include "TreeReader.h"
#include <TBranch.h>
#include <TDataType.h>
#include <TLeaf.h>
#include <TObjArray.h>

// Constructor
TreeReader::TreeReader()
    : m_file(nullptr), m_tree(nullptr), m_entries(0), m_isOpen(kFALSE) {}
// Contructor with parameters
TreeReader::TreeReader(const char *filename, const char *treename)
    : m_file(nullptr), m_tree(nullptr), m_entries(0), m_isOpen(kFALSE) {
  OpenFile(filename, treename);
}

// Destructor
TreeReader::~TreeReader() { CloseFile(); }

// Open ROOT file and get tree
Bool_t TreeReader::OpenFile(const char *filename, const char *treename) {
  CloseFile();

  m_filename = filename;
  m_treename = treename;

  // Open the file
  m_file = TFile::Open(m_filename, "READ");
  if (!m_file || m_file->IsZombie()) {
    std::cerr << "Error: Cannot open file " << m_filename << std::endl;
    m_isOpen = kFALSE;
    return kFALSE;
  }

  // Get the tree
  m_tree = (TTree *)m_file->Get(m_treename);
  if (!m_tree) {
    std::cerr << "Error: Cannot find tree " << m_treename << " in file "
              << m_filename << std::endl;
    m_file->Close();
    delete m_file;
    m_file = nullptr;
    m_isOpen = kFALSE;
    return kFALSE;
  }

  // Number of entries in the tree
  m_entries = m_tree->GetEntries();
  m_isOpen = kTRUE;

  // Initialize branch information
  InitializeBranches();

  std::cout << "Successfully opened file: " << m_filename << std::endl;
  std::cout << "Tree: " << m_treename << " with " << m_entries << " entries"
            << std::endl;
  std::cout << "Found " << GetNBranches() << " branches" << std::endl;

  return kTRUE;
}

// Initialize branch information
Bool_t TreeReader::InitializeBranches() {
  if (!m_isOpen || !m_tree)
    return kFALSE;

  m_branchNames.clear();
  m_branches.clear();
  m_branchTypes.clear();

  TObjArray *branches = m_tree->GetListOfBranches();
  if (!branches)
    return kFALSE;

  for (Int_t i = 0; i < branches->GetEntries(); i++) {
    TBranch *branch = (TBranch *)branches->At(i);
    TString branchName = branch->GetName();

    m_branchNames.push_back(branchName);
    m_branches[branchName] = branch;

    // Get branch type
    TLeaf *leaf = branch->GetLeaf(branchName);
    if (leaf) {
      m_branchTypes[branchName] = leaf->GetTypeName();
    } else {
      m_branchTypes[branchName] = "unknown";
    }
  }

  return kTRUE;
}

// Get specific branch
TBranch *TreeReader::GetBranch(const char *branchName) {
  TString name(branchName);
  if (m_branches.find(name) != m_branches.end()) {
    return m_branches[name];
  }
  return nullptr;
}

// Get branch type
TString TreeReader::GetBranchType(const char *branchName) {
  TString name(branchName);
  if (m_branchTypes.find(name) != m_branchTypes.end()) {
    return m_branchTypes[name];
  }
  return "not found";
}

// Close the file
void TreeReader::CloseFile() {
  if (m_file) {
    m_file->Close();
    delete m_file;
    m_file = nullptr;
  }
  m_tree = nullptr;
  m_entries = 0;
  m_isOpen = kFALSE;
  m_branchNames.clear();
  m_branches.clear();
  m_branchTypes.clear();
}

// Print basic information
void TreeReader::PrintInfo() const {
  if (!m_isOpen) {
    std::cout << "No file open." << std::endl;
    return;
  }

  std::cout << "\n=== File Information ===" << std::endl;
  std::cout << "File: " << m_filename << std::endl;
  std::cout << "Tree: " << m_treename << std::endl;
  std::cout << "Entries: " << m_entries << std::endl;
  std::cout << "Branches: " << GetNBranches() << std::endl;

  if (m_tree) {
    std::cout << "Total size: " << m_tree->GetTotBytes() << " bytes"
              << std::endl;
    std::cout << "Zipped size: " << m_tree->GetZipBytes() << " bytes"
              << std::endl;
  }
}

// List all branches in the tree
void TreeReader::ListBranches() const {
  if (!m_isOpen || !m_tree) {
    std::cout << "No tree available." << std::endl;
    return;
  }

  std::cout << "\n=== Branches in tree ===" << std::endl;
  for (size_t i = 0; i < m_branchNames.size(); i++) {
    TString name = m_branchNames[i];
    TString type = m_branchTypes.at(name);
    std::cout << "Branch " << i << ": " << name << " [" << type << "]"
              << std::endl;
  }
}

// Print detailed information about a specific branch
void TreeReader::PrintBranchInfo(const char *branchName) const {
  if (!m_isOpen || !m_tree) {
    std::cout << "No tree available." << std::endl;
    return;
  }

  TString name(branchName);
  if (m_branches.find(name) == m_branches.end()) {
    std::cout << "Branch '" << branchName << "' not found." << std::endl;
    return;
  }

  TBranch *branch = m_branches.at(name);
  std::cout << "\n=== Branch Information ===" << std::endl;
  std::cout << "Name: " << branch->GetName() << std::endl;
  std::cout << "Title: " << branch->GetTitle() << std::endl;
  std::cout << "Type: " << m_branchTypes.at(name) << std::endl;
  std::cout << "Entries: " << branch->GetEntries() << std::endl;
  std::cout << "Total size: " << branch->GetTotBytes() << " bytes" << std::endl;
  std::cout << "Zipped size: " << branch->GetZipBytes() << " bytes"
            << std::endl;
}

// Branch value reading methods
Bool_t TreeReader::GetBranchValue(const char *branchName, Int_t &value) {
  return GetBranchValue<Int_t>(branchName, value);
}

Bool_t TreeReader::GetBranchValue(const char *branchName, Float_t &value) {
  return GetBranchValue<Float_t>(branchName, value);
}

Bool_t TreeReader::GetBranchValue(const char *branchName, Double_t &value) {
  return GetBranchValue<Double_t>(branchName, value);
}

Bool_t TreeReader::GetBranchValue(const char *branchName, Bool_t &value) {
  return GetBranchValue<Bool_t>(branchName, value);
}

Bool_t TreeReader::GetBranchValue(const char *branchName, TString &value) {
  return GetBranchValue<TString>(branchName, value);
}

Bool_t TreeReader::GetBranchValue(const char *branchName, Char_t &value) {
  return GetBranchValue<Char_t>(branchName, value);
}
// Template method for reading branch values
template <typename T>
Bool_t TreeReader::GetBranchValue(const char *branchName, T &value) {
  if (!m_isOpen || !m_tree)
    return kFALSE;

  TBranch *branch = GetBranch(branchName);
  if (!branch) {
    std::cerr << "Error: Branch '" << branchName << "' not found." << std::endl;
    return kFALSE;
  }

  // Set branch address and get entry
  m_tree->SetBranchAddress(branchName, &value);
  branch->GetEntry(m_tree->GetReadEntry());

  return kTRUE;
}

// Process a single event
void TreeReader::ProcessEvent(Long64_t entry) {
  if (!m_isOpen)
    return;

  m_tree->GetEntry(entry);

  // Print the entry number every 1000 events
  // if (entry % 1000 == 0) {
  //     std::cout << "Processing entry " << entry << std::endl;
  // }
}

// Loop over all events
void TreeReader::Loop() {
  if (!m_isOpen) {
    std::cout << "No file open for looping." << std::endl;
    return;
  }

  std::cout << "Starting loop over " << m_entries << " entries..." << std::endl;

  for (Long64_t i = 0; i < m_entries; i++) {
    ProcessEvent(i);
  }

  std::cout << "Loop completed." << std::endl;
}

// Main analysis function
void TreeReader::Analyze() {
  if (!m_isOpen) {
    std::cout << "No file open for analysis." << std::endl;
    return;
  }

  PrintInfo();
  ListBranches();
}

// Read configuration file
void TreeReader::init_config(const char *config_file) {
  if (config.ReadFile(config_file, kEnvUser) < 0) {
    MSG(ERR, "[MAT] Config file not loaded --> exit.");
    exit(0);
  } else {
    RootFileName = config.GetValue("RootFileName", "");
    TreeObjectName = config.GetValue("TreeObjectName", "");
    HistTitle = config.GetValue("HistTitle", "");
    HistOutputName = config.GetValue("HistOutputName", "");
    SaveCanvas = config.GetValue("SaveCanvas", 0.0);
  }
}

// To be overridden by user
void TreeReader::SetBranchAddresses() {
  // This method is overriden in DoseCalculation.cxx to set branch addresses for
  // their specific variables
}
