#ifndef TREEREADER_H
#define TREEREADER_H

#include "TEnv.h"
#include "cout_msg.h"
#include <TBranch.h>
#include <TFile.h>
#include <TString.h>
#include <TStyle.h>
#include <TTree.h>
#include <iostream>
#include <map>
#include <vector>

class TreeReader {
private:
  TFile *m_file;      // ROOT file pointer
  TTree *m_tree;      // Tree pointer
  TString m_filename; // Input file name
  TString m_treename; // Tree name
  Long64_t m_entries; // Number of entries in the tree
  Bool_t m_isOpen;    // File open status

  // Branch information
  std::vector<TString> m_branchNames;
  std::map<TString, TBranch *> m_branches;
  std::map<TString, TString> m_branchTypes;

public:
  // Constructor
  TreeReader();
  TreeReader(const char *filename, const char *treename = "tree");

  // Destructor
  virtual ~TreeReader();

  // Main methods
  Bool_t OpenFile(const char *filename, const char *treename = "tree");
  void CloseFile();

  // Getters
  TFile *GetFile() const { return m_file; }
  TTree *GetTree() const { return m_tree; }
  Long64_t GetEntries() const { return m_entries; }
  Bool_t IsOpen() const { return m_isOpen; }
  const char *GetFileName() const { return m_filename.Data(); }
  const char *GetTreeName() const { return m_treename.Data(); }

  // Branch management
  Bool_t InitializeBranches();
  std::vector<TString> GetBranchNames() const { return m_branchNames; }
  TBranch *GetBranch(const char *branchName);
  TString GetBranchType(const char *branchName);
  Int_t GetNBranches() const { return m_branchNames.size(); }

  // Utility methods
  void PrintInfo() const;
  void ListBranches() const;
  void PrintBranchInfo(const char *branchName) const;

  // Branch value reading methods
  Bool_t GetBranchValue(const char *branchName, Int_t &value);
  Bool_t GetBranchValue(const char *branchName, Float_t &value);
  Bool_t GetBranchValue(const char *branchName, Double_t &value);
  Bool_t GetBranchValue(const char *branchName, Bool_t &value);
  Bool_t GetBranchValue(const char *branchName, TString &value);
  Bool_t GetBranchValue(const char *branchName, Char_t &value);

  // Template method for any type
  template <typename T> Bool_t GetBranchValue(const char *branchName, T &value);

  // Analysis methods
  virtual void ProcessEvent(Long64_t entry);
  virtual void Loop();
  virtual void Analyze();

  // Parameters to be read from config file
  std::string RootFileName{""};    // Path to the input root file.
  const char *TreeObjectName;      // Name of the tree object in the root
  std::string HistTitle = "";      // Name of the histogram title
  std::string HistOutputName = ""; // Name of the histogram output name
  Int_t SaveCanvas = 0;            // Whether to save canvas or not

  // Configuration object
  TEnv config;

  // Method for reading configuration file
  void init_config(const char *config_file);


protected:
  virtual void SetBranchAddresses(); // To be overridden by user
};

#endif
