#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <iostream>
#include <cmath>
#include <map>

void plot_strip_edep()
{
    TFile *file = TFile::Open("/media/sf_vm-share/Neutron source output/rootoutput/testoutput.root");
    TTree *tree = (TTree*)file->Get("HeliumEdep");

    Int_t fEvent;
    Double_t Edep;
    Double_t fX;

    double W = 2.5;  // strip width
    double threshold = 1e-3;  // threshold to ignore zero depositions

    tree->SetBranchAddress("fEvent", &fEvent);
    tree->SetBranchAddress("Edep", &Edep);
    tree->SetBranchAddress("fX", &fX);

    Long64_t nentries = tree->GetEntries();

    double xmin = 85;
    double xmax = 115;
    int nbins = (xmax - xmin)/W;

    int currentEvent = -1;
    TH1D *hStrip = nullptr;
    TCanvas *c1 = new TCanvas("c1","Strip Edep",800,600);

    // Store mean X per event
    std::map<int, double> meanX_per_event;
    std::map<int,double> firstX_per_event;

    for(Long64_t i = 0; i < nentries; i++)
    {
        tree->GetEntry(i);
        if(Edep < threshold) continue;

        // New event
        if(fEvent != currentEvent)
        {
            // Draw previous event
            if(hStrip)
            {
                //hStrip->Draw("HIST");
                //c1->Update();

                // Calculate mean X
                double meanX = hStrip->GetMean();
                meanX_per_event[currentEvent] = meanX;
                //std::cout << "Event " << currentEvent << " mean X = " << meanX << std::endl;

                //std::cout << "Press Enter to continue...";
                //std::cin.ignore();

                delete hStrip;
            }

            currentEvent = fEvent;
            // store first X of the event
            firstX_per_event[currentEvent] = fX;
            hStrip = new TH1D(Form("hStrip_evt%d", currentEvent),
                               Form("Energy deposition per strip, Event %d;X position;Edep", currentEvent),
                               nbins, xmin, xmax);
        }

        int stripID = floor((fX-xmin) / W);
        double xpos = (stripID + 0.5)*W + xmin;
        hStrip->Fill(xpos, Edep);
    }

    // Draw last event
    if(hStrip)
    {
        hStrip->Draw("HIST");
        c1->Update();

        double meanX = hStrip->GetMean();
        meanX_per_event[currentEvent] = meanX;
        //std::cout << "Event " << currentEvent << " mean X = " << meanX << std::endl;

        //std::cout << "Press Enter to finish...";
        //std::cin.ignore();
        delete hStrip;
    }

    delete c1;

    // --- Create histogram of mean X values ---
    int nEvents = meanX_per_event.size();
    double hMin = xmin;   // min X
    double hMax = xmax;   // max X
    int nBinsMean =100;  // same bin width as strips

    TH1D *hMeanX = new TH1D("hMeanX","Distribution of Mean X per Event;Mean X;Counts", nBinsMean, hMin, hMax);

    for(auto const& entry : meanX_per_event)
        hMeanX->Fill(entry.second);

    // Draw mean X histogram
    TCanvas *c2 = new TCanvas("c2","Mean X Distribution",800,600);
    hMeanX->Draw("HIST");
    c2->Update();
    c2->SaveAs("meanX_distribution.png");

    std::cout << "\nHistogram of mean X values created and saved as meanX_distribution.png" << std::endl;

    TH1D *hDeltaX = new TH1D("hDeltaX",
                         "Mean X - First Hit X;#DeltaX;Counts",
                         3000, -15, 15);

    for(auto const& entry : meanX_per_event)
    {
        int eventID = entry.first;
        double meanX = entry.second;
        double firstX = firstX_per_event[eventID];

        double deltaX = meanX - firstX;

        hDeltaX->Fill(deltaX);
    }

        // Draw delta X histogram
        TCanvas *c3 = new TCanvas("c3","Delta X Distribution",800,600);
        hDeltaX->Draw("HIST");
        c3->Update();
        c3->SaveAs("deltaX_distribution.png");

        std::cout << "\nHistogram of Delta X values created and saved as deltaX_distribution.png" << std::endl;
    }