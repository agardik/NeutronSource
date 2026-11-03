void plotXY_Edep2D() {

    // Open file
    TFile *f = TFile::Open("testoutput2.root");
    if (!f || f->IsZombie()) { printf("Cannot open file\n"); return; }
    TTree *t = (TTree*)f->Get("Hits");

    // Remove ROOT stats box
    gStyle->SetOptStat(0);

    // Canvas
    TCanvas *c1 = new TCanvas("c1", "Energy Deposition XY Map", 900, 700);
    c1->SetRightMargin(0.15);

    // Color palette
    gStyle->SetPalette(kBird);
    gStyle->SetNumberContours(100);

    // 2D histogram with units
    TH2D *h2 = new TH2D("h2",
        "Energy Deposition in Electronics;X [mm];Y [mm]",
        9, 0, 200,
        6, -120, 0);

    // The new tree has no fX/fY branches, only pre/post step positions.
    // Use the step midpoint as the hit position, weighted by Edep.
    // Optional: restrict to a volume, e.g. "Edep*(fVolumeName==\"Electronics\")"
    t->Draw("(fYpre+fYpost)/2.:(fXpre+fXpost)/2.>>h2", "Edep", "goff");

    // Draw with color scale
    h2->Draw("colz");
    h2->GetZaxis()->SetTitle("Deposited Energy [MeV]");

    // Cosmetics
    gPad->SetGrid();
    c1->Update();
}

 plotXY_Edep2D();