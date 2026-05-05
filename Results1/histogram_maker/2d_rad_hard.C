void plotXY_Edep2D() {

    // Open file
    TFile *f = new TFile("Electronics_Full.root");
    TTree *t = (TTree*)f->Get("Electronics");

    // Remove ROOT stats box
    gStyle->SetOptStat(0);

    // Canvas
    TCanvas *c1 = new TCanvas("c1", "Energy Deposition XY Map", 900, 700);
    c1->SetRightMargin(0.15);

    // Color palette
    gStyle->SetPalette(kBird);
    gStyle->SetNumberContours(100);

    // Create 2D histogram with units
    TH2D *h2 = new TH2D("h2",
        "Energy Deposition in Electronics;X [mm];Y [mm]",
        60, -30, 230,
        30, -140, 20);

    // Fill histogram (Edep as weight)
    t->Draw("fY:fX>>h2", "Edep", "goff");

    // Optional smoothing
    h2->Smooth();

    // Draw with color scale
    h2->Draw("colz");

    // Label Z axis (energy)
    h2->GetZaxis()->SetTitle("Deposited Energy [MeV]");

    // Cosmetics
    gPad->SetGrid();

    c1->Update();
}

 plotXY_Edep2D();