void plotXY_Edep2D() {

    // Open file
    TFile *f = new TFile("Electronics_Top_Strip_8_test_kill.root");
    TTree *t = (TTree*)f->Get("Hits");

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
        50, 0, 200,
        30, -120, 0);

    // Fill histogram (Edep as weight)
    t->Draw("fY:fX>>h2", "Edep*(fParticleName==\"gamma\")", "goff");

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