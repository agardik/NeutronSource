void plotXY_Edep2D() {

    // Open file
    TFile *f = new TFile("/media/sf_vm-share/Neutron source output/rootoutput/Electronics8.root");
    TTree *t = (TTree*)f->Get("Electronics");

    // Canvas
    TCanvas *c1 = new TCanvas("c1", "Edep XY Map", 900, 700);
    c1->SetRightMargin(0.15);

    // Better color palette
    gStyle->SetPalette(kBird);   // nice blue→red gradient
    gStyle->SetNumberContours(100);

    // Create 2D histogram
    TH2D *h2 = new TH2D("h2", "Energy Deposition Map;X;Y",
                        500, -100, 300,
                        500, -150, 100);

    // Fill with weight = Edep
    t->Draw("fY:fX>>h2", "Edep",  "colz");

    // Draw with color
    h2->Draw("colz");

    // Optional: smooth a bit
    h2->Smooth();

    // Grid + cosmetics
    gPad->SetGrid();
    c1->Update();
}

 plotXY_Edep2D()