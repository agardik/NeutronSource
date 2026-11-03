void plotSphere3D_pretty() {

    // -------------------------
    // Open file & get tree
    // -------------------------
    TFile *f = TFile::Open("combine_sphere.root");
    if (!f || f->IsZombie()) {
        Error("plotSphere3D_pretty", "Cannot open file!");
        return;
    }

    TTree *t = (TTree*)f->Get("Sphere");
    if (!t) {
        Error("plotSphere3D_pretty", "Tree 'Sphere' not found!");
        return;
    }

    // -------------------------
    // Global style tweaks
    // -------------------------
    gStyle->SetOptStat(0);
    gStyle->SetTitleFontSize(0.04);
    gStyle->SetPalette(kBird); // fallback palette

    // -------------------------
    // Canvas
    // -------------------------
    TCanvas *c1 = new TCanvas("c1", "3D Energy Deposition", 1000, 800);
    c1->SetRightMargin(0.20);
    c1->SetLeftMargin(0.10);
    c1->SetBottomMargin(0.10);

    // -------------------------
    // Custom gradient palette
    // -------------------------
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;

    Double_t stops[NRGBs] = {0.00, 0.34, 0.61, 0.84, 1.00};
    Double_t red[NRGBs]   = {0.00, 0.00, 0.87, 1.00, 0.51};
    Double_t green[NRGBs] = {0.00, 0.81, 1.00, 0.20, 0.00};
    Double_t blue[NRGBs]  = {0.51, 1.00, 0.12, 0.00, 0.00};

    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);

    // -------------------------
    // Frame histogram (axes)
    // -------------------------
    TH3F *frame = new TH3F("frame",
        "3D Energy Deposition;X [mm];Y [mm];Z [mm]",
        10, -400, 400,
        10, -400, 400,
        10, -400, 400
    );

    frame->GetXaxis()->CenterTitle();
    frame->GetYaxis()->CenterTitle();
    frame->GetZaxis()->CenterTitle();
    frame->Draw();

    // -------------------------
    // Draw 3D scatter with color (Edep)
    // -------------------------
    t->Draw("fZ:fY:fX:Edep", "", "glcol same");

    // -------------------------
    // Visual improvements
    // -------------------------
    gPad->SetGrid(1,1);
    gPad->SetTheta(25);
    gPad->SetPhi(40);

    // Slightly smoother rendering (if OpenGL is enabled)
    gStyle->SetCanvasPreferGL(kTRUE);

    // -------------------------
    // Final update
    // -------------------------
    c1->Modified();
    c1->Update();
}

 plotSphere3D_pretty();