void plotXY_Dose2D()
{
    // Open file
    TFile *f = TFile::Open("testoutput2.root");   // new file
    if (!f || f->IsZombie()) { printf("Cannot open file\n"); return; }
    TTree *t = (TTree*)f->Get("Hits");

    gStyle->SetOptStat(0);

    TCanvas *c1 = new TCanvas("c1", "Dose XY Map", 900, 700);
    c1->SetRightMargin(0.15);

    gStyle->SetPalette(kBird);
    gStyle->SetNumberContours(100);

    // -----------------------------
    // Geometry + material
    // -----------------------------
    double rho = 2.33e-3;   // g/mm^3 (Si)
    double thickness = 2.0; // mm

    int nx = 18, ny = 12;

    double xmin = 0, xmax = 200;
    double ymin = -120, ymax = 0;

    double binArea = ((xmax - xmin) / nx) * ((ymax - ymin) / ny); // mm^2
    double volume  = binArea * thickness;                          // mm^3
    double mass_g  = rho * volume;
    double mass_kg = mass_g * 1e-3;

    double MeV_to_J = 1.602e-13;

    // Gy per MeV
    double doseFactor_Gy = MeV_to_J / mass_kg;

    // Convert Gy -> urad (1 Gy = 100 rad, 1 rad = 1e6 urad)
    double Gy_to_urad = 1e8;

    double doseFactor_urad = doseFactor_Gy * Gy_to_urad;

    // -----------------------------
    // Histogram (DOSE in urad)
    // -----------------------------
    TH2D *h2 = new TH2D("h2",
        "Dose in Electronics;X [mm];Y [mm]",
        nx, xmin, xmax,
        ny, ymin, ymax);

    // New tree has no fX/fY: use the step midpoint (pre/post average).
    // %g keeps full precision of the small conversion factor.
    // Optional volume cut: Form("Edep * %g * (fVolumeName==\"YourVolume\")", ...)
    t->Draw("(fYpre+fYpost)/2.:(fXpre+fXpost)/2.>>h2",
            Form("Edep * %g", doseFactor_urad), "goff");

    h2->Draw("colz");
    h2->GetZaxis()->SetTitle("Dose [#murad]");

    gPad->SetGrid();
    c1->Update();
}