void plotXY_Dose2D()
{
    // Open file
    TFile *f = new TFile("Electronics.root");
    TTree *t = (TTree*)f->Get("Electronics");

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

    int nx = 20, ny = 10;

    double xmin = -30, xmax = 230;
    double ymin = -140, ymax = 20;

    double binArea = ((xmax - xmin) / nx) * ((ymax - ymin) / ny); // mm^2
    double volume  = binArea * thickness;                          // mm^3
    double mass_g  = rho * volume;
    double mass_kg = mass_g * 1e-3;

    double MeV_to_J = 1.602e-13;

    // Gy per MeV
    double doseFactor_Gy = MeV_to_J / mass_kg;

    // Convert Gy → μrad
    double Gy_to_urad = 1e8; // 1 Gy = 100 rad, 1 rad = 1e6 μrad

    double doseFactor_urad = doseFactor_Gy * Gy_to_urad;

    // -----------------------------
    // Histogram (DOSE in μrad)
    // -----------------------------
    TH2D *h2 = new TH2D("h2",
        "Dose in Electronics;X [mm];Y [mm]",
        nx, xmin, xmax,
        ny, ymin, ymax);

    // Fill: Edep → μrad
    t->Draw("fY:fX>>h2", Form("Edep * %f", doseFactor_urad), "goff");

    h2->Smooth();
    h2->Draw("colz");

    h2->GetZaxis()->SetTitle("Dose [#murad]");

    gPad->SetGrid();
    c1->Update();
}

 plotXY_Dose2D();

