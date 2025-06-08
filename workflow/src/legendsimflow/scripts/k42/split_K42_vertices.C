struct mini_shroud {
    double height;
    double radius;
    double center_x;
    double center_y;
    double center_z;
};

// volume of lar-cylinder-K42 simulated region (estimated by MaGe): 3.72495e+06 cm^3
// TODO: update values
double in_fraction = 0.00195451;
double out_fraction = 0.99804549;

std::vector<mini_shroud> mss = {
    {956  , 51,  109   ,  188.79, -25},
    {960  , 48,  198.22,   95.44, -27},
    {652  , 44,  220   ,   0    , 127},
    {898  , 48,  198.22,  -95.44,   4},
    {543  , 54,  109   , -188.79, 181.5},
    {795  , 51, -109   , -188.79,  55.5},
    {959.5, 48, -198.22,  -95.44, -26.75},
    {966.5, 44, -220   ,    0   , -30.25},
    {951.5, 48, -198.22,   95.44, -22.75},
    {913  , 51, -109   ,  188.79,  -3.5},
};


bool point_is_inside_MS(double x, double y, double z, mini_shroud ms) {
    // translate coordinate system
    double tr_x = x - ms.center_x;
    double tr_y = y - ms.center_y;
    double tr_z = z - ms.center_z;

    // height inside MS AND radius inside MS
    return (abs(tr_z) <= ms.height/2) and (sqrt(tr_x * tr_x + tr_y * tr_y) <= ms.radius);
}


bool point_is_inside_any_MS(double x, double y, double z) {
    for (auto& ms: mss) {
        if (point_is_inside_MS(x, y, z, ms)) {
            return true;
        }
    }
    return false;
}


void usage() {
    std::cerr << "\n"
              << "USAGE: split_K42_vertices [options] input_file output_file_in output_file_out\n"
              << "\n"
              << "options:\n"
              << "  --help|-h : print this help message and exit\n"
              << "\n";
    gSystem->Exit(1);
}


void split_K42_vertices(std::string infile, std::string outfile_in, std::string outfile_out) {
    // open file
    auto file = TFile::Open(infile.c_str(), "READ");
    auto fTree = dynamic_cast<TTree*>(file->Get("fTree"));
    auto nev_obj = dynamic_cast<TNamed*>(file->Get("NumberOfEvents"));
    if (!nev_obj) {
        std::cerr << "ERROR: NumberOfEvents not found!" << std::endl;
        gSystem->Exit(1);
    }
    auto n_sim_ev = std::stoul(nev_obj->GetTitle());

    // get vertex position
    MGTMCEventSteps* evt = nullptr;
    fTree->SetBranchAddress("eventPrimaries", &evt);

    // create output files and copy original tree twice
    auto file_of_out = TFile::Open(outfile_out.c_str(), "RECREATE");
    auto fTree_out = dynamic_cast<TTree*>(fTree->CloneTree(0));
    unsigned long n_sim_ev_out = n_sim_ev * out_fraction;
    auto nev_out = new TNamed("NumberOfEvents", std::to_string(n_sim_ev_out));
    auto file_of_in = TFile::Open(outfile_in.c_str(), "RECREATE");
    auto fTree_in  = dynamic_cast<TTree*>(fTree->CloneTree(0));
    unsigned long n_sim_ev_in = n_sim_ev * in_fraction;
    auto nev_in = new TNamed("NumberOfEvents", std::to_string(n_sim_ev_in));

    int nevents = fTree->GetEntries();
    int n_inside = 0;

    // loop events and write them in the right tree
    for(int e = 0; e < nevents; e++) {
        fTree->GetEntry(e);

        if (!evt) std::cerr << "ERROR: eventPrimaries == nullptr" << std::endl;

        if (evt->GetNSteps() != 1) {
            std::cerr << "ERROR: number of steps in eventPrimaries != 1, skipping event" << std::endl;
            continue;
        }

        auto v = evt->GetStep(0)->GetPositionVector();

        if (point_is_inside_any_MS(v.X(), v.Y(), v.Z())) {
            fTree_in->Fill();
            n_inside++;
        }
        else {
            fTree_out->Fill();
        }
    }

    std::cout << "INFO: total entries: " << nevents << " of which "
              << n_inside << " inside the mini-shrouds" << std::endl;

    file->Close();

    file_of_out->cd();
    fTree_out->Write();
    nev_out->Write();
    file_of_out->Close();

    file_of_in->cd();
    fTree_in->Write();
    nev_in->Write();
    file_of_in->Close();
}
