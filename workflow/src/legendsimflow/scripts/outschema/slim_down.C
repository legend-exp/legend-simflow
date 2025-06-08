void usage() {
    std::cerr << "\n"
              << "USAGE: slim_down.C [options] input_file output_file\n"
              << "\n"
              << "options:\n"
              << "  --help|-h : print this help message and exit\n"
              << "\n";
    gSystem->Exit(1);
}


void slim_down(std::string infile, std::string outfile) {
    // open file
    auto file = TFile::Open(infile.c_str(), "READ");
    auto fTree = dynamic_cast<TTree*>(file->Get("fTree"));
    auto nev_obj = dynamic_cast<TNamed*>(file->Get("NumberOfEvents"));
    if (!nev_obj) {
        std::cerr << "ERROR: NumberOfEvents not found!" << std::endl;
        gSystem->Exit(1);
    }

    // deactivate some branches
    std::vector<std::string> bad_branches = {
        "eventPrimaries.fSteps.fPhysVolName",
        "eventPrimaries.fSteps.fProcessName",
        "eventPrimaries.fSteps.fTotalTrackLength",
        "eventPrimaries.fSteps.fPx",
        "eventPrimaries.fSteps.fPy",
        "eventPrimaries.fSteps.fPz",
        "eventSteps.fSteps.fPhysVolName",
        "eventSteps.fSteps.fProcessName",
        "eventSteps.fSteps.fTotalTrackLength",
        "eventSteps.fSteps.fPx",
        "eventSteps.fSteps.fPy",
        "eventSteps.fSteps.fPz",
    };

    for (auto br : bad_branches) {
        fTree->SetBranchStatus(br.c_str(), 0);
    }

    // create output files and copy original tree without branches above
    auto file_out = TFile::Open(outfile.c_str(), "RECREATE");
    auto fTree_out = fTree->CloneTree();

    // write out
    file->Close();

    file_out->cd();
    fTree_out->Write();
    file_out->WriteObject(nev_obj, "NumberOfEvents");
    file_out->Close();
}
