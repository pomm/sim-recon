TFile * thisFile;

TH1I * Get1DHistogram(const char * plugin, const char * directoryName, const char * name){
    TH1I * histogram;
    TString fullName = TString(plugin) + "/" + TString(directoryName) + "/" + TString(name);
    thisFile->GetObject(fullName, histogram);
    if (histogram == 0){
        cout << "Unable to find histogram " << fullName.Data() << endl;
        return NULL;
    }
    return histogram;
}

TH2I * Get2DHistogram(const char * plugin, const char * directoryName, const char * name){
    TH2I * histogram;
    TString fullName = TString(plugin) + "/" + TString(directoryName) + "/" + TString(name);
    histogram = thisFile->GetObject(fullName, histogram);
    if (histogram == 0){
        cout << "Unable to find histogram " << fullName.Data() << endl;
        return NULL;
    }
    return histogram;
}

void GetCCDBConstants(TString path, Int_t run, TString variation, vector<double>& container, Int_t column = 1){
    char command[256];

    sprintf(command, "ccdb dump %s:%i:%s", path.Data(), run, variation.Data());
    FILE* inputPipe = gSystem->OpenPipe(command, "r");
    if(inputPipe == NULL)
        return 0;
    //get the first (comment) line
    char buff[1024];
    if(fgets(buff, sizeof(buff), inputPipe) == NULL)
        return 0;
    //get the remaining lines
    double entry;
    int counter = 0;
    while(fgets(buff, sizeof(buff), inputPipe) != NULL){
        istringstream locConstantsStream(buff);
        while (locConstantsStream >> entry){
            counter++;
            if (counter % column == 0) container.push_back(entry);
        }
    }
    //Close the pipe
    gSystem->ClosePipe(inputPipe);
}

//Overload this function to handle the base time offsets
void GetCCDBConstants1(TString path, Int_t run, TString variation, double& constant1){
    char command[256];
    sprintf(command, "ccdb dump %s:%i:%s", path.Data(), run, variation.Data());
    FILE* inputPipe = gSystem->OpenPipe(command, "r");
    if(inputPipe == NULL)
        return 0;
    //get the first (comment) line
    char buff[1024];
    if(fgets(buff, sizeof(buff), inputPipe) == NULL)
        return 0;
    //get the line containing the values
    while(fgets(buff, sizeof(buff), inputPipe) != NULL){
        istringstream locConstantsStream(buff);
        locConstantsStream >> constant1;
    }
    //Close the pipe
    gSystem->ClosePipe(inputPipe);
}

void GetCCDBConstants2(TString path, Int_t run, TString variation, double& constant1, double& constant2){
    char command[256];
    sprintf(command, "ccdb dump %s:%i:%s", path.Data(), run, variation.Data());
    FILE* inputPipe = gSystem->OpenPipe(command, "r");
    if(inputPipe == NULL)
        return 0;
    //get the first (comment) line
    char buff[1024];
    if(fgets(buff, sizeof(buff), inputPipe) == NULL)
        return 0;
    //get the line containing the values
    while(fgets(buff, sizeof(buff), inputPipe) != NULL){
        istringstream locConstantsStream(buff);
        locConstantsStream >> constant1 >> constant2;
    }
    //Close the pipe
    gSystem->ClosePipe(inputPipe);
}

int GetCCDBIndexTAGM(unsigned int column, unsigned int row){
    int CCDBIndex = column + row;
    if (column > 9) CCDBIndex += 5;
    if (column > 27) CCDBIndex += 5;
    if (column > 81) CCDBIndex += 5;
    if (column > 99) CCDBIndex += 5;

    return CCDBIndex;
}

void ExtractTrackBasedTiming(TString fileName = "hd_root.root", int runNumber = 2931, TString variation = "default", bool verbose = false,TString prefix = ""){

    // set "prefix" in case you want to ship the txt files elsewhere...
    cout << "Performing TDC/ADC timing fits for File: " << fileName.Data() << " Run: " << runNumber << " Variation: " << variation.Data() << endl;

    thisFile = TFile::Open( fileName , "UPDATE");
    if (thisFile == 0) {
        cout << "Unable to open file " << fileName.Data() << "...Exiting" << endl;
        return;
    }

    //We need the existing constants, The best we can do here is just read them from the file.
    vector<double> sc_tdc_time_offsets;
    vector<double> sc_fadc_time_offsets;
    vector<double> tof_tdc_time_offsets;
    vector<double> tof_fadc_time_offsets;
    vector<double> tagm_tdc_time_offsets;
    vector<double> tagm_fadc_time_offsets;
    vector<double> tagh_tdc_time_offsets;
    vector<double> tagh_fadc_time_offsets;

    double sc_t_base_fadc, sc_t_base_tdc;
    double tof_t_base_fadc, tof_t_base_tdc;
    double bcal_t_base_fadc, bcal_t_base_tdc;
    double tagm_t_base_fadc, tagm_t_base_tdc;
    double tagh_t_base_fadc, tagh_t_base_tdc;
    double fcal_t_base;
    double cdc_t_base;

    cout << "Grabbing CCDB constants..." << endl;
    // Base times
    GetCCDBConstants1("/CDC/base_time_offset" ,runNumber, variation, cdc_t_base);
    GetCCDBConstants1("/FCAL/base_time_offset",runNumber, variation, fcal_t_base);
    //GetCCDBConstants2("/FDC/base_time_offset" ,runNumber, variation, fdc_base_time_adc, fdc_base_time_tdc);
    GetCCDBConstants2("/BCAL/base_time_offset" ,runNumber, variation, bcal_t_base_fadc, bcal_t_base_tdc);
    GetCCDBConstants2("/PHOTON_BEAM/microscope/base_time_offset" ,runNumber, variation, tagm_t_base_fadc, tagm_t_base_tdc);
    GetCCDBConstants2("/PHOTON_BEAM/hodoscope/base_time_offset" ,runNumber, variation, tagh_t_base_fadc, tagh_t_base_tdc);
    GetCCDBConstants2("/START_COUNTER/base_time_offset" ,runNumber, variation, sc_t_base_fadc, sc_t_base_tdc);
    GetCCDBConstants2("/TOF/base_time_offset" ,runNumber, variation, tof_t_base_fadc, tof_t_base_tdc);
    // Per channel
    //GetCCDBConstants("/BCAL/TDC_offsets"    ,runNumber, variation, bcal_tdc_offsets);
    //GetCCDBConstants("/FCAL/timing_offsets" ,runNumber, variation, fcal_adc_offsets);
    GetCCDBConstants("/START_COUNTER/adc_timing_offsets" ,runNumber, variation, sc_fadc_time_offsets);
    GetCCDBConstants("/START_COUNTER/tdc_timing_offsets" ,runNumber, variation, sc_tdc_time_offsets);
    GetCCDBConstants("/PHOTON_BEAM/microscope/fadc_time_offsets" ,runNumber, variation, tagm_fadc_time_offsets,3);// Interested in 3rd column
    GetCCDBConstants("/PHOTON_BEAM/microscope/tdc_time_offsets"  ,runNumber, variation, tagm_tdc_time_offsets,3);
    GetCCDBConstants("/PHOTON_BEAM/hodoscope/fadc_time_offsets"  ,runNumber, variation, tagh_fadc_time_offsets,2);// Interested in 2nd column
    GetCCDBConstants("/PHOTON_BEAM/hodoscope/tdc_time_offsets"   ,runNumber, variation, tagh_tdc_time_offsets,2);
    GetCCDBConstants("/TOF/adc_timing_offsets",runNumber, variation, tof_fadc_time_offsets);
    GetCCDBConstants("/TOF/timing_offsets",runNumber, variation, tof_tdc_time_offsets);

    cout << "CDC base times = " << cdc_t_base << endl;
    cout << "FCAL base times = " << fcal_t_base << endl;
    //cout << "FDC base times = " << fdc_t_base_fadc << ", " << fdc_t_base_tdc << endl;
    cout << "BCAL base times = " << bcal_t_base_fadc << ", " << bcal_t_base_tdc << endl;
    cout << "SC base times = " << sc_t_base_fadc << ", " << sc_t_base_tdc << endl;
    cout << "TOF base times = " << tof_t_base_fadc << ", " << tof_t_base_tdc << endl;
    cout << "TAGH base times = " << tagh_t_base_fadc << ", " << tagh_t_base_tdc << endl;
    cout << "TAGH base times = " << tagm_t_base_fadc << ", " << tagm_t_base_tdc << endl;
    cout << endl;

    cout << "Done grabbing CCDB constants...Entering fits..." << endl;

    // Do our final step in the timing alignment with tracking

    //When the RF is present we can try to simply pick out the correct beam bucket for each of the runs
    //First just a simple check to see if we have the appropriate data
    bool useRF = false;
    double RF_Period = 4.0080161;
    TH1I *testHist = Get1DHistogram("HLDetectorTiming", "TAGH_TDC_RF_Compare","Counter ID 001");
    if (testHist != NULL){ // Not great since we rely on channel 1 working, but can be craftier later.
        cout << "Using RF Times for Calibration" << endl;
        useRF = true;
    }
    ofstream outFile;
    TH2I *thisHist; 
    thisHist = Get2DHistogram("HLDetectorTiming", "TRACKING", "TAGM - SC Target Time");
    if (useRF) thisHist = Get2DHistogram("HLDetectorTiming", "TRACKING", "TAGM - RFBunch Time");
    if (thisHist != NULL){
        //Statistics on these histograms are really quite low we will have to rebin and do some interpolation
        outFile.open(prefix + "tagm_tdc_timing_offsets.txt", ios::out | ios::trunc);
        outFile.close(); // clear file
        outFile.open(prefix + "tagm_adc_timing_offsets.txt", ios::out | ios::trunc);
        outFile.close(); // clear file
        int nBinsX = thisHist->GetNbinsX();
        int nBinsY = thisHist->GetNbinsY();
        TH1D * selectedTAGMOffset = new TH1D("selectedTAGMOffset", "Selected TAGM Offset; Column; Offset [ns]", nBinsX, 0.5, nBinsX + 0.5);
        TH1I * TAGMOffsetDistribution = new TH1I("TAGMOffsetDistribution", "TAGM Offset; TAGM Offset [ns]; Entries", 500, -250, 250);
        for (int i = 1 ; i <= nBinsX; i++){ 
            TH1D *projY = thisHist->ProjectionY("temp", i, i);
            // Scan over the histogram
            //chose the correct number of bins based on the histogram
            float nsPerBin = (projY->GetBinCenter(projY->GetNbinsX()) - projY->GetBinCenter(1)) / projY->GetNbinsX();
            float timeWindow = 3; //ns (Full Width)
            int binWindow = int(timeWindow / nsPerBin);
            double maxEntries = 0;
            double maxMean = 0;
            for (int j = 1 ; j <= projY->GetNbinsX();j++){
                int minBin = j;
                int maxBin = (j + binWindow) <= projY->GetNbinsX() ? (j + binWindow) : projY->GetNbinsX();
                double sum = 0, nEntries = 0;
                for (int bin = minBin; bin <= maxBin; bin++){
                    sum += projY->GetBinContent(bin) * projY->GetBinCenter(bin);
                    nEntries += projY->GetBinContent(bin);
                    if (bin == maxBin){
                        if (nEntries > maxEntries) {
                            maxMean = sum / nEntries;
                            maxEntries = nEntries;
                        }
                    } 
                }
            }
            //In the case there is RF, our job is to pick just the number of the correct beam bunch, so that's really all we need.
            if(useRF) {
                int beamBucket = int((maxMean / RF_Period) + 0.5); // +0.5 to handle rounding correctly
                selectedTAGMOffset->SetBinContent(i, beamBucket);
                TAGMOffsetDistribution->Fill(beamBucket);
            }
            else{
                selectedTAGMOffset->SetBinContent(i, maxMean);
                TAGMOffsetDistribution->Fill(maxMean);
            }
        }
        double meanOffset = TAGMOffsetDistribution->GetMean();
        // This might be in units of beam bunches, so we need to convert
        if (useRF) meanOffset *= RF_Period;
        if (verbose) {
            cout << "Dumping TAGM results...\n=======================================" << endl;
            cout << "TAGM mean Offset = " << meanOffset << endl;
            cout << "fADC Offsets" << endl;
        }

        outFile.open(prefix + "tagm_adc_timing_offsets.txt", ios::out);
        //for (int i = 1 ; i <= nBinsX; i++){
        // Loop over rows
        if (verbose) cout << "Column\tRow\tvalueToUse\toldValue\tmeanOffset\tTotal" << endl;
        for (unsigned int column = 1; column <= 102; column++){
            int index = GetCCDBIndexTAGM(column, 0);
            double valueToUse = selectedTAGMOffset->GetBinContent(index);
            if (useRF) valueToUse *= RF_Period;

            //if (valueToUse == 0) valueToUse = meanOffset;
            outFile << "0 " << column << " " << valueToUse + tagm_fadc_time_offsets[index-1] - meanOffset<< endl;
            if (verbose) printf("0\t%i\t%.3f\t\t%.3f\t\t%.3f\t\t%.3f\n", column, valueToUse, tagm_fadc_time_offsets[index-1], meanOffset, 
                    valueToUse + tagm_fadc_time_offsets[index-1] - meanOffset);
            if (column == 9 || column == 27 || column == 81 || column == 99){
                for (unsigned int row = 1; row <= 5; row++){
                    index = GetCCDBIndexTAGM(column, row);
                    valueToUse = selectedTAGMOffset->GetBinContent(index);
                    if (useRF) valueToUse *= RF_Period;
                    //if (valueToUse == 0) valueToUse = meanOffset;
                    outFile << row << " " << column << " " << valueToUse + tagm_fadc_time_offsets[index-1] - meanOffset<< endl;
                    if (verbose) printf("%i\t%i\t%.3f\t\t%.3f\t\t%.3f\t\t%.3f\n", row, column, valueToUse, tagm_fadc_time_offsets[index-1], meanOffset,
                            valueToUse + tagm_fadc_time_offsets[index-1] - meanOffset);
                }
            }
        }
        outFile.close();

        if (verbose) {
            cout << "TDC Offsets" << endl;
            cout << "Column\tRow\tvalueToUse\toldValue\tmeanOffset\tTotal" << endl;
        }
        outFile.open(prefix + "tagm_tdc_timing_offsets.txt", ios::out);
        //for (int i = 1 ; i <= nBinsX; i++){
        // Loop over rows
        for (unsigned int column = 1; column <= 102; column++){
            int index = GetCCDBIndexTAGM(column, 0);
            double valueToUse = selectedTAGMOffset->GetBinContent(index);
            if (useRF) valueToUse *= RF_Period;
            //if (valueToUse == 0) valueToUse = meanOffset;
            outFile << "0 " << column << " " << valueToUse + tagm_tdc_time_offsets[index-1] - meanOffset << endl;
            if (verbose) printf("0\t%i\t%.3f\t\t%.3f\t\t%.3f\t\t%.3f\n", column, valueToUse, tagm_tdc_time_offsets[index-1], meanOffset,
                    valueToUse + tagm_tdc_time_offsets[index-1] - meanOffset);
            if (column == 9 || column == 27 || column == 81 || column == 99){
                for (unsigned int row = 1; row <= 5; row++){
                    index = GetCCDBIndexTAGM(column, row);
                    valueToUse = selectedTAGMOffset->GetBinContent(index);
                    if (useRF) valueToUse *= RF_Period;
                    //if (valueToUse == 0) valueToUse = meanOffset;
                    outFile << row << " " << column << " " << valueToUse + tagm_tdc_time_offsets[index-1] - meanOffset << endl;
                    if (verbose) printf("%i\t%i\t%.3f\t\t%.3f\t\t%.3f\t\t%.3f\n", row, column, valueToUse, tagm_tdc_time_offsets[index-1], meanOffset,
                            valueToUse + tagm_tdc_time_offsets[index-1] - meanOffset);
                }
            }
        }
        outFile.close();
        outFile.open(prefix + "tagm_base_time.txt", ios::out);
        if (verbose) {
            printf("TAGM ADC Base = %f - (%f) = %f\n", tagm_t_base_fadc, meanOffset, tagm_t_base_fadc - meanOffset);
            printf("TAGM TDC Base = %f - (%f) = %f\n", tagm_t_base_tdc, meanOffset, tagm_t_base_tdc - meanOffset);
        }
        outFile << tagm_t_base_fadc - meanOffset << " " << tagm_t_base_tdc - meanOffset << endl;
        outFile.close();

    }

    thisHist = Get2DHistogram("HLDetectorTiming", "TRACKING", "TAGH - SC Target Time");
    if (useRF) thisHist = Get2DHistogram("HLDetectorTiming", "TRACKING", "TAGH - RFBunch Time");
    if(thisHist != NULL){
        outFile.open(prefix + "tagh_tdc_timing_offsets.txt", ios::out | ios::trunc);
        outFile.close(); // clear file
        outFile.open(prefix + "tagh_adc_timing_offsets.txt", ios::out | ios::trunc);
        outFile.close(); // clear file

        int nBinsX = thisHist->GetNbinsX();
        int nBinsY = thisHist->GetNbinsY();
        TH1D * selectedTAGHOffset = new TH1D("selectedTAGHOffset", "Selected TAGH Offset; ID; Offset [ns]", nBinsX, 0.5, nBinsX + 0.5);
        TH1I * TAGHOffsetDistribution = new TH1I("TAGHOffsetDistribution", "TAGH Offset; TAGH Offset [ns]; Entries", 500, -250, 250);
        for (int i = 1 ; i <= nBinsX; i++){
            TH1D *projY = thisHist->ProjectionY("temp", i, i);
            // Scan over the histogram
            //chose the correct number of bins based on the histogram
            float nsPerBin = (projY->GetBinCenter(projY->GetNbinsX()) - projY->GetBinCenter(1)) / projY->GetNbinsX();
            float timeWindow = 2; //ns (Full Width)
            int binWindow = int(timeWindow / nsPerBin);

            double maxEntries = 0;
            double maxMean = 0;
            for (int j = 1 ; j <= projY->GetNbinsX();j++){
                int minBin = j;
                int maxBin = (j + binWindow) <= projY->GetNbinsX() ? (j + binWindow) : projY->GetNbinsX();
                double sum = 0; 
                double nEntries = 0;
                for (int bin = minBin; bin <= maxBin; bin++){
                    sum += projY->GetBinContent(bin) * projY->GetBinCenter(bin);
                    nEntries += projY->GetBinContent(bin);
                    if (bin == maxBin){
                        if (nEntries > maxEntries){
                            maxMean = sum / nEntries;
                            maxEntries = nEntries;
                        }
                    }
                }
            }

            if(useRF) {
                int beamBucket;
                if (maxMean >= 0) beamBucket = int((maxMean / RF_Period) + 0.5); // +0.5 to handle rounding correctly
                else beamBucket = int((maxMean / RF_Period) - 0.5);
                selectedTAGHOffset->SetBinContent(i, beamBucket);
                TAGHOffsetDistribution->Fill(beamBucket);
            }
            else{
                selectedTAGHOffset->SetBinContent(i, maxMean);
            }
        }

        double meanOffset = TAGHOffsetDistribution->GetMean();
        if (useRF) meanOffset *= RF_Period;
        if (verbose) {
            cout << "Dumping TAGH results...\n=======================================" << endl;
            cout << "TAGH mean Offset = " << meanOffset << endl;
            cout << "Type\tChannel\tvalueToUse\toldValue\tmeanOffset\tTotal" << endl;
        }
        for (int i = 1 ; i <= nBinsX; i++){
            valueToUse = selectedTAGHOffset->GetBinContent(i);
            if (useRF) valueToUse *= RF_Period;
            //if (valueToUse == 0) valueToUse = meanOffset;
            outFile.open(prefix + "tagh_tdc_timing_offsets.txt", ios::out | ios::app);
            outFile << i << " " << valueToUse + tagh_tdc_time_offsets[i-1] - meanOffset << endl;
            if (verbose) printf("TDC\t%i\t%.3f\t\t%.3f\t\t%.3f\t\t%.3f\n", i, valueToUse, tagh_tdc_time_offsets[i-1], meanOffset,
                    valueToUse + tagh_tdc_time_offsets[i-1] - meanOffset);
            outFile.close();
            outFile.open(prefix + "tagh_adc_timing_offsets.txt", ios::out | ios::app);
            outFile << i << " " << valueToUse + tagh_fadc_time_offsets[i-1] - meanOffset << endl;
            if (verbose) printf("ADC\t%i\t%.3f\t\t%.3f\t\t%.3f\t\t%.3f\n", i, valueToUse, tagh_fadc_time_offsets[i-1], meanOffset,
                    valueToUse + tagh_fadc_time_offsets[i-1] - meanOffset);
            outFile.close();
        }

        outFile.open(prefix + "tagh_base_time.txt", ios::out);
        outFile << tagh_t_base_fadc - meanOffset << " " << tagh_t_base_tdc - meanOffset << endl;
        if (verbose) {
            printf("TAGH ADC Base = %f - (%f) = %f\n", tagh_t_base_fadc, meanOffset, tagh_t_base_fadc - meanOffset);
            printf("TAGH TDC Base = %f - (%f) = %f\n", tagh_t_base_tdc, meanOffset, tagh_t_base_tdc - meanOffset);
        }
        outFile.close();
    }

    // We can use the RF time to calibrate the SC time (Experimental for now)
    double meanSCOffset = 0.0; // In case we change the time of the SC, we need this in this scope
    if(useRF){
        TH1F * selectedSCSectorOffset = new TH1F("selectedSCSectorOffset", "Selected TDC-RF offset;Sector; Time", 30, 0.5, 30.5);
        TH1F * selectedSCSectorOffsetDistribution = new TH1F("selectedSCSectorOffsetDistribution", "Selected TDC-RF offset;Time;Entries", 100, -3.0, 3.0);
        TF1* f = new TF1("f","pol0(0)+gaus(1)", -3.0, 3.0);
        for (int sector = 1; sector <= 30; sector++){
            TH1I *scRFHist = Get1DHistogram("HLDetectorTiming", "SC_Target_RF_Compare", Form("Sector %.2i", sector));
            if (scRFHist == NULL) continue;
            //Do the fit
            TFitResultPtr fr = scRFHist->Fit("pol0", "SQ", "", -2, 2);
            double p0 = fr->Parameter(0);

            f->FixParameter(0,p0);
            f->SetParLimits(2, -2, 2);
            f->SetParLimits(3, 0, 2);
            f->SetParameter(1, 10);
            f->SetParameter(2, scRFHist->GetBinCenter(scRFHist->GetMaximumBin()));
            f->SetParameter(3, 0);

            fr = scRFHist->Fit(f, "SQ", "", -2, 2);
            double SCOffset = fr->Parameter(2);
            selectedSCSectorOffset->SetBinContent(sector, SCOffset);
            selectedSCSectorOffsetDistribution->Fill(SCOffset);
        }
        // Now write out the offsets
        meanSCOffset = selectedSCSectorOffsetDistribution->GetMean();
        if (verbose){
            cout << "Dumping SC results...\n=======================================" << endl;
            cout << "SC mean Offset = " << meanSCOffset << endl;
            cout << "TDC Offsets" << endl;
            cout << "Sector\toldValue\tValueToUse\tmeanOffset\tTotal" << endl;
        }
        outFile.open(prefix + "sc_tdc_timing_offsets.txt");
        for (int sector = 1; sector <= 30; sector++){
            outFile << sc_tdc_time_offsets[sector-1] + selectedSCSectorOffset->GetBinContent(sector) - meanSCOffset << endl;
            if (verbose) printf("%i\t%.3f\t\t%.3f\t\t%.3f\t\t%.3f\n",sector, sc_tdc_time_offsets[sector-1], selectedSCSectorOffset->GetBinContent(sector), meanSCOffset,
                    sc_tdc_time_offsets[sector-1] + selectedSCSectorOffset->GetBinContent(sector) - meanSCOffset);
        }
        outFile.close();
        if (verbose){
            cout << "ADC Offsets" << endl;
            cout << "Sector\tvalueToUse\toldValue\tmeanOffset\tTotal" << endl;
        }
        outFile.open(prefix + "sc_adc_timing_offsets.txt");
        for (int sector = 1; sector <= 30; sector++){
            outFile << sc_fadc_time_offsets[sector-1] + selectedSCSectorOffset->GetBinContent(sector) - meanSCOffset << endl;
            if (verbose) printf("%i\t%.3f\t\t%.3f\t\t%.3f\t\t%.3f\n",sector,sc_fadc_time_offsets[sector-1], selectedSCSectorOffset->GetBinContent(sector), meanSCOffset,
                    sc_fadc_time_offsets[sector-1] + selectedSCSectorOffset->GetBinContent(sector) - meanSCOffset);
        }
        outFile.close();

        outFile.open(prefix + "sc_base_time.txt");
        outFile << sc_t_base_fadc - meanSCOffset << " " << sc_t_base_tdc - meanSCOffset << endl;
        if (verbose) {
            printf("SC ADC Base = %f - (%f) = %f\n", sc_t_base_fadc, meanSCOffset, sc_t_base_fadc - meanSCOffset);
            printf("SC TDC Base = %f - (%f) = %f\n", sc_t_base_tdc, meanSCOffset, sc_t_base_tdc - meanSCOffset);
        }
        outFile.close();
    }

    TH1I *this1DHist = Get1DHistogram("HLDetectorTiming", "TRACKING", "TOF - SC Target Time");
    if(this1DHist != NULL){
        //Gaussian
        Double_t maximum = this1DHist->GetBinCenter(this1DHist->GetMaximumBin());
        TFitResultPtr fr = this1DHist->Fit("gaus", "S", "", maximum - 1.5, maximum + 1.5);
        float mean = fr->Parameter(1);
        outFile.open(prefix + "tof_base_time.txt");
        if (verbose) {
            printf("TOF ADC Base = %f - (%f) - (%f) = %f\n", tof_t_base_fadc, mean, meanSCOffset, tof_t_base_fadc - mean - meanSCOffset);
            printf("TOF TDC Base = %f - (%f) - (%f) = %f\n", tof_t_base_tdc, mean, meanSCOffset, tof_t_base_tdc - mean - meanSCOffset);
        }
        outFile << tof_t_base_fadc - mean - meanSCOffset<< " " << tof_t_base_tdc - mean - meanSCOffset<< endl;
        outFile.close();
    }

    this1DHist = Get1DHistogram("HLDetectorTiming", "TRACKING", "BCAL - SC Target Time");
    if(this1DHist != NULL){
        //Gaussian
        Double_t maximum = this1DHist->GetBinCenter(this1DHist->GetMaximumBin());
        TFitResultPtr fr = this1DHist->Fit("gaus", "S", "", maximum - 5, maximum + 5);
        float mean = fr->Parameter(1);
        outFile.open(prefix + "bcal_base_time.txt");
        if (verbose) {
            printf("BCAL ADC Base = %f - (%f) - (%f) = %f\n", bcal_t_base_fadc, mean, meanSCOffset, bcal_t_base_fadc - mean - meanSCOffset);
            printf("BCAL TDC Base = %f - (%f) - (%f) = %f\n", bcal_t_base_tdc, mean, meanSCOffset, bcal_t_base_tdc - mean - meanSCOffset);
        }
        outFile << bcal_t_base_fadc - mean - meanSCOffset << " " << bcal_t_base_tdc - mean - meanSCOffset << endl; // TDC info not used
        outFile.close();
    }

    this1DHist = Get1DHistogram("HLDetectorTiming", "TRACKING", "FCAL - SC Target Time");
    if(this1DHist != NULL){
        //Gaussian
        Double_t maximum = this1DHist->GetBinCenter(this1DHist->GetMaximumBin());
        TFitResultPtr fr = this1DHist->Fit("gaus", "S", "", maximum - 5, maximum + 5);
        float mean = fr->Parameter(1);
        outFile.open(prefix + "fcal_base_time.txt");
        if (verbose) {
            printf("BCAL ADC Base = %f - (%f) - (%f) = %f\n",fcal_t_base, mean, meanSCOffset, fcal_t_base - mean - meanSCOffset);
        }
        outFile << fcal_t_base - mean - meanSCOffset<< endl; 
        outFile.close();
    }

    this1DHist = Get1DHistogram("HLDetectorTiming", "TRACKING", "Earliest CDC Time Minus Matched SC Time");
    if(this1DHist != NULL){
        //Gaussian
        Double_t maximum = this1DHist->GetBinCenter(this1DHist->GetMaximumBin());
        TFitResultPtr fr = this1DHist->Fit("gaus", "S", "", maximum - 15, maximum + 10);
        float mean = fr->Parameter(1);
        outFile.open(prefix + "cdc_base_time.txt");
        if (verbose) {
            printf("BCAL ADC Base = %f - (%f) - (%f) = %f\n",cdc_t_base, mean, meanSCOffset, cdc_t_base - mean - meanSCOffset);
        }
        outFile << cdc_t_base - mean - meanSCOffset << endl;
        outFile.close();
    }
    thisFile->Write();
    return;
    }
