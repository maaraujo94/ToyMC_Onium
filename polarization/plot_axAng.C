// code to plot profiles of the rapidity relations
// obtained in bins of pT/M, as a function of y
// also stored in ROOT file

void plot_axAng()
{
  /////////////////////////////////////////
  // part 0 : auxiliary variables
  
  string loc = "/home/mariana/Documents/2020_PhD_work/Phenom/ToyMC_Onium/MC/";
  
  // associate the chosen state and sqrt(s) with its mass and tree name
  map<string, double> stateMass;
  stateMass["jpsi"] = 3.097;
  stateMass["psi2"] = 3.686;
  stateMass["ups1"] = 9.46;
  stateMass["ups2"] = 10.023;
  stateMass["ups3"] = 10.355;
  map<string, string> ssName;
  ssName["7"] = "qqbar_0";
  ssName["13"] = "qqbar_1";
    
  /////////////////////////////////////////
  // part 1 : variables to be changed on each run

  // type of MC files
  string type = "angle_kT/";
  
  // choose state we're plotting
  const int n_states = 1;
  string dataNames[n_states] = {"jpsi"};   // "jpsi", "psi2", "ups1", "ups2", "ups3"
  // choose sqrt(s) we're plotting
  const int n_sqs = 1;
  string sqsNames[n_sqs] = {"7"};   // "7", "13"

  /////////////////////////////////////////
  // part 2 : vars defined by above or constant  
  
  for(int i_s = 0; i_s < n_states; i_s++) {
    string dataName = dataNames[i_s];
    cout << "running state " << dataName << endl;

    string fName = Form("%s%sMC_res_%s", loc.c_str(), type.c_str(), dataName.c_str());   // file to open

    TFile *outfile = new TFile(Form("%saxAng_plots/plots_%s.root", type.c_str(), dataName.c_str()), "RECREATE");
    outfile->Close();
  
    for(int i_sqs = 0; i_sqs < n_sqs; i_sqs++)
      {
	string sqsName = sqsNames[i_sqs];
	cout << "running sqrt(s) = " << sqsName << " TeV" << endl;
	
	double mass = stateMass[dataName];
      
	/////////////////////////////////////////
	// part 4 : making the MC histograms
	
	// variables to open the MC
	Double_t xi; // to split in pT/M bins
	Double_t y; // to split in y bins
	Double_t w_gg, w_cos; // plotting weight - fixed
	Double_t w_HX_tr, w_ggHX_tr, w_HX_lg, w_ggHX_lg; // plotting weight - polarizations
	Double_t axAng;

	int y_bins = 8;
	double ylims[y_bins+1];
	for(int i = 0; i <= y_bins; i++)
	  ylims[i] = i*0.5;

	int pt_bins = 5;
	double ptlims[pt_bins+1];
	for(int i = 0; i <= pt_bins; i++)
	  ptlims[i] = 1.+i;
	
	string fr[4] = {"HX", "ggHX", "HX", "ggHX"};
	string pol[4] = {"tr", "tr", "lg", "lg"};

	TProfile ***ang_prof = new TProfile**[4];
	for(int i = 0; i < 4; i++) {
	  ang_prof[i] = new TProfile*[pt_bins];
	  for(int j = 0; j < pt_bins; j++) {
	    ang_prof[i][j] = new TProfile(Form("%s%s_ang_prof%d_y%d", dataName.c_str(), sqsName.c_str(), i, j), Form("%s TeV %.0f < #xi < %.0f HH-ggHX angle %s %s", sqsName.c_str(), ptlims[j], ptlims[j+1], fr[i].c_str(), pol[i].c_str()), y_bins, ylims, 0, 1.);
	    ang_prof[i][j]->Sumw2();
	  }
	}
	
	Int_t nentries;
	
	TFile *fin = new TFile(Form("%s_beta2.root", fName.c_str()));
	TTree *tree = (TTree*)fin->Get(ssName[sqsName].c_str());
	
	tree->SetBranchAddress("xi", &xi);
	tree->SetBranchAddress("y", &y);
	tree->SetBranchAddress("axAngle", &axAng);
	tree->SetBranchAddress("w_gg", &w_gg);
	tree->SetBranchAddress("w_cos", &w_cos);
	tree->SetBranchAddress("w_HX_tr", &w_HX_tr);
	tree->SetBranchAddress("w_ggHX_tr", &w_ggHX_tr);
	tree->SetBranchAddress("w_HX_lg", &w_HX_lg);
	tree->SetBranchAddress("w_ggHX_lg", &w_ggHX_lg);
	
	nentries = (Int_t)tree->GetEntries();
	cout << nentries << " events in tree " << endl;
	int chk = nentries / 100;
	
	for( Int_t i = 0; i < nentries; i++)
	  {
	    tree->GetEntry(i);
	    
	    // xi part (forward y)
	    for (int k = 0; k < pt_bins; k++) 
	      if(xi > ptlims[k] && xi < ptlims[k+1]) {
		ang_prof[0][k]->Fill(abs(y), axAng, w_gg*w_cos*w_HX_tr);
		
		ang_prof[1][k]->Fill(abs(y), axAng, w_gg*w_cos*w_ggHX_tr);
		
		ang_prof[2][k]->Fill(abs(y), axAng, w_gg*w_cos*w_HX_lg);
		
		ang_prof[3][k]->Fill(abs(y), axAng, w_gg*w_cos*w_ggHX_lg);
	      }
	    if((i+1)%chk == 0) {
	      cout << (i+1)/chk << "% | " << flush;
	    }
	  }
	cout << endl;
	fin->Close();
      
	
	/////////////////////////////////////////
	// part 5 : plotting

	TFile *outfile2 = new TFile(Form("%saxAng_plots/plots_%s.root", type.c_str(), dataName.c_str()), "UPDATE");
	
	TCanvas *can = new TCanvas("", "", 700, 700);
	
	for(int j = 0; j < 4; j++) {
	  for(int i = 0; i < pt_bins; i++) {
	    int color = (i==4? i+2 : i+1);
	    
	    ang_prof[j][i]->SetLineColor(color);
	    ang_prof[j][i]->SetMarkerColor(color);
	    
	    ang_prof[j][i]->SetLineStyle(j%2 == 0 ? kSolid : kDashed);
	    
	    ang_prof[j][i]->SetStats(0);
	    ang_prof[j][i]->Write();
	    
	  }
	}
	outfile2->Close();
	
      }
  }
}

