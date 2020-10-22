// code to plot profiles of the rapidity relations
// obtained in bins of pT/M, as a function of y
// also stored in ROOT file

void plot_y()
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
  string type = "rho2delta0/";
  
  // choose state we're plotting
  const int n_states = 3;
  string dataNames[n_states] = {"jpsi", "psi2", "ups1"};   // "jpsi", "psi2", "ups1", "ups2", "ups3"
  // choose sqrt(s) we're plotting
  const int n_sqs = 1;
  string sqsNames[n_sqs] = {"7"};   // "7", "13"

  /////////////////////////////////////////
  // part 2 : vars defined by above or constant  
  
  for(int i_s = 0; i_s < n_states; i_s++) {
    string dataName = dataNames[i_s];
    cout << "running state " << dataName << endl;

    string fName = Form("%s%sMC_res_%s", loc.c_str(), type.c_str(), dataName.c_str());   // file to open

    TFile *outfile = new TFile(Form("%sy_profs/plots_%s.root", type.c_str(), dataName.c_str()), "RECREATE");
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
	Double_t y, y_gg; // to split in y bins
	Double_t w_gg, w_cos; // plotting weight - fixed
	Double_t w_HX_tr, w_ggHX_tr, w_HX_lg, w_ggHX_lg; // plotting weight - polarizations

	int xi_bins = 4;
	double xilims[xi_bins+1];
	for(int i = 0; i <= xi_bins; i++)
	  xilims[i] = i*10.;
	xilims[0] = 1.;
	
	int y_bins = 36;
	double ylims[y_bins+1];
	for(int i = 0; i <= y_bins; i++)
	  ylims[i] = -4.5+0.25*i;;
	
	string fr[5] = {"HX", "ggHX", "HX", "ggHX", ""};
	string pol[5] = {"tr", "tr", "lg", "lg", "unp"};

	TProfile ***y_prof = new TProfile**[5];
	for(int i = 0; i < 5; i++) {
	  y_prof[i] = new TProfile*[xi_bins];
	  for(int j = 0; j < xi_bins; j++) {
	    y_prof[i][j] = new TProfile(Form("%s_y_prof%d_xi%d", dataName.c_str(), i, j), Form("%s TeV %.0f < #xi < %.0f #hat{y}/y %s %s", sqsName.c_str(), xilims[j], xilims[j+1], fr[i].c_str(), pol[i].c_str()), y_bins, ylims, 0, 0.2);
	    y_prof[i][j]->Sumw2();
	  }
	}
	
	Int_t nentries;
	
	TFile *fin = new TFile(Form("%s_beta2.root", fName.c_str()));
	TTree *tree = (TTree*)fin->Get(ssName[sqsName].c_str());
	
	tree->SetBranchAddress("xi", &xi);
	tree->SetBranchAddress("y", &y);
	tree->SetBranchAddress("y_gg", &y_gg);
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
	    for (int k = 0; k < xi_bins; k++) 
	      if(xi > xilims[k] && xi < xilims[k+1]) {
		y_prof[0][k]->Fill(y, (y-y_gg)/y, w_gg*w_cos*w_HX_tr);
		
		y_prof[1][k]->Fill(y, (y-y_gg)/y, w_gg*w_cos*w_ggHX_tr);
		
		y_prof[2][k]->Fill(y, (y-y_gg)/y, w_gg*w_cos*w_HX_lg);
		
		y_prof[3][k]->Fill(y, (y-y_gg)/y, w_gg*w_cos*w_ggHX_lg);

		y_prof[4][k]->Fill(y, (y-y_gg)/y, w_gg*w_cos);
	      }
	    if((i+1)%chk == 0) {
	      cout << (i+1)/chk << "% | " << flush;
	    }
	  }
	cout << endl;
	fin->Close();
      
	
	/////////////////////////////////////////
	// part 5 : plotting

	TFile *outfile2 = new TFile(Form("%sy_profs/plots_%s.root", type.c_str(), dataName.c_str()), "UPDATE");
	
	TCanvas *can = new TCanvas("", "", 700, 700);
	
	for(int j = 0; j < 5; j++) {
	  for(int i = 0; i < xi_bins; i++) {
	    int color = (i==4? i+2 : i+1);
	    
	    y_prof[j][i]->SetLineColor(color);
	    y_prof[j][i]->SetMarkerColor(color);
	    
	    //y_prof[j][i]->SetLineStyle(j%2 == 0 ? kSolid : kDashed);
	    
	    y_prof[j][i]->SetStats(0);
	    y_prof[j][i]->Write();
	    
	  }
	}
	outfile2->Close();
	
      }
  }
}

