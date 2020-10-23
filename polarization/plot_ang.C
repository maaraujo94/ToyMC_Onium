// code to plot angular distributions in bins of y, for all pT

void plot_ang()
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
  const int n_sqs = 2;
  string sqsNames[n_sqs] = {"7", "13"};   // "7", "13"

  /////////////////////////////////////////
  // part 2 : vars defined by above or constant  
  
  for(int i_s = 0; i_s < n_states; i_s++) {
    string dataName = dataNames[i_s];
    cout << "running state " << dataName << endl;

    string fName = Form("%s%sMC_res_%s", loc.c_str(), type.c_str(), dataName.c_str());   // file to open

    TFile *outfile = new TFile(Form("%sang_plots/plots_%s.root", type.c_str(), dataName.c_str()), "RECREATE");
    outfile->Close();
  
    for(int i_sqs = 0; i_sqs < n_sqs; i_sqs++)
      {
	string sqsName = sqsNames[i_sqs];
	cout << "running sqrt(s) = " << sqsName << " TeV" << endl;
	
	double mass = stateMass[dataName];
      
	/////////////////////////////////////////
	// part 4 : making the MC histograms
	
	// variables to open the MC
	Double_t costh_HX, phi_HX, phith_HX; // angular variables
	Double_t xi; // to plot lambda(pT)
	Double_t y; // to split in y bins
	Double_t w_gg, w_cos; // plotting weight - fixed
	Double_t w_HX_tr, w_ggHX_tr, w_HX_lg, w_ggHX_lg; // plotting weight - polarizations

	int y_bins = 5;
	int pt_bins = 10;
	double ylims[y_bins+1];
	for(int i = 0; i <= y_bins; i++)
	  ylims[i] = 1.*i;

	TH1F ***costh = new TH1F**[4];
	TH1F ***phi = new TH1F**[4];
	TH1F ***phi_t = new TH1F**[4];
	for(int i = 0; i < 4; i++) {
	  costh[i] = new TH1F*[y_bins];
	  phi[i] = new TH1F*[y_bins];
	  phi_t[i] = new TH1F*[y_bins];
	  for(int j = 0; j < y_bins; j++) {
	    costh[i][j] = new TH1F(Form("%s%s_ct%d_y%d", dataName.c_str(), sqsName.c_str(), i, j), "cos#theta_{HX}", 100, -1, 1);
	    phi[i][j] = new TH1F(Form("%s%s_phi%d_y%d", dataName.c_str(), sqsName.c_str(), i, j), "#phi_{HX}", 100, -180, 180);
	    phi_t[i][j] = new TH1F(Form("%s%s_pt%d_y%d", dataName.c_str(), sqsName.c_str(), i, j), "#tilde{phi}_{HX}", 100, -180, 180);
	  }
	}
	
	Int_t nentries;
	
	TFile *fin = new TFile(Form("%s_beta2.root", fName.c_str()));
	TTree *tree = (TTree*)fin->Get(ssName[sqsName].c_str());
	
	tree->SetBranchAddress("costh_HX", &costh_HX);
	tree->SetBranchAddress("phi_HX", &phi_HX);
	tree->SetBranchAddress("phith_HX", &phith_HX);
	tree->SetBranchAddress("xi", &xi);
	tree->SetBranchAddress("y", &y);
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
	    for (int k = 0; k < y_bins; k++) 
	      if(abs(y) > ylims[k] && abs(y) < ylims[k+1] && xi > 1 && xi < 50) {
		costh[0][k]->Fill(costh_HX, w_gg*w_cos*w_HX_tr);
		phi[0][k]->Fill(phi_HX, w_gg*w_cos*w_HX_tr);
		phi_t[0][k]->Fill(phith_HX, w_gg*w_cos*w_HX_tr);
		
		costh[1][k]->Fill(costh_HX, w_gg*w_cos*w_ggHX_tr);
		phi[1][k]->Fill(phi_HX, w_gg*w_cos*w_ggHX_tr);
		phi_t[1][k]->Fill(phith_HX, w_gg*w_cos*w_ggHX_tr);
		
		costh[2][k]->Fill(costh_HX, w_gg*w_cos*w_HX_lg);
		phi[2][k]->Fill(phi_HX, w_gg*w_cos*w_HX_lg);
		phi_t[2][k]->Fill(phith_HX, w_gg*w_cos*w_HX_lg);
		
		costh[3][k]->Fill(costh_HX, w_gg*w_cos*w_ggHX_lg);
		phi[3][k]->Fill(phi_HX, w_gg*w_cos*w_ggHX_lg);
		phi_t[3][k]->Fill(phith_HX, w_gg*w_cos*w_ggHX_lg);
	      }
	    if((i+1)%chk == 0) {
	      cout << (i+1)/chk << "% | " << flush;
	    }
	  }
	cout << endl;
	fin->Close();
      
	
	/////////////////////////////////////////
	// part 5 : plotting

	TFile *outfile2 = new TFile(Form("%sang_plots/plots_%s.root", type.c_str(), dataName.c_str()), "UPDATE");
	
	TCanvas *can = new TCanvas("", "", 700, 700);
	string fr[4] = {"HX", "ggHX", "HX", "ggHX"};
	string pol[4] = {"tr", "tr", "lg", "lg"};
	
	for(int j = 0; j < 4; j++) {
	  for(int i = 0; i < y_bins; i++) {
	    costh[j][i]->SetLineColor(i+1);
	    phi[j][i]->SetLineColor(i+1);
	    phi_t[j][i]->SetLineColor(i+1);
	    
	    costh[j][i]->SetLineStyle(j%2 == 0 ? kSolid : kDashed);
	    phi[j][i]->SetLineStyle(j%2 == 0 ? kSolid : kDashed);
	    phi_t[j][i]->SetLineStyle(j%2 == 0 ? kSolid : kDashed);
	    
	    costh[j][i]->Write();
	    phi[j][i]->Write();
	    phi_t[j][i]->Write();
	    
	    costh[j][i]->SetMinimum(0);
	    costh[j][i]->SetMaximum(48e3);
	    phi[j][i]->SetMinimum(0);
	    phi[j][i]->SetMaximum(35e3);
	    phi_t[j][i]->SetMinimum(0);
	    phi_t[j][i]->SetMaximum(36e3);
	  }
	}
	outfile2->Close();

	// draw the costh plots
	costh[0][0]->Draw("hist");
	for(int i = 0; i < y_bins; i++) {
	  if(i!=0) costh[0][i]->Draw("hist same");
	  costh[1][i]->Draw("hist same");
	}
	
	can->SaveAs(Form("%sang_plots/%s_costh_trans_%s.pdf", type.c_str(), dataName.c_str(), sqsName.c_str()));
	can->Clear();
	
	costh[2][0]->Draw("hist");
	for(int i = 0; i < y_bins; i++) {
	  if(i!=0) costh[2][i]->Draw("hist same");
	  costh[3][i]->Draw("hist same");
	}
	
	can->SaveAs(Form("%sang_plots/%s_costh_long_%s.pdf", type.c_str(), dataName.c_str(), sqsName.c_str()));
	can->Clear();
	
	// draw the phi plots
	phi[0][0]->Draw("hist");
	for(int i = 0; i < y_bins; i++) {
	  if(i!=0) phi[0][i]->Draw("hist same");
	  phi[1][i]->Draw("hist same");
	}
	
	can->SaveAs(Form("%sang_plots/%s_phi_trans_%s.pdf", type.c_str(), dataName.c_str(), sqsName.c_str()));
	can->Clear();
	
	phi[2][0]->Draw("hist");
	for(int i = 0; i < y_bins; i++) {
	  if(i!=0) phi[2][i]->Draw("hist same");
	  phi[3][i]->Draw("hist same");
	}
	
	can->SaveAs(Form("%sang_plots/%s_phi_long_%s.pdf", type.c_str(), dataName.c_str(), sqsName.c_str()));
	can->Clear();
	
	// draw the phi_t plots
	phi_t[0][0]->Draw("hist");
	for(int i = 0; i < y_bins; i++) {
	  if(i!=0) phi_t[0][i]->Draw("hist same");
	  phi_t[1][i]->Draw("hist same");
	}
	
	can->SaveAs(Form("%sang_plots/%s_phi_t_trans_%s.pdf", type.c_str(), dataName.c_str(), sqsName.c_str()));
	can->Clear();
	
	phi_t[2][0]->Draw("hist");
	for(int i = 0; i < y_bins; i++) {
	  if(i!=0) phi_t[2][i]->Draw("hist same");
	  phi_t[3][i]->Draw("hist same");
	}
	
	can->SaveAs(Form("%sang_plots/%s_phi_t_long_%s.pdf", type.c_str(), dataName.c_str(), sqsName.c_str()));
	can->Clear();
      }
  }
}

