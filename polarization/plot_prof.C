// code to plot profiles of the necessary angular quantities for calculating lambdas
// obtained in bins of y, as a function of pT/M
// also stored in ROOT file

void plot_prof()
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
  string type = "comp_j/";
  
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

    TFile *outfile = new TFile(Form("%sang_profs/plots_%s.root", type.c_str(), dataName.c_str()), "RECREATE");
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

	int y_bins = 4;
	double ylims[y_bins+1];
	for(int i = 0; i <= y_bins; i++)
	  ylims[i] = i;

	int pt_bins = 20;
	double ptlims[pt_bins+1];
	for(int i = 0; i < 9; i++)
	  ptlims[i] = 1.+i;
	for(int i = 0; i < 5; i++)
	  ptlims[i+9] = 10.+2.*i;
	for(int i = 0; i < 4; i++)
	  ptlims[i+14] = 20+2.5*i;
	for(int i = 0; i <= 2; i++)
	  ptlims[i+18] = 30+5.*i;

	string fr[4] = {"HX", "ggHX", "HX", "ggHX"};
	string pol[4] = {"tr", "tr", "lg", "lg"};

	TProfile ***costh = new TProfile**[4];
	TProfile ***phi = new TProfile**[4];
	TProfile ***phith = new TProfile**[4];
	for(int i = 0; i < 4; i++) {
	  costh[i] = new TProfile*[y_bins];
	  phi[i] = new TProfile*[y_bins];
	  phith[i] = new TProfile*[y_bins];
	  for(int j = 0; j < y_bins; j++) {
	    costh[i][j] = new TProfile(Form("h_ct%d_y%d", i, j), Form("%s TeV %.1f < |y| < %.1f cos^{2}#theta_{HX} %s %s", sqsName.c_str(), ylims[j], ylims[j+1], fr[i].c_str(), pol[i].c_str()), pt_bins, ptlims, 0, 1);
	    costh[i][j]->Sumw2();
	    phi[i][j] = new TProfile(Form("h_phi%d_y%d", i, j), Form("%s TeV %.1f < |y| < %.1f cos(2#phi_{HX}) %s %s", sqsName.c_str(), ylims[j], ylims[j+1], fr[i].c_str(), pol[i].c_str()), pt_bins, ptlims, -1, 1);
	    phi[i][j]->Sumw2();
	    phith[i][j] = new TProfile(Form("h_pt%d_y%d", i, j), Form("%s TeV %.1f < |y| < %.1f sin(2#theta_{HX})cos#phi_{HX} %s %s", sqsName.c_str(), ylims[j], ylims[j+1], fr[i].c_str(), pol[i].c_str()), pt_bins, ptlims, -1, 1);
	    phith[i][j]->Sumw2();
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
		costh[0][k]->Fill(xi, costh_HX*costh_HX, w_gg*w_cos*w_HX_tr);
		phi[0][k]->Fill(xi, cos(2.*phi_HX), w_gg*w_cos*w_HX_tr);
		phith[0][k]->Fill(xi, 2.*sqrt(1-costh_HX*costh_HX)*costh_HX*cos(phi_HX), w_gg*w_cos*w_HX_tr);
		
		costh[1][k]->Fill(xi, costh_HX*costh_HX, w_gg*w_cos*w_ggHX_tr);
		phi[1][k]->Fill(xi, cos(2.*phi_HX), w_gg*w_cos*w_ggHX_tr);
		phith[1][k]->Fill(xi, 2.*sqrt(1-costh_HX*costh_HX)*costh_HX*cos(phi_HX), w_gg*w_cos*w_ggHX_tr);
		
		costh[2][k]->Fill(xi, costh_HX*costh_HX, w_gg*w_cos*w_HX_lg);
		phi[2][k]->Fill(xi, cos(2.*phi_HX), w_gg*w_cos*w_HX_lg);
		phith[2][k]->Fill(xi, 2.*sqrt(1-costh_HX*costh_HX)*costh_HX*cos(phi_HX), w_gg*w_cos*w_HX_lg);
		
		costh[3][k]->Fill(xi, costh_HX*costh_HX, w_gg*w_cos*w_ggHX_lg);
		phi[3][k]->Fill(xi, cos(2.*phi_HX), w_gg*w_cos*w_ggHX_lg);
		phith[3][k]->Fill(xi, 2.*sqrt(1-costh_HX*costh_HX)*costh_HX*cos(phi_HX), w_gg*w_cos*w_ggHX_lg);
	      }
	    if((i+1)%chk == 0) {
	      cout << (i+1)/chk << "% | " << flush;
	    }
	  }
	cout << endl;
	fin->Close();
      
	
	/////////////////////////////////////////
	// part 5 : plotting

	TFile *outfile2 = new TFile(Form("%sang_profs/plots_%s.root", type.c_str(), dataName.c_str()), "UPDATE");
	
	TCanvas *can = new TCanvas("", "", 700, 700);
	
	for(int j = 0; j < 4; j++) {
	  for(int i = 0; i < y_bins; i++) {
	    int color = (i==4? i+2 : i+1);
	    
	    costh[j][i]->SetLineColor(color);
	    phi[j][i]->SetLineColor(color);
	    phith[j][i]->SetLineColor(color);
	    costh[j][i]->SetMarkerColor(color);
	    phi[j][i]->SetMarkerColor(color);
	    phith[j][i]->SetMarkerColor(color);
	    
	    costh[j][i]->SetLineStyle(j%2 == 0 ? kSolid : kDashed);
	    phi[j][i]->SetLineStyle(j%2 == 0 ? kSolid : kDashed);
	    phith[j][i]->SetLineStyle(j%2 == 0 ? kSolid : kDashed);
	    
	    costh[j][i]->Write();
	    phi[j][i]->Write();
	    phith[j][i]->Write();

	    costh[j][i]->GetYaxis()->SetRangeUser(0,1);
	    phi[j][i]->GetYaxis()->SetRangeUser(-1,1);
	    phith[j][i]->GetYaxis()->SetRangeUser(-1,1);
	    
	  }
	}
	outfile2->Close();

	// draw the costh plots
	costh[0][0]->Draw("error");
	for(int i = 0; i < y_bins; i++) {
	  if(i!=0) costh[0][i]->Draw("same");
	  costh[1][i]->Draw("same");
	}
	TLine *ct1 = new TLine(1, 0.2, 40, 0.2);
	ct1->SetLineStyle(kDotted);
	ct1->Draw();
	TLine *ct2 = new TLine(1, 0.4, 40, 0.4);
	ct2->SetLineStyle(kDotted);
	ct2->Draw();
	
	can->SaveAs(Form("%sang_profs/%s_cos2th_trans_%s.pdf", type.c_str(), dataName.c_str(), sqsName.c_str()));
	can->Clear();
	
	costh[2][0]->Draw("error");
	for(int i = 0; i < y_bins; i++) {
	  if(i!=0) costh[2][i]->Draw("same");
	  costh[3][i]->Draw("same");
	}
	ct1->Draw();
	ct2->Draw();
	
	can->SaveAs(Form("%sang_profs/%s_cos2th_long_%s.pdf", type.c_str(), dataName.c_str(), sqsName.c_str()));
	can->Clear();
	
	// draw the phi plots
	phi[0][0]->Draw("error");
	for(int i = 0; i < y_bins; i++) {
	  if(i!=0) phi[0][i]->Draw("same");
	  phi[1][i]->Draw("same");
	}
	
	TLine *ph1 = new TLine(1, -0.5, 40, -0.5);
	ph1->SetLineStyle(kDotted);
	ph1->Draw();
	TLine *ph2 = new TLine(1, 0.5, 40, 0.5);
	ph2->SetLineStyle(kDotted);
	ph2->Draw();

	can->SaveAs(Form("%sang_profs/%s_cos2phi_trans_%s.pdf", type.c_str(), dataName.c_str(), sqsName.c_str()));
	can->Clear();
	
	phi[2][0]->Draw("error");
	for(int i = 0; i < y_bins; i++) {
	  if(i!=0) phi[2][i]->Draw("same");
	  phi[3][i]->Draw("same");
	}
	
	ph1->Draw();
	ph2->Draw();

	can->SaveAs(Form("%sang_profs/%s_cos2phi_long_%s.pdf", type.c_str(), dataName.c_str(), sqsName.c_str()));
	can->Clear();
	
	// draw the phith plots
	phith[0][0]->Draw("error");
	for(int i = 0; i < y_bins; i++) {
	  if(i!=0) phith[0][i]->Draw("same");
	  phith[1][i]->Draw("same");
	}

	TLine *pt1 = new TLine(1, -0.4, 40, -0.4);
	pt1->SetLineStyle(kDotted);
	pt1->Draw();
	TLine *pt2 = new TLine(1, 0.4, 40, 0.4);
	pt2->SetLineStyle(kDotted);
	pt2->Draw();

	can->SaveAs(Form("%sang_profs/%s_phitheta_trans_%s.pdf", type.c_str(), dataName.c_str(), sqsName.c_str()));
	can->Clear();
	
	phith[2][0]->Draw("error");
	for(int i = 0; i < y_bins; i++) {
	  if(i!=0) phith[2][i]->Draw("same");
	  phith[3][i]->Draw("same");
	}

	pt1->Draw();
	pt2->Draw();

	
	can->SaveAs(Form("%sang_profs/%s_phitheta_long_%s.pdf", type.c_str(), dataName.c_str(), sqsName.c_str()));
	can->Clear();
      }
  }
}

