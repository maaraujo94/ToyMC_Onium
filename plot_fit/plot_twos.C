/* code to plot toy MC histos vs data graphs
   - make TGraphs (external code)
   - make THistos
   - plot both together
   process is repeated for both xi plots in y bins and y plots in xi bins
*/

#import "plot_codes.C"
#import "pull_codes.C"

void plot_twos()
{
  /////////////////////////////////////////
  // part 0 : auxiliary variables
  
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
  
  double y_min_ups[7] = {0, 0.6, 2, 2.5, 3, 3.5, 4};
  double y_max_ups[7] = {0.6, 1.2, 2.5, 3, 3.5, 4, 4.5};
  double y_min_psi[9] = {0, 0.3, 0.6, 0.9, 2, 2.5, 3, 3.5, 4};
  double y_max_psi[9] = {0.3, 0.6, 0.9, 1.2, 2.5, 3, 3.5, 4, 4.5};

  double xi_min_global = 1;
  
  /////////////////////////////////////////
  // part 1 : variables to be changed on each run
  
  // choose state we're plotting
  string dataName = "jpsi";   // "jpsi", "psi2", "ups1", "ups2", "ups3"
  // choose sqrt(s) we're plotting
  string sqsNames[2] = {"7", "13"};   // "7", "13"

  string fName = Form("MC_res_%s", dataName.c_str());   // file to open
  double xi_bin_width = 0.1;   // constant for easier normalization

  double y_plot_min = 0, y_plot_max = 5;
  
  /////////////////////////////////////////
  // part 2 : vars defined by above or constant

  double normFactor_b[2] = {1.e4, 1.e4}; // normalization to be applied on MC

  if(dataName == "jpsi") {
    normFactor_b[0] = 5.e3;
    normFactor_b[1] = 4.e4;
  }
  else if(dataName == "psi2") {
    normFactor_b[0] = 2e3;
    normFactor_b[1] = 1.e4;
  }
  else if(dataName == "ups1") {
    normFactor_b[0] = 2.e3;
    normFactor_b[1] = 9.e3;
  }
  else if(dataName == "ups2") {
    normFactor_b[0] = 1.4e3;
    normFactor_b[1] = 5.e3;
  }
  else if(dataName == "ups3") {
    normFactor_b[0] = 1.e3;
    normFactor_b[1] = 2.5e3;
  }

  // extra part: create ROOT file and save results

  TFile *tf = new TFile(Form("MC_vs_Data_%s.root", dataName.c_str()), "RECREATE", "MC histograms and Data graphs");
  TTree *norm_factors = new TTree("norms", "norm factors");
  
  TBranch *nf_b1 = norm_factors->Branch("normFactor_b1", &normFactor_b[0]);
  TBranch *nf_b2 = norm_factors->Branch("normFactor_b2", &normFactor_b[1]);

  int xi_n = 1, y_n = 1;
  
  norm_factors->Fill();
  tf->Write();
  tf->Close();
  
  for(int i_sqs = 0; i_sqs < 2; i_sqs++)
    {
      string sqsName = sqsNames[i_sqs];

      // determine xi_n from the data used
      if(dataName == "jpsi" || dataName == "psi2")
	xi_n = 9;
      else if(dataName == "ups1" || dataName == "ups2" || dataName == "ups3")
	xi_n = 7;

      // determine y_n from the data used (for now same in both, but keep separate bc it might change
      if(sqsName == "7") {
	if(dataName == "jpsi") y_n = 7;
	else if(dataName == "psi2") y_n = 6;
	else if(dataName == "ups1" || dataName == "ups2" || dataName == "ups3") y_n = 5;
      }
      else if(sqsName == "13") { // just doing LHCb bc no intersect
	if(dataName == "jpsi") y_n = 7;
	else if(dataName == "psi2") y_n = 6;
	else if(dataName == "ups1" || dataName == "ups2" || dataName == "ups3") y_n = 5;
      }
      
      // data TGraphs (xi and y dists)
      TGraphAsymmErrors **xi_g = new TGraphAsymmErrors*[xi_n];
      TGraphAsymmErrors **y_g = new TGraphAsymmErrors*[y_n];

      TH1F ***xi_h = new TH1F**[2];
      for(int i = 0; i < 2; i++) xi_h[i] = new TH1F*[xi_n]; // xi histogram, defined before the data cycle
      double y_min[xi_n], y_max[xi_n]; // y bins of each xi dist
      double xi_min[xi_n], xi_max[xi_n]; // limits for plotting the xi dists (x-axis)
      double xi_yax_min[xi_n], xi_yax_max[xi_n]; // limits for plotting the xi dists (y-axis)

      double y_lim[xi_n+2]; // binning for the y dist
      double y_xi_min[y_n], y_xi_max[y_n]; // pT bins for the y dists
      double y_yax_min_lin, y_yax_max_lin; // limits for plotting the y dists (y-axis, lin)
      double y_yax_min_log, y_yax_max_log; // limits for plotting the y dists (y-axis, log)
      
      double mass = stateMass[dataName];

      /////////////////////////////////////////
      // part 3 : making the data TGraphs

      // y limits are the same for both psi states in both sqrt(s)
      if(dataName == "jpsi" || dataName == "psi2") {
	y_lim[0] = y_min_psi[0];
	for(int i = 0; i < 4; i++)
	  y_lim[i+1] = y_max_psi[i];
	y_lim[5] = y_min_psi[4];
	for(int i = 4; i < xi_n; i++)
	  y_lim[i+2] = y_max_psi[i];
	
	for(int i = 0; i < xi_n; i++) {
	  y_min[i] = y_min_psi[i];
	  y_max[i] = y_max_psi[i]; 
	}
      }
      // y limits are the same for the 3 ups states in both sqrt(s)
      else if (dataName == "ups1" || dataName == "ups2" || dataName == "ups3") {
	y_lim[0] = y_min_ups[0];
	for(int i = 0; i < 2; i++)
	  y_lim[i+1] = y_max_ups[i];
	y_lim[3] = y_min_ups[2];
	for(int i = 2; i < xi_n; i++)
	  y_lim[i+2] = y_max_ups[i];
	
	for(int i = 0; i < xi_n; i++) {
	  y_min[i] = y_min_ups[i];
	  y_max[i] = y_max_ups[i]; 
	}
      }
      
      // one if condition for each state and sqrt(s)
      if(sqsName == "7") {
	// filling in the j/psi 7 TeV values
	if (dataName == "jpsi") {

	  // xi plots
	  for(int i = 0; i < xi_n; i++) {
	    xi_min[i] = -1;
	    xi_max[i] = i < 4 ? 35 : 6;
	    xi_yax_min[i] = i < 4 ? 1e-5 : 1;
	    xi_yax_max[i] = i < 4 ? 1e2 : 1e4;
	  }
	  
	  xi_jpsi_7(xi_g, xi_h);

	  // y plots
	  for(int i = 0; i < y_n; i++) {
	    y_xi_min[i] = (i+7.) / mass;
	    y_xi_max[i] = (i+8.) / mass;
	  }
	  
	  y_yax_min_lin = 0;
	  y_yax_max_lin = 400;
	  y_yax_min_log = 1;
	  y_yax_max_log = 1e4;

	  y_jpsi_7(y_g, y_n);

	}
	
	// filling in the psi(2S) 7 TeV values
	else if (dataName == "psi2") {
	  
	  // xi plots
	  for(int i = 0; i < xi_n; i++) {
	    xi_min[i] = 0;
	    xi_max[i] = i < 4 ? 30 : 5;
	    xi_yax_min[i] = i < 4 ? 1e-4 : 1e-1;
	    xi_yax_max[i] = i < 4 ? 1e2 : 1e4;
	  }
	  
	  xi_psi2_7(xi_g, xi_h);
	  
	  // y plots
	  for(int i = 0; i < y_n; i++) {
	    y_xi_min[i] = (i+8.) / mass;
	    y_xi_max[i] = (i+9.) / mass;
	  }
	  
	  y_yax_min_lin = 0;
	  y_yax_max_lin = 70;
	  y_yax_min_log = 1e-1;
	  y_yax_max_log = 1e3;
	  
	  y_psi2_7(y_g, y_n);
	  
	}
	
	else if (dataName == "ups1") {
	  
	  // xi plots
	  for(int i = 0; i < xi_n; i++) {
	    xi_min[i] = -1;
	    xi_max[i] = i < 2 ? 14 : 5;
	    xi_yax_min[i] = i < 2 ? 1e-4 : 1e-2;
	    xi_yax_max[i] = i < 2 ? 1e2 : 1e3;
	  }
	  
	  xi_ups1_7(xi_g, xi_h);
	  
	  // y plots
	  for(int i = 0; i < y_n; i++) {
	    y_xi_min[i] = (2.*i+10.)/mass;
	    y_xi_max[i] = (2.*i+12.)/mass;	
	  }
	  
	  y_yax_min_lin = 0;
	  y_yax_max_lin = 40;
	  y_yax_min_log = 1e-1;
	  y_yax_max_log = 1e2;
	  
	  y_ups1_7(y_g, y_n);
	  
	}
	else if (dataName == "ups2") {
	  
	  // xi plots
	  for(int i = 0; i < xi_n; i++) {
	    xi_min[i] = -1;
	    xi_max[i] = i < 2 ? 14 : 5;
	    xi_yax_min[i] = i < 2 ? 1e-4 : 1e-2;
	    xi_yax_max[i] = i < 2 ? 1e2 : 1e3;
	  }
	  
	  xi_ups2_7(xi_g, xi_h);
	  
	  // y plots
	  for(int i = 0; i < y_n; i++) {
	    y_xi_min[i] = (2.*i+10.)/mass;
	    y_xi_max[i] = (2.*i+12.)/mass;	
	  }
	  
	  y_yax_min_lin = 0;
	  y_yax_max_lin = 30;
	  y_yax_min_log = 1e-1;
	  y_yax_max_log = 1e2;
	  
	  y_ups2_7(y_g, y_n);
	  
	}
	else if (dataName == "ups3") {
	  
	  // xi plots
	  for(int i = 0; i < xi_n; i++) {
	    xi_min[i] = -1;
	    xi_max[i] = i < 2 ? 14 : 5;
	    xi_yax_min[i] = i < 2 ? 1e-4 : 1e-2;
	    xi_yax_max[i] = i < 2 ? 1e1 : 1e2;
	  }
	  
	  xi_ups3_7(xi_g, xi_h);
	  
	  // y plots
	  for(int i = 0; i < y_n; i++) {
	    y_xi_min[i] = (2.*i+10.)/mass;
	    y_xi_max[i] = (2.*i+12.)/mass;	
	  }
	  
	  y_yax_min_lin = 0;
	  y_yax_max_lin = 10;
	  y_yax_min_log = 1e-1;
	  y_yax_max_log = 1e2;
	  
	  y_ups3_7(y_g, y_n);
	  
	}
      }
      else if(sqsName == "13") {
	if (dataName == "jpsi") {
	  
	  // xi plots
	  for(int i = 0; i < xi_n; i++) {
	    xi_min[i] = -1;
	    xi_max[i] = i < 4 ? 40 : 6;	
	    xi_yax_min[i] = i < 4 ? 1e-5 : 1;
	    xi_yax_max[i] = i < 4 ? 1e2 : 1e4;
	  }
	  
	  xi_jpsi_13(xi_g, xi_h);
	  
	  // y plots
	  for(int i = 0; i < y_n; i++) {
	    y_xi_min[i] = (i+7.) / mass;
	    y_xi_max[i] = (i+8.) / mass;
	  }
	  
	  y_yax_min_lin = 0;
	  y_yax_max_lin = 600;
	  y_yax_min_log = 1;
	  y_yax_max_log = 1e4;	

	  y_jpsi_13(y_g, y_n);
	  
	}
	else if (dataName == "psi2") {
	  
	  // xi plots
	  for(int i = 0; i < xi_n; i++) {
	    xi_min[i] = 0;
	    xi_max[i] = i < 4 ? 30 : 6;
	    xi_yax_min[i] = i < 4 ? 1e-4 : 1e-1;
	    xi_yax_max[i] = i < 4 ? 1e2 : 1e4;
	  }
	  
	  xi_psi2_13(xi_g, xi_h);
	  
	  // y plots
	  for(int i = 0; i < y_n; i++) {
	    y_xi_min[i] = (i+8.) / mass;
	    y_xi_max[i] = (i+9.) / mass;
	  }
	  
	  y_yax_min_lin = 0;
	  y_yax_max_lin = 100;
	  y_yax_min_log = 5e-1;
	  y_yax_max_log = 5e3;
	  
	  y_psi2_13(y_g, y_n);
	  
	}
      	else if (dataName == "ups1") {
	  
	  // xi plots
	  for(int i = 0; i < xi_n; i++) {
	    xi_min[i] = -1;
	    xi_max[i] = i < 2 ? 14 : 5;
	    xi_yax_min[i] = i < 2 ? 1e-4 : 1e-2;
	    xi_yax_max[i] = i < 2 ? 1e2 : 1e3;
	  }
	  
	  xi_ups1_13(xi_g, xi_h);
	  
	  // y plots
	  for(int i = 0; i < y_n; i++) {
	    y_xi_min[i] = (2.*i+10.)/mass;
	    y_xi_max[i] = (2.*i+12.)/mass;	
	  }
	  
	  y_yax_min_lin = 0;
	  y_yax_max_lin = 60;
	  y_yax_min_log = 1e-1;
	  y_yax_max_log = 1e3;
	  
	  y_ups1_13(y_g, y_n);
	  
	}
	else if (dataName == "ups2") {
	  
	  // xi plots
	  for(int i = 0; i < xi_n; i++) {
	    xi_min[i] = -1;
	    xi_max[i] = i < 2 ? 14 : 5;
	    xi_yax_min[i] = i < 2 ? 1e-4 : 1e-2;
	    xi_yax_max[i] = i < 2 ? 1e2 : 1e3;
	  }
	  
	  xi_ups2_13(xi_g, xi_h);
	  
	  // y plots
	  for(int i = 0; i < y_n; i++) {
	    y_xi_min[i] = (2.*i+10.)/mass;
	    y_xi_max[i] = (2.*i+12.)/mass;	
	  }
	  
	  y_yax_min_lin = 0;
	  y_yax_max_lin = 40;
	  y_yax_min_log = 1e-1;
	  y_yax_max_log = 1e3;
	  
	  y_ups2_13(y_g, y_n);
	  
	}
	else if (dataName == "ups3") {
	  
	  // xi plots
	  for(int i = 0; i < xi_n; i++) {
	    xi_min[i] = -1;
	    xi_max[i] = i < 2 ? 14 : 5;
	    xi_yax_min[i] = i < 2 ? 1e-4 : 1e-2;
	    xi_yax_max[i] = i < 2 ? 1e1 : 1e2;
	  }
	  
	  xi_ups3_13(xi_g, xi_h);
	  
	  // y plots
	  for(int i = 0; i < y_n; i++) {
	    y_xi_min[i] = (2.*i+10.)/mass;
	    y_xi_max[i] = (2.*i+12.)/mass;	
	  }
	  
	  y_yax_min_lin = 0;
	  y_yax_max_lin = 20;
	  y_yax_min_log = 1e-1;
	  y_yax_max_log = 1e3;
	  
	  y_ups3_13(y_g, y_n);
	  
	}
      }
      
      /////////////////////////////////////////
      // part 4 : making the MC histograms

      int y_nbin = xi_n+1;
      TH1F ***y_h = new TH1F**[2];
      for(int ib = 0; ib < 2 ; ib++) {
	y_h[ib] = new TH1F*[y_n];
	for(int j = 0; j < y_n; j++)
	  y_h[ib][j] = new TH1F(Form("y_xi%d_b%d", j, ib), Form("y_xi%d_b%d", j, ib), y_nbin, y_lim);
      }
      
      // variables to open the MC
      Double_t xi;
      Double_t w_gg;
      Double_t w_cos;
      Double_t y;

      Int_t nentries[2];

      for(int ib = 0; ib < 2; ib++)
	{
	  TFile *fin = new TFile(Form("%s_beta%d.root", fName.c_str(), ib+1));
	  TTree *tree = (TTree*)fin->Get(ssName[sqsName].c_str());
      
	  tree->SetBranchAddress("w_gg", &w_gg);
	  tree->SetBranchAddress("w_cos", &w_cos);
	  tree->SetBranchAddress("xi", &xi);
	  tree->SetBranchAddress("y", &y);
      
	  nentries[ib] = (Int_t)tree->GetEntries();
	  cout << nentries[ib] << " events in tree " << endl;
	  int chk = nentries[ib] / 100;

	  for( Int_t i = 0; i < nentries[ib]; i++)
	    {
	      tree->GetEntry(i);
	      
	      // xi part (forward y)
	      for (int k = 0; k < xi_n; k++) 
		if(y < y_max[k] && y > y_min[k] && xi > xi_min_global && xi < 50) {
		  xi_h[ib][k]->Fill(xi, w_gg*w_cos);
		}
	      
	      // y part (also forward y)
	      for (int k = 0; k < y_n; k++)
		if(xi < y_xi_max[k] && xi > y_xi_min[k])
		  y_h[ib][k]->Fill(y, w_gg * w_cos);
	      
	      if((i+1)%chk == 0) {
		cout << (i+1)/chk << "% | " << flush;
	      }
	    }
	  cout << endl;
	  fin->Close();
	}
      
      /////////////////////////////////////////
      // part 5 : plotting
      
      TFile *tf_h = new TFile(Form("MC_vs_Data_%s.root", dataName.c_str()), "UPDATE", "MC histograms and Data graphs");
      
      TCanvas *can = new TCanvas("", "", 700, 700);
      can->SetLogy();
      
      // xi plots
      
      // scaling the histos to nr events + y bin width
      for(int ib = 0; ib < 2 ; ib++) {
	double n = nentries[ib];
	for(int j = 0; j < xi_n; j++) {
	  xi_h[ib][j]->Scale(normFactor_b[ib]/(n*(y_max[j]-y_min[j])));
	  // scaling each xi bin by its width
	  int nb = xi_h[ib][j]->GetNbinsX();
	  for(int bin = 0; bin < nb; bin++) {
	    double wd = xi_h[ib][j]->GetBinWidth(bin+1) * mass;
	    double vl = xi_h[ib][j]->GetBinContent(bin+1);
	    xi_h[ib][j]->SetBinContent(bin+1, vl/wd);
	    if(xi_h[ib][j]->GetBinLowEdge(bin+1) < xi_min_global)
	      xi_h[ib][j]->SetBinContent(bin+1, 0);
	  }
	}
      }

      // cycle for all the y bins (histo + graph)
      for(int j_y = 0; j_y < xi_n; j_y++) {
	// plot frame + title etc
	TH1F *fc = can->DrawFrame(xi_min[j_y], xi_yax_min[j_y], xi_max[j_y], xi_yax_max[j_y]);
	fc->SetXTitle("#xi");
	fc->SetYTitle("d#sigma / d#xidy (nb/GeV)");
	fc->GetYaxis()->SetTitleOffset(1.3);
	fc->SetTitle(Form("%s TeV %.1f < y < %.1f", sqsName.c_str(), y_min[j_y], y_max[j_y]));
	can->Modified();
	can->SetTitle("");

	// draw legend
	TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
	leg->SetTextSize(0.03);
	
	// draw histogram
	for(int ib = 0; ib < 2; ib++) {
	  xi_h[ib][j_y]->SetStats(0);
	  xi_h[ib][j_y]->SetLineColor(ib+1);
	  xi_h[ib][j_y]->Draw("hist same");
		
	  //save in ROOT fle
	  xi_h[ib][j_y]->SetName(Form("%s TeV %.1f < y < %.1f beta = %d", sqsName.c_str(), y_min[j_y], y_max[j_y], ib+1));
	  xi_h[ib][j_y]->Write();
	  
	  leg->AddEntry(xi_h[ib][j_y], Form("MC #beta=%d", ib+1), "l");
	}
	
	// draw data graph
	xi_g[j_y]->SetMarkerStyle(20);
	xi_g[j_y]->SetMarkerSize(.75);
	xi_g[j_y]->Draw("p");
	// save in ROOT file
	xi_g[j_y]->SetName(Form("%s TeV %.1f < y < %.1f", sqsName.c_str(), y_min[j_y], y_max[j_y]));
	xi_g[j_y]->Write();
	
	leg->AddEntry(xi_g[j_y], "data", "pl");
	
	// save plot
	leg->Draw();
	can->SaveAs(Form("plots_%s/xi_%s_y%d.pdf", dataName.c_str(), sqsName.c_str(), j_y+1));
	can->Clear();
      }
      
      // y plots
      
      // scaling the histos to nr events
      for(int ib = 0; ib < 2; ib++) {
	double n = nentries[ib];
	for(int j = 0; j < y_n; j++) { 
	  y_h[ib][j]->Scale(normFactor_b[ib]/(n*(y_xi_max[j]-y_xi_min[j])*mass));
	  // scaling each y bin to its width
	  int nb = y_h[ib][j]->GetNbinsX();
	  for(int bin = 0; bin < nb; bin++) {
	    double wd = y_h[ib][j]->GetBinWidth(bin+1);
	    double vl = y_h[ib][j]->GetBinContent(bin+1);
	    y_h[ib][j]->SetBinContent(bin+1, vl/wd);
	  }
	}
      }
      
      // draw frame for the y plot (log-plot)
      TH1F *fc_log = can->DrawFrame(0, y_yax_min_log, 5, y_yax_max_log);
      fc_log->SetXTitle("y");
      fc_log->SetYTitle("d#sigma / d#xidy (nb/GeV)");
      fc_log->GetYaxis()->SetTitleOffset(1.3);
      fc_log->SetTitle(Form("%s TeV", sqsName.c_str()));
      can->Modified();
      can->SetTitle("");
      
      // plot the single y plot        
      TLegend *leg = new TLegend(0.65, 0.6, 0.9, 0.9);
      leg->SetTextSize(0.03);
      
      // plot the histograms
      for(int ib = 0; ib < 2 ; ib++) {
	for(int j_xi = 0; j_xi < y_n; j_xi++) {
	  y_h[ib][j_xi]->SetLineColor(j_xi+1);
	  if(j_xi+1 == 5)  y_h[ib][j_xi]->SetLineColor(kYellow+1);
	  y_h[ib][j_xi]->Draw("hist same");
	  leg->AddEntry(y_h[ib][j_xi], Form("%.0f < p_{T} < %.0f beta%d", y_xi_min[j_xi]*mass, y_xi_max[j_xi]*mass, ib+1), "l");
	  y_h[ib][j_xi]->SetName(Form("%s TeV %.0f < p_{T} < %.0f beta = %d", sqsName.c_str(), y_xi_min[j_xi]*mass, y_xi_max[j_xi]*mass, ib+1));
	y_h[ib][j_xi]->Write();
	}
      }

      // plot the graphs
      for(int j_xi = 0; j_xi < y_n; j_xi++) {
	y_g[j_xi]->SetMarkerStyle(20);
	y_g[j_xi]->SetMarkerSize(.75);
	y_g[j_xi]->SetMarkerColor(j_xi+1);
	y_g[j_xi]->SetLineColor(j_xi+1);
	if(j_xi+1 == 5)  {
	  y_g[j_xi]->SetMarkerColor(kYellow+1);
	  y_g[j_xi]->SetLineColor(kYellow+1);
	}
	y_g[j_xi]->Draw("p");
	y_g[j_xi]->SetName(Form("%s TeV %.0f < p_{T} < %.0f", sqsName.c_str(), y_xi_min[j_xi]*mass, y_xi_max[j_xi]*mass));
	y_g[j_xi]->Write();
      }
      
      leg->Draw();
      
      can->SaveAs(Form("plots_%s/y_dist_log_%s.pdf", dataName.c_str(), sqsName.c_str()));
      
      can->Clear();
      
      // same for linear scale
      can->SetLogy(0);

      // frame with lin limits
      TH1F *fc_lin = can->DrawFrame(0, y_yax_min_lin, 5, y_yax_max_lin);
      fc_lin->SetXTitle("y");
      fc_lin->SetYTitle("d#sigma / d#xidy (nb/GeV)");
      fc_lin->GetYaxis()->SetTitleOffset(1.3);
      fc_lin->SetTitle(Form("%s TeV", sqsName.c_str()));
      can->Modified();
      can->SetTitle("");

      // histograms
      for(int ib = 0; ib < 2; ib++) {
	for(int j_xi = 0; j_xi < y_n; j_xi++) {
	  y_h[ib][j_xi]->SetLineColor(j_xi+1);
	  if(j_xi+1 == 5)  y_h[ib][j_xi]->SetLineColor(kYellow+1);
	  y_h[ib][j_xi]->Draw("hist same");
	}
      }

      // data graphs
      for(int j_xi = 0; j_xi < y_n; j_xi++) {
	y_g[j_xi]->SetMarkerStyle(20);
	y_g[j_xi]->SetMarkerSize(.75);
	y_g[j_xi]->SetMarkerColor(j_xi+1);
	y_g[j_xi]->SetLineColor(j_xi+1);
	if(j_xi+1 == 5)  {
	  y_g[j_xi]->SetMarkerColor(kYellow+1);
	  y_g[j_xi]->SetLineColor(kYellow+1);
	}

	y_g[j_xi]->Draw("p");
      }

      // no need to redo the legend, same entries
      leg->Draw();
      can->SaveAs(Form("plots_%s/y_dist_%s.pdf", dataName.c_str(), sqsName.c_str()));
      can->Clear();
      can->Destructor();

      tf_h->Write();
      tf_h->Close();

    }
  
}
