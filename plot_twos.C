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
  map<string, string> ssName;
  ssName["7"] = "qqbar_0";
  ssName["13"] = "qqbar_1";
  
  double y_min_1C[6] = {0, 2, 2.5, 3, 3.5, 4};
  double y_max_1C[6] = {1.2, 2.5, 3, 3.5, 4, 4.5};
  double y_min_4C[9] = {0, 0.3, 0.6, 0.9, 2, 2.5, 3, 3.5, 4};
  double y_max_4C[9] = {0.3, 0.6, 0.9, 1.2, 2.5, 3, 3.5, 4, 4.5};

  /////////////////////////////////////////
  // part 1 : variables to be changed on each run
  
  // choose state we're plotting
  string dataName = "psi2";   // "jpsi", "psi2", "ups1"
  // choose sqrt(s) we're plotting
  string sqsNames[2] = {"7", "13"};   // "7", "13"

  string fName = "MC_res.root";   // file to open
  double xi_bin_width = 0.1;   // constant for easier normalization

  double y_plot_min = 0, y_plot_max = 5;
  
  /////////////////////////////////////////
  // part 2 : vars defined by above or constant
  double xi_normFactor = 1.5e4;   // overall normalization btw data, MC
  double y_normFactor = 1.5e4;

  // extra part: create ROOT file and save results

  TFile *tf = new TFile("MC_vs_Data.root", "RECREATE", "MC histograms and Data graphs");
  TTree *norm_factors = new TTree("norms", "norm factors");
  
  TBranch *xi_nf = norm_factors->Branch("xi_normFactor", &xi_normFactor);
  TBranch *y_nf = norm_factors->Branch("y_normFactor", &y_normFactor);

  norm_factors->Fill();
  //norm_factors->Write();
  tf->Write();
  tf->Close();
  
  for(int i_sqs = 0; i_sqs < 2; i_sqs++)
    {
      string sqsName = sqsNames[i_sqs];
  
      const int xi_n = dataName == "psi2" && sqsName == "7" ? 9 : 6;
      const int y_n = sqsName == "7" ? dataName == "jpsi" ? 7 : dataName == "psi2" ? 6 : dataName == "ups1" ? 3 : 0 : sqsName == "13" ? dataName == "psi2" ? 6 : 0 : 0;
      TGraphAsymmErrors **xi_g = new TGraphAsymmErrors*[xi_n];
      TGraphAsymmErrors **y_g = new TGraphAsymmErrors*[y_n];

      double y_min[xi_n]; 
      double y_max[xi_n];
      double y_lim[xi_n+2];
      double xi_min[xi_n], xi_max[xi_n];
      double xi_yax_min[xi_n], xi_yax_max[xi_n];
      
      double y_xi_min[y_n], y_xi_max[y_n];
      double y_yax_min_lin, y_yax_max_lin;
      double y_yax_min_log, y_yax_max_log;
      
      double mass = stateMass[dataName];

      TH1F **xi_h = new TH1F*[xi_n]; // histograms

      /////////////////////////////////////////
      // part 3 : making the data TGraphs
  
      if (dataName == "jpsi" && sqsName == "7")
	{
	  y_lim[0] = y_min_1C[0];
	  y_lim[1] = y_max_1C[0];
	  for(int i = 0; i < xi_n-1; i++)
	      y_lim[i+2] = y_min_1C[i+1];
	  y_lim[xi_n+1] = y_max_1C[xi_n-1];

	  // xi plots
	  for(int i = 0; i < xi_n; i++) {
	    y_min[i] = y_min_1C[i];
	    y_max[i] = y_max_1C[i];
	    xi_min[i] = -1;
	    xi_max[i] = i == 0 ? 40 : 11;
	    xi_yax_min[i] = 1e-5;
	    xi_yax_max[i] = 1e4;
	  }
	  
	  xi_jpsi_7(xi_g, xi_h, 1.);
      
	  // y plots
	  for(int i = 0; i < y_n; i++) {
	    y_xi_min[i] = (i+7.) / mass;
	    y_xi_max[i] = (i+8.) / mass;
	  }
	  
	  y_yax_min_lin = 0;
	  y_yax_max_lin = 0.4;
	  y_yax_min_log = 1e-3;
	  y_yax_max_log = 2;
	  
	  y_jpsi_7(y_g, 1., y_n);
	  
	}
      if (dataName == "psi2" && sqsName == "7")
	{
	  y_lim[0] = y_min_4C[0];
	  for(int i = 0; i < 4; i++)
	      y_lim[i+1] = y_max_4C[i];
	  y_lim[5] = y_min_4C[4];
	  for(int i = 4; i < xi_n; i++)
	    y_lim[i+2] = y_max_4C[i];
	  
	  // xi plots
	  for(int i = 0; i < xi_n; i++) {
	    y_min[i] = y_min_4C[i];
	    y_max[i] = y_max_4C[i];
	    xi_min[i] = 0;
	    xi_max[i] = i < 4 ? 20 : 6;
	    xi_yax_min[i] = i < 4 ? 1e-3 : 1e-1;
	    xi_yax_max[i] = i < 4 ? 1e2 : 1e4;
	  }
	  
	  xi_psi2_7(xi_g, xi_h, 1.);
	  
	  // y plots
	  for(int i = 0; i < y_n; i++) {
	    y_xi_min[i] = (i+8.) / mass;
	    y_xi_max[i] = (i+9.) / mass;
	  }
	  
	  y_yax_min_lin = 0;
	  y_yax_max_lin = 70;
	  y_yax_min_log = 1e-1;
	  y_yax_max_log = 1e3;
	  
	  y_psi2_7(y_g, 1., y_n);
	  
	}
      if (dataName == "psi2" && sqsName == "13")
	{
	  y_lim[0] = y_min_1C[0];
	  y_lim[1] = y_max_1C[0];
	  for(int i = 0; i < xi_n-1; i++)
	      y_lim[i+2] = y_min_1C[i+1];
	  y_lim[xi_n+1] = y_max_1C[xi_n-1];
	  
	  // xi plots
	  for(int i = 0; i < xi_n; i++) {
	    y_min[i] = y_min_1C[i];
	    y_max[i] = y_max_1C[i];
	    xi_min[i] = 0;
	    xi_max[i] = 6;
	    xi_yax_min[i] = 1e-1;
	    xi_yax_max[i] = 1e4;
	  }
	  
	  xi_psi2_13(xi_g, xi_h, 1.);
	  
	  // y plots
	  for(int i = 0; i < y_n; i++) {
	    y_xi_min[i] = (i+8.) / mass;
	    y_xi_max[i] = (i+9.) / mass;
	  }
	  
	  y_yax_min_lin = 0;
	  y_yax_max_lin = 100;
	  y_yax_min_log = 5e-1;
	  y_yax_max_log = 5e3;
	  
	  y_psi2_13(y_g, 1., y_n);
	  
	}
      if (dataName == "ups1" && sqsName == "7")
	{
	  y_lim[0] = y_min_1C[0];
	  y_lim[1] = y_max_1C[0];
	  for(int i = 0; i < xi_n-1; i++)
	    y_lim[i+2] = y_min_1C[i+1];
	  y_lim[xi_n+1] = y_max_1C[xi_n-1];
	  
	  // xi plots
	  for(int i = 0; i < xi_n; i++) {
	    y_min[i] = y_min_1C[i];
	    y_max[i] = y_max_1C[i];
	    xi_min[i] = 0;
	    xi_max[i] = 6;
	    xi_yax_min[i] = 1e-3;
	    xi_yax_max[i] = 1e2;
	  }
	  
	  xi_ups1_7(xi_g, xi_h, 1.);
	  
	  // y plots
	  for(int i = 0; i < y_n; i++) {
	    y_xi_min[i] = (i+10.)/mass;
	    y_xi_max[i] = (i+11.)/mass;	
	  }
	  
	  y_yax_min_lin = 0;
	  y_yax_max_lin = 0.4;
	  y_yax_min_log = 1e-3;
	  y_yax_max_log = 2;
	  
	  y_ups1_7(y_g, 1., y_n);
	  
	}
      
      /////////////////////////////////////////
      // part 4 : making the MC histograms

      int y_nbin = xi_n+1;
      TH1F **y_h = new TH1F*[y_n];
      for(int j = 0; j < y_n; j++)
	y_h[j] = new TH1F(Form("y_xi%d", j), Form("y_xi%d", j), y_nbin, y_lim); 
      
      // variables to open the MC
      Double_t xi;
      Double_t w_gg;
      Double_t w_cos;
      Double_t y;
      
      TFile *fin = new TFile(Form("%s", fName.c_str()));
      TTree *tree = (TTree*)fin->Get(ssName[sqsName].c_str());
      
      tree->SetBranchAddress("w_gg", &w_gg);
      tree->SetBranchAddress("w_cos", &w_cos);
      tree->SetBranchAddress("xi", &xi);
      tree->SetBranchAddress("y", &y);
      
      Int_t nentries = (Int_t)tree->GetEntries();
      cout << nentries << " events in tree " << endl;
      int chk = nentries / 100;

      for( Int_t i = 0; i < nentries; i++)
	{
	  tree->GetEntry(i);

	  // xi part (only forward (pos) y)
	  for (int k = 0; k < xi_n; k++) 
	    if(y < y_max[k] && y > y_min[k] && xi > 2 && xi < 50) {
	      xi_h[k]->Fill(xi, w_gg*w_cos);
	    }
	  
	  // y part (always pos value)
	  for (int k = 0; k < y_n; k++)
	    if(xi < y_xi_max[k] && xi > y_xi_min[k])
	      y_h[k]->Fill(y, w_gg * w_cos);
	  
	  if((i+1)%chk == 0) {
	    cout << (i+1)/chk << "% | " << flush;
	  }
	}
      cout << endl;
      fin->Close();
      
      /////////////////////////////////////////
      // part 5 : plotting

      TFile *tf_h = new TFile("MC_vs_Data.root", "UPDATE", "MC histograms and Data graphs");
      
      TCanvas *can = new TCanvas("", "", 700, 700);
      can->SetLogy();
      
      // xi plots
      
      // scaling the histos to the first integral
      double n = nentries;//xi_h[0]->Integral("width");
      for(int j = 0; j < xi_n; j++) {
	xi_h[j]->Scale(xi_normFactor/(n*(y_max[j]-y_min[j])));
	int nb = xi_h[j]->GetNbinsX();
	for(int bin = 0; bin < nb; bin++) {
	  double wd = xi_h[j]->GetBinWidth(bin+1) * mass;
	  double vl = xi_h[j]->GetBinContent(bin+1);
	  xi_h[j]->SetBinContent(bin+1, vl/wd);
	}
      }
      
      // cycle for all the y bins (histo + graph)
      for(int j_y = 0; j_y < xi_n; j_y++) {
	TH1F *fc = can->DrawFrame(xi_min[j_y], xi_yax_min[j_y], xi_max[j_y], xi_yax_max[j_y]);
	fc->SetXTitle("#xi");
	fc->SetYTitle("d#sigma / d#xidy (nb/GeV)");
	fc->GetYaxis()->SetTitleOffset(1.3);
	fc->SetTitle(Form("%s %.1f < %s < %.1f", j_y < 1 ? "CMS" : "LHCb", y_min[j_y], "y", y_max[j_y]));
	can->Modified();
	can->SetTitle("");
	
	xi_h[j_y]->SetStats(0);
	xi_h[j_y]->SetLineColor(1);
	//	xi_h[j_y]->GetXaxis()->SetRangeUser(xi_min[j_y], xi_max[j_y]);
	//xi_h[j_y]->GetYaxis()->SetRangeUser(xi_yax_min[j_y], xi_yax_max[j_y]);
	xi_h[j_y]->Draw("hist same");
	
	TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
	leg->SetTextSize(0.03);
	leg->AddEntry(xi_h[j_y], Form("MC"), "l");
	
	xi_g[j_y]->SetMarkerStyle(20);
	xi_g[j_y]->SetMarkerSize(.75);
	xi_g[j_y]->Draw("p");
	leg->AddEntry(xi_g[j_y], "data", "pl");

	if(xi_h[j_y]->GetNbinsX() > 1) {
	  xi_h[j_y]->SetName(Form("%s TeV %.1f < y < %.1f", sqsName.c_str(), y_min[j_y], y_max[j_y]));
	  xi_h[j_y]->Write();
	  xi_g[j_y]->SetName(Form("%s TeV %.1f < y < %.1f", sqsName.c_str(), y_min[j_y], y_max[j_y]));
	  xi_g[j_y]->Write();
	}
	
	leg->Draw();
	can->SaveAs(Form("plots/xi_%s_%s.pdf", sqsName == "7" ? j_y < 4? Form("CMS_y%d", j_y+1) : Form("LHCb_y%d", j_y-3) : j_y < 1 ? "CMS" : Form("LHCb_y%d", j_y), sqsName.c_str()));
	can->Clear();
      }
      
      // y plots
      
      // scaling the histos to the CMS integral
      n = nentries;//y_h[0]->Integral("width");
      for(int j = 0; j < y_n; j++) { 
	y_h[j]->Scale(y_normFactor/n);
	int nb = y_h[j]->GetNbinsX();
	for(int bin = 0; bin < nb; bin++) {
	  double wd = y_h[j]->GetBinWidth(bin+1);
	  double vl = y_h[j]->GetBinContent(bin+1);
	  y_h[j]->SetBinContent(bin+1, vl/wd);
	}

      }

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
      
      for(int j_xi = 0; j_xi < y_n; j_xi++) {
	y_h[j_xi]->SetLineColor(j_xi+1);
	if(j_xi+1 == 5)  y_h[j_xi]->SetLineColor(kYellow+1);
	y_h[j_xi]->Draw("hist same");
	leg->AddEntry(y_h[j_xi], Form("%.0f < p_{T} < %.0f", y_xi_min[j_xi]*mass, y_xi_max[j_xi]*mass), "l");
	y_h[j_xi]->SetName(Form("%s TeV %.0f < p_{T} < %.0f", sqsName.c_str(), y_xi_min[j_xi]*mass, y_xi_max[j_xi]*mass));
	y_h[j_xi]->Write();
      }
      
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
      
      can->SaveAs(Form("plots/y_dist_log_%s.pdf", sqsName.c_str()));
      
      can->Clear();
      
      // same for linear scale
      can->SetLogy(0);

      TH1F *fc_lin = can->DrawFrame(0, y_yax_min_lin, 5, y_yax_max_lin);
      fc_lin->SetXTitle("y");
      fc_lin->SetYTitle("d#sigma / d#xidy (nb/GeV)");
      fc_lin->GetYaxis()->SetTitleOffset(1.3);
      fc_lin->SetTitle(Form("%s TeV", sqsName.c_str()));
      can->Modified();
      can->SetTitle("");

      
      // double cycle for all 4 plots of n_xi = 8 histos each  
      y_h[0]->SetStats(0);
      y_h[0]->SetLineColor(1);
      y_h[0]->Draw("hist same");
      
      /*TLegend *leg = new TLegend(0.65, 0.6, 0.9, 0.9);
	leg->SetTextSize(0.03);
	leg->AddEntry(y_h[0], Form("%.0f < p_{T} < %.0f", y_xi_min[0]*mass, y_xi_max[0]*mass), "l");*/
      
      for(int j_xi = 1; j_xi < y_n; j_xi++) {
	y_h[j_xi]->SetLineColor(j_xi+1);
	if(j_xi+1 == 5)  y_h[j_xi]->SetLineColor(kYellow+1);
	y_h[j_xi]->Draw("hist same");
	//leg->AddEntry(y_h[j_xi], Form("%.0f < p_{T} < %.0f", y_xi_min[j_xi]*mass, y_xi_max[j_xi]*mass), "l");
      }
      
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
      
      leg->Draw();
      can->SaveAs(Form("plots/y_dist_%s.pdf", sqsName.c_str()));
      can->Clear();

      // FINAL PART: PLOT THE PULLS

      TGraphAsymmErrors **xi_p = new TGraphAsymmErrors*[xi_n];
      TGraphAsymmErrors **y_p = new TGraphAsymmErrors*[y_n];
      xi_Pulls(xi_g, xi_h, xi_p, xi_n);
      y_Pulls(y_g, y_h, y_p, y_n);

      // xi pulls
      
      TH1F *fc_pull = can->DrawFrame(0, -6, 6, 6);
      fc_pull->SetXTitle("#xi");
      fc_pull->SetYTitle("pulls");
      fc_pull->GetYaxis()->SetTitleOffset(1.3);
      fc_pull->SetTitle(Form("#xi pulls %s TeV", sqsName.c_str()));
      can->Modified();
      can->SetTitle("");

      TLegend *leg_xip = new TLegend(0.1, 0.6, 0.35, 0.9);
      leg_xip->SetTextSize(0.03);

      for(int i = 0; i < xi_n; i++) {
	if(sqsName == "13" && i == 0) continue;
	xi_p[i]->SetMarkerColor(i+1);
	xi_p[i]->SetMarkerStyle(20);
	xi_p[i]->SetMarkerSize(.75);
	xi_p[i]->SetLineColor(i+1);
	if(i+1 == 5)
	  {
	    xi_p[i]->SetLineColor(kYellow+1);
	    xi_p[i]->SetMarkerColor(kYellow+1);
	  }
	xi_p[i]->Draw("pl same");
	leg_xip->AddEntry(xi_p[i], Form("%.1f < y < %.1f", y_min[i], y_max[i]), "l");

      }

      leg_xip->Draw();
 
      TF1* line = new TF1("line", "0", 0, 6);
      line->SetLineStyle(kDashed);
      line->SetLineColor(kBlack);
      line->Draw("same");

      TLine* vert = new TLine(2, -6, 2, 6);
      vert->SetLineStyle(kDashed);
      vert->SetLineColor(kBlack);
      vert->Draw("same");
      
      can->SaveAs(Form("plots/xi_pulls_%s.pdf", sqsName.c_str()));
      can->Clear();

      // y pulls
      
      TH1F *fc_y_pull = can->DrawFrame(0, -6, 5, 6);
      fc_y_pull->SetXTitle("y");
      fc_y_pull->SetYTitle("pulls");
      fc_y_pull->GetYaxis()->SetTitleOffset(1.3);
      fc_y_pull->SetTitle(Form("y pulls %s TeV", sqsName.c_str()));
      can->Modified();
      can->SetTitle("");

      TLegend *leg_yp = new TLegend(0.1, 0.6, 0.35, 0.9);
      leg_yp->SetTextSize(0.03);
      
      for(int i = 0; i < y_n; i++) {
	y_p[i]->SetMarkerColor(i+1);
	y_p[i]->SetLineColor(i+1);	
	if(i+1 == 5)  {
	  y_p[i]->SetMarkerColor(kYellow+1);
	  y_p[i]->SetLineColor(kYellow+1);
	}
	y_p[i]->SetMarkerStyle(20);
	y_p[i]->SetMarkerSize(.75);
	y_p[i]->Draw("pl same");
	leg_yp->AddEntry(y_p[i], Form("%.0f < p_{T} < %.0f", y_xi_min[i]*mass, y_xi_max[i]*mass), "l");
		
      }

      leg_yp->Draw();

      line->SetLineStyle(kDashed);
      line->SetLineColor(kBlack);
      line->Draw("same");
      
      can->SaveAs(Form("plots/y_pulls_%s.pdf", sqsName.c_str()));
      
      can->Destructor();
      tf_h->Write();
      tf_h->Close();

    }

}

