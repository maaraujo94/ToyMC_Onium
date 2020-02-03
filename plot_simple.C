/* code to plot toy MC histos vs data graphs
   - make TGraphs (external code)
   - make THistos
   - plot both together
   process is repeated for both xi plots in y bins and y plots in xi bins
*/

#import "plot_codes.C"

void plot_simple()
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

  /////////////////////////////////////////
  // part 1 : variables to be changed on each run
  
  // choose state we're plotting
  string dataName = "psi2";   // "jpsi", "psi2", "ups1"
  // choose sqrt(s) we're plotting
  string sqsName = "7";   // "7", "13"

  string fName = "MC_res.root";   // file to open
  double xi_normFactor = 2.5;   // overall normalization btw data, MC
  double xi_bin_width = 0.1;   // constant for easier normalization

  double y_plot_min = 0, y_plot_max = 5;
  double y_normFactor = 1.9;
  
  /////////////////////////////////////////
  // part 2 : vars defined by above or constant
  const int xi_n = 6;
  const int y_n = sqsName == "7" ? dataName == "jpsi" ? 7 : dataName == "psi2" ? 6 : dataName == "ups1" ? 3 : 0 : 0;
  TGraphAsymmErrors **xi_g = new TGraphAsymmErrors*[xi_n];
  TGraphAsymmErrors **y_g = new TGraphAsymmErrors*[y_n];

  double y_min[xi_n] = {0, 2, 2.5, 3, 3.5, 4};
  double y_max[xi_n] = {1.2, 2.5, 3, 3.5, 4, 4.5};
  double xi_min[xi_n], xi_max[xi_n];
  double xi_yax_min[xi_n], xi_yax_max[xi_n];

  double y_xi_min[y_n], y_xi_max[y_n];
  double y_yax_min_lin, y_yax_max_lin;
  double y_yax_min_log, y_yax_max_log;

  double mass = stateMass[dataName];

  /////////////////////////////////////////
  // part 3 : making the data TGraphs
  
  if (dataName == "jpsi" && sqsName == "7")
    {
      // xi plots
      for(int i = 0; i < xi_n; i++) {
	xi_min[i] = -1;
	xi_max[i] = i == 0 ? 40 : 11;
	xi_yax_min[i] = 1e-5;
	xi_yax_max[i] = 1e4;
      }
      
      xi_jpsi_7(xi_g, xi_normFactor);

      // y plots
      for(int i = 0; i < y_n; i++) {
	y_xi_min[i] = (i+7.) / mass;
	y_xi_max[i] = (i+8.) / mass;
      }
      
      y_yax_min_lin = 0;
      y_yax_max_lin = 0.4;
      y_yax_min_log = 1e-3;
      y_yax_max_log = 2;

      y_jpsi_7(y_g, y_normFactor, y_n);

    }
  if (dataName == "psi2" && sqsName == "7")
    {
      // xi plots
      for(int i = 0; i < xi_n; i++) {
	xi_min[i] = 0;
	xi_max[i] = 6;
	xi_yax_min[i] = 1e-3;
	xi_yax_max[i] = 1e2;
      }
      
      xi_psi2_7(xi_g, xi_normFactor);

      // y plots
      for(int i = 0; i < y_n; i++) {
	y_xi_min[i] = (i+8.) / mass;
	y_xi_max[i] = (i+9.) / mass;
      }

      y_yax_min_lin = 0;
      y_yax_max_lin = 0.4;
      y_yax_min_log = 1e-3;
      y_yax_max_log = 2;

      y_psi2_7(y_g, y_normFactor, y_n);

    }
  if (dataName == "ups1" && sqsName == "7")
    {
      // xi plots
      for(int i = 0; i < xi_n; i++) {
	xi_min[i] = 0;
	xi_max[i] = 6;
	xi_yax_min[i] = 1e-3;
	xi_yax_max[i] = 1e2;
      }
      
      xi_ups1_7(xi_g, xi_normFactor);

      // y plots
      for(int i = 0; i < y_n; i++) {
	y_xi_min[i] = (i+10.)/mass;
	y_xi_max[i] = (i+11.)/mass;	
      }

      y_yax_min_lin = 0;
      y_yax_max_lin = 0.4;
      y_yax_min_log = 1e-3;
      y_yax_max_log = 2;

      y_ups1_7(y_g, y_normFactor, y_n);

    }

  /////////////////////////////////////////
  // part 4 : making the MC histograms

  int xi_n_bins[xi_n]; // bin number
  for(int i = 0; i < xi_n; i++) {
    xi_n_bins[i] = (int)((xi_max[i] - xi_min[i]) / xi_bin_width);
  }
  
  TH1F **xi_h = new TH1F*[xi_n]; // histograms
  for(int j = 0; j < xi_n; j++)
    xi_h[j] = new TH1F(Form("xi_y%d", j), Form("%s %.1f < |y| < %.1f", j < 1 ? "CMS" : "LHCb", y_min[j], y_max[j]), xi_n_bins[j], xi_min[j], xi_max[j]); 
  TH1F **y_h = new TH1F*[y_n];
  for(int j = 0; j < y_n; j++)
    y_h[j] = new TH1F(Form("y_xi%d", j), Form("y_xi%d", j), 50, y_plot_min, y_plot_max); 

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

      for (int k = 0; k < xi_n; k++)
	if(abs(y) < y_max[k] && abs(y) > y_min[k] && xi > 2 && xi < 50)
	  xi_h[k]->Fill(xi, w_gg * w_cos);

      for (int k = 0; k < y_n; k++)
	if(xi < y_xi_max[k] && xi > y_xi_min[k])
	  y_h[k]->Fill(abs(y), w_gg * w_cos);
      
      if((i+1)%chk == 0) {
	cout << (i+1)/chk << "% | " << flush;
      }
    }
  cout << endl;

  /////////////////////////////////////////
  // part 5 : plotting

  TCanvas *can = new TCanvas("", "", 700, 700);
  can->SetLogy();

  // xi plots
  
  // scaling the histos to the first integral
  double n = xi_h[0]->Integral("width");
  for(int j = 0; j < xi_n; j++) {
    xi_h[j]->Scale(1./(2.*n*(y_max[j]-y_min[j])));
  }
  
  // cycle for all the y bins (histo + graph)
  for(int j_y = 0; j_y < xi_n; j_y++) {    
    xi_h[j_y]->SetStats(0);
    xi_h[j_y]->SetLineColor(1);
    xi_h[j_y]->GetXaxis()->SetTitle("#xi");
    xi_h[j_y]->GetYaxis()->SetTitle(Form("Events/%.1f", xi_bin_width));
    xi_h[j_y]->GetYaxis()->SetTitleOffset(1.3);
    xi_h[j_y]->GetYaxis()->SetRangeUser(xi_yax_min[j_y], xi_yax_max[j_y]);
    xi_h[j_y]->Draw("hist");
    
    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->SetTextSize(0.03);
    leg->AddEntry(xi_h[j_y], Form("MC"), "l");
    
    xi_g[j_y]->SetMarkerStyle(20);
    xi_g[j_y]->SetMarkerSize(.75);
    xi_g[j_y]->Draw("p");
    leg->AddEntry(xi_g[j_y], "data", "pl");
    
    leg->Draw();
    can->SaveAs(Form("plots/xi_%s.pdf", j_y < 1? "CMS" : Form("LHCb_y%d", j_y)));
    can->Clear();
  }

  // y plots
  
  // scaling the histos to the CMS integral
  n = y_h[0]->Integral("width");
  for(int j = 0; j < y_n; j++) 
    y_h[j]->Scale(1./n);
  
  // plot the single y plot  
  y_h[0]->SetStats(0);
  y_h[0]->SetLineColor(1);
  y_h[0]->GetXaxis()->SetTitle("y");
  y_h[0]->GetYaxis()->SetTitle("Events");
  y_h[0]->SetTitle("");
  y_h[0]->GetYaxis()->SetTitleOffset(1.3);
  y_h[0]->GetYaxis()->SetRangeUser(y_yax_min_log, y_yax_max_log);
  y_h[0]->Draw("hist"); 
  
  TLegend *leg = new TLegend(0.65, 0.6, 0.9, 0.9);
  leg->SetTextSize(0.03);
  leg->AddEntry(y_h[0], Form("%.0f < p_{T} < %.0f", y_xi_min[0]*mass, y_xi_max[0]*mass), "l");

  for(int j_xi = 1; j_xi < y_n; j_xi++) {
    y_h[j_xi]->SetLineColor(j_xi+1);
    y_h[j_xi]->Draw("hist same");
    leg->AddEntry(y_h[j_xi], Form("%.0f < p_{T} < %.0f", y_xi_min[j_xi]*mass, y_xi_max[j_xi]*mass), "l");
    
  }
  
  for(int j_xi = 0; j_xi < y_n; j_xi++) {
    y_g[j_xi]->SetMarkerStyle(20);
    y_g[j_xi]->SetMarkerSize(.75);
    y_g[j_xi]->SetMarkerColor(j_xi+1);
    y_g[j_xi]->SetLineColor(j_xi+1);
    y_g[j_xi]->Draw("p");
  }
  
  leg->Draw();
  
  can->SaveAs("plots/y_dist_log.pdf");

  can->Clear();

  // same for linear scale
  can->SetLogy(0);

  // double cycle for all 4 plots of n_xi = 8 histos each  
  y_h[0]->SetStats(0);
  y_h[0]->SetLineColor(1);
  y_h[0]->GetXaxis()->SetTitle("y");
  y_h[0]->GetYaxis()->SetTitle("Events");
  y_h[0]->SetTitle("");
  y_h[0]->GetYaxis()->SetTitleOffset(1.3);
  y_h[0]->GetYaxis()->SetRangeUser(y_yax_min_lin, y_yax_max_lin);
  y_h[0]->Draw("hist");

  /*TLegend *leg = new TLegend(0.65, 0.6, 0.9, 0.9);
  leg->SetTextSize(0.03);
  leg->AddEntry(y_h[0], Form("%.0f < p_{T} < %.0f", y_xi_min[0]*mass, y_xi_max[0]*mass), "l");*/
    
  for(int j_xi = 1; j_xi < y_n; j_xi++) {
    y_h[j_xi]->SetLineColor(j_xi+1);
    y_h[j_xi]->Draw("hist same");
    //leg->AddEntry(y_h[j_xi], Form("%.0f < p_{T} < %.0f", y_xi_min[j_xi]*mass, y_xi_max[j_xi]*mass), "l");
  }
  
  for(int j_xi = 0; j_xi < y_n; j_xi++) {
    y_g[j_xi]->SetMarkerStyle(20);
    y_g[j_xi]->SetMarkerSize(.75);
    y_g[j_xi]->SetMarkerColor(j_xi+1);
    y_g[j_xi]->SetLineColor(j_xi+1);
    y_g[j_xi]->Draw("p");
  }
  
  leg->Draw();
  can->SaveAs("plots/y_dist.pdf");
  can->Clear();

}
