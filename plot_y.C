// code to compare y distributions for beta values
// dists used: (rho, delta) = (2, 0); 2 < xi < 50, |y| < 1.2
// CT14NNLO PDF, 7 TeV J/psi; 1e8 events generated
// beta = {1.5, 2, 2.5, 3}

void plot_y()
{
  int dataPsi = 0;
  // part 0: variables
  const int n_histo = 1;
  const int n_xi =  dataPsi == 1 ? 7 : 3;
  double mass = dataPsi == 1 ? 3.097 : 9.46;
  
  double y_max = 5;
  string varName = "";
  double val[n_histo] = {2.};
  string name[n_histo] = {"MC_res_Ups_lowxi.root"};
  double dataN = 4.5;

  double xi_min[n_xi], xi_max[n_xi];
  for(int i = 0; i < n_xi; i++) {
    if(dataPsi == 1) {
      xi_min[i] = (i+7.) / mass;
      xi_max[i] = (i+8.) / mass;}
    else if(dataPsi == 0) {
      xi_min[i] = (i+10.)/mass;
      xi_max[i] = (i+11.)/mass;
    }
  }

  double yax_max[n_histo] = {2};
  double yax_min[n_histo] = {2e-3};

  double yax_max_lin[n_histo] = {0.4};
  double yax_min_lin[n_histo] = {0};

  int base_MC = 1;
  
  // part 1: get the MC information and fill histograms for each y bin
  
  // making the histograms

  TH1F ***y_h = new TH1F**[n_histo];
  for(int i = 0; i < n_histo; i++) {
    y_h[i] = new TH1F*[n_xi];
    for(int j = 0; j < n_xi; j++)
      y_h[i][j] = new TH1F(Form("y_xi%d_var%d", j, i), Form("y_xi%d_var%d", j, i), 100, 0, y_max); 
  }

  // opening the files

  // name of tree branches
  Double_t xi;
  Double_t w_qg;
  Double_t w_gg;
  Double_t w_cos;
  Double_t w_sstar;
  Double_t y;

  Double_t w_xi;
  Double_t w_sig;
  
  // cycle over all files
  for(int j = 0; j < n_histo; j++)
    {
      TFile *fin = new TFile(Form("%s", name[j].c_str()));
      TTree *tree = (TTree*)fin->Get("qqbar");

      if(base_MC == 1) {
	tree->SetBranchAddress("w_gg", &w_gg);
	tree->SetBranchAddress("w_qg", &w_qg);
	tree->SetBranchAddress("w_cos", &w_cos);
	tree->SetBranchAddress("w_sstar", &w_sstar);
      }
      else if (base_MC == 0) {
	tree->SetBranchAddress("w_sig", &w_sig);
	tree->SetBranchAddress("w_xi", &w_xi);
      }
      tree->SetBranchAddress("xi", &xi);
      tree->SetBranchAddress("y", &y);
      
      Int_t nentries = (Int_t)tree->GetEntries();
      cout << nentries << " events in tree " << j+1 << endl;
      int chk = nentries / 100;
      for( Int_t i = 0; i < nentries; i++)
	{
	  tree->GetEntry(i);
	  
	  for (int k = 0; k < n_xi; k++)
	    if(xi < xi_max[k] && xi > xi_min[k]) {
	      if(base_MC == 0) y_h[j][k]->Fill(abs(y), w_xi*w_sig);
	      else if(base_MC == 1) y_h[j][k]->Fill(abs(y), w_gg*w_cos);
	    }
	  
	  if((i+1)%chk == 0) {
	    cout << (i+1)/chk << "% | " << flush;
	  }
	}
      cout << endl << endl;
    }

  // part 2: read the data to plot
  TGraphErrors **g = new TGraphErrors*[n_xi];

  int n_val[n_xi];
  for (int i = 0; i < n_xi; i++) 
    n_val[i] = i == 6 ? 4 : i == 3 ? 6 : 5;
    
  if (dataPsi == 1)
    {
      ifstream file;
      string data;
      int n_pts[6] = {31, 12, 12, 12, 11, 9};
      double pt_lo[6][n_pts[0]], pt_hi[6][n_pts[0]], y_v[6][n_pts[0]], dy[6][n_pts[0]], sig[6][n_pts[0]], dsig[6][n_pts[0]];
      double nums[12];
  
      
      string filelist[6] = {"data/CMS_jpsi_7_cs_y1.txt",
			    "data/LHCb_jpsi_7_cs_y1.txt",
			    "data/LHCb_jpsi_7_cs_y2.txt",
			    "data/LHCb_jpsi_7_cs_y3.txt",
			    "data/LHCb_jpsi_7_cs_y4.txt",
			"data/LHCb_jpsi_7_cs_y5.txt"};
      for(int k = 0; k < 6; k++) {
	file.open(filelist[k].c_str());
	
	for(int j = 0; j < n_pts[k]; j++) {
	  for(int i = 0; i < 12; i++)
	    file >> nums[i];
	  
	  y_v[k][j] = 0.5 * (nums[1] + nums[0]);
	  dy[k][j] = 0.5 * (nums[1] - nums[0]);
	  
	  pt_lo[k][j] = nums[3];
	  pt_hi[k][j] = nums[4];
	  
	  double signorm = k < 1 ? 5.961e-2*1e3 : 1;
	  sig[k][j] = nums[5]*mass/(signorm);
	  dsig[k][j] = sqrt(nums[6]*nums[6]+nums[8]*nums[8])*mass/(signorm);
	}
	file.close();
	
	// normalize to first point
	double aux;
	if(k < 1) aux = sig[0][0]*dataN;
	for(int i = 0; i < n_pts[k]; i++) {
	  sig[k][i] /= aux;
	  dsig[k][i] /= aux;
	}
	
      }
      
      for(int i = 0; i < n_xi; i++) {
	double y_av[n_val[i]], sig_av[n_val[i]], dy_av[n_val[i]], dsig_av[n_val[i]];
	for(int j = 0; j < n_val[i]; j++) {
	  y_av[j] = i < 3 ? y_v[j+1][0] : y_v[j][0];
	  dy_av[j] = i < 3 ? dy[j+1][0] : dy[j][0];
	  
	  sig_av[j] = i < 3 ? sig[j+1][i+5] : sig[j][j < 1 ? i-3 : i+5];
	  dsig_av[j] = i < 3 ? dsig[j+1][i+5] : dsig[j][j < 1 ? i-3 : i+5];
	  
	}
	g[i] = new TGraphErrors(n_val[i], y_av, sig_av, dy_av, dsig_av); 
      }
    }
  
  else if (dataPsi == 0)
    {
      ifstream file;
      string data;
      int n_pts[6] = {22, 24, 25, 24, 21, 15};
      double pt_lo[6][n_pts[2]], pt_hi[6][n_pts[2]], y_v[6][n_pts[2]], dy[6][n_pts[2]], sig[6][n_pts[2]], dsig[6][n_pts[2]];
      double nums[12];
       
      string filelist[6] = {"data/CMS_ups1_7_cs_y1.txt",
			    "data/LHCb_ups1_7_cs_y1.txt",
			    "data/LHCb_ups1_7_cs_y2.txt",
			    "data/LHCb_ups1_7_cs_y3.txt",
			    "data/LHCb_ups1_7_cs_y4.txt",
			    "data/LHCb_ups1_7_cs_y5.txt"};

      for(int k = 0; k < 6; k++) {
	file.open(filelist[k].c_str());
	
	for(int j = 0; j < n_pts[k]; j++) {
	  for(int i = 0; i < 12; i++)
	    file >> nums[i];
	  
	  y_v[k][j] = 0.5 * (nums[1] + nums[0]);
	  dy[k][j] = 0.5 * (nums[1] - nums[0]);
	  
	  pt_lo[k][j] = nums[3];
	  pt_hi[k][j] = nums[4];
	  
	  double signorm = k < 1 ? 2.4*2.48e-2*1e6 : 2.48e-2*1e3;
	  sig[k][j] = nums[5]*mass/(signorm);
	  dsig[k][j] = sqrt(nums[6]*nums[6]+nums[8]*nums[8])*mass/(signorm);
	}
	file.close();
	
	// normalize to first point
	double aux;
	if(k < 1) aux = sig[0][0]*dataN;
	cout << k << " " << aux << endl;;
	for(int i = 0; i < n_pts[k]; i++) {
	  sig[k][i] /= aux;
	  dsig[k][i] /= aux;
	}
	
      }
      
      for(int i = 0; i < n_xi; i++) {
	double y_av[n_val[i]], sig_av[n_val[i]], dy_av[n_val[i]], dsig_av[n_val[i]];
	for(int j = 0; j < n_val[i]; j++) {
	  y_av[j] = y_v[j+1][0];
	  dy_av[j] = dy[j+1][0];
	  
	  sig_av[j] = sig[j+1][i+10];
	  cout << i << " " << j << " " << sig[j+1][i+10] << endl;
	  dsig_av[j] = dsig[j+1][i+10];

	}
	g[i] = new TGraphErrors(n_val[i], y_av, sig_av, dy_av, dsig_av);
      }
    }
    
  // part 3: draw the results (1 canvas for each beta)
  
  TCanvas *can = new TCanvas("", "", 700, 700);
  can->SetLogy();

  // scaling the histos to the CMS integral
  for(int i = 0; i < n_histo; i++) {
    double n = y_h[i][0]->Integral("width");
    cout << n << endl;
    for(int j = 0; j < n_xi; j++) 
      y_h[i][j]->Scale(1./n);
  }
  
  // double cycle for all 4 plots of n_xi = 8 histos each  
  for(int i_b = 0; i_b < n_histo; i_b++) {    
    y_h[i_b][0]->SetStats(0);
    y_h[i_b][0]->SetLineColor(1);
    y_h[i_b][0]->GetXaxis()->SetTitle("y");
    y_h[i_b][0]->GetYaxis()->SetTitle("Events");
    if(n_histo > 1)
      y_h[i_b][0]->SetTitle(Form("%.s = %.1f", varName.c_str(), val[i_b]));
    else
      y_h[i_b][0]->SetTitle("");
    y_h[i_b][0]->GetYaxis()->SetTitleOffset(1.3);
    y_h[i_b][0]->GetYaxis()->SetRangeUser(yax_min[i_b], yax_max[i_b]);
    y_h[i_b][0]->Draw("hist"); 
    
    TLegend *leg = new TLegend(0.65, 0.6, 0.9, 0.9);
    leg->SetTextSize(0.03);
    leg->AddEntry(y_h[i_b][0], Form("%.0f < p_{T} < %.0f", xi_min[0]*mass, xi_max[0]*mass), "l");

    for(int j_xi = 1; j_xi < n_xi; j_xi++) {
      y_h[i_b][j_xi]->SetLineColor(j_xi+1);
      y_h[i_b][j_xi]->Draw("hist same");
      leg->AddEntry(y_h[i_b][j_xi], Form("%.0f < p_{T} < %.0f", xi_min[j_xi]*mass, xi_max[j_xi]*mass), "l");
      
    }

    for(int j_xi = 0; j_xi < n_xi; j_xi++) {
      g[j_xi]->SetMarkerStyle(20);
      g[j_xi]->SetMarkerSize(.75);
      g[j_xi]->SetMarkerColor(j_xi+1);
      g[j_xi]->SetLineColor(j_xi+1);
      g[j_xi]->Draw("p");
    }
    
    leg->Draw();

    cout << "leg draw" << endl;
    //if(n_histo < 2)
      can->SaveAs("plots/y_dist_log.pdf");
    //else
    //  can->SaveAs(Form("plots/%s_%d_log.pdf", varName.c_str(), i_b));

    can->Clear();

  }

  // same cycle for linear scale
  can->SetLogy(0);

  // double cycle for all 4 plots of n_xi = 8 histos each  
  for(int i_b = 0; i_b < n_histo; i_b++) {    
    y_h[i_b][0]->SetStats(0);
    y_h[i_b][0]->SetLineColor(1);
    y_h[i_b][0]->GetXaxis()->SetTitle("y");
    y_h[i_b][0]->GetYaxis()->SetTitle("Events");
    if(n_histo > 1)
      y_h[i_b][0]->SetTitle(Form("%.s = %.1f", varName.c_str(), val[i_b]));
    else
      y_h[i_b][0]->SetTitle("");
y_h[i_b][0]->GetYaxis()->SetTitleOffset(1.3);
    y_h[i_b][0]->GetYaxis()->SetRangeUser(yax_min_lin[i_b], yax_max_lin[i_b]);
    y_h[i_b][0]->Draw("hist");

    TLegend *leg = new TLegend(0.65, 0.6, 0.9, 0.9);
    leg->SetTextSize(0.03);
    leg->AddEntry(y_h[i_b][0], Form("%.0f < p_{T} < %.0f", xi_min[0]*mass, xi_max[0]*mass), "l");
    
    for(int j_xi = 1; j_xi < n_xi; j_xi++) {
      y_h[i_b][j_xi]->SetLineColor(j_xi+1);
      y_h[i_b][j_xi]->Draw("hist same");
      leg->AddEntry(y_h[i_b][j_xi], Form("%.0f < p_{T} < %.0f", xi_min[j_xi]*mass, xi_max[j_xi]*mass), "l");
     }

    for(int j_xi = 0; j_xi < n_xi; j_xi++) {
      g[j_xi]->SetMarkerStyle(20);
      g[j_xi]->SetMarkerSize(.75);
      g[j_xi]->SetMarkerColor(j_xi+1);
      g[j_xi]->SetLineColor(j_xi+1);
      g[j_xi]->Draw("p");
    }
    
    leg->Draw();
    if(n_histo < 2)
      can->SaveAs("plots/y_dist.pdf");
    else
      can->SaveAs(Form("plots/%s_%d.pdf", varName.c_str(), i_b));
    can->Clear();
  }

}
