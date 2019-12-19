// code to compare xi distributions for beta values
// dists used: (rho, delta) = (2, 0); 2 < xi < 50, |y| < 1.2
// CT14NNLO PDF, 7 TeV J/psi; 1e8 events generated
// beta = {1.5, 2, 2.5, 3}

void plot_xi()
{
  // part 0 : variables
  const int n_y = 6;
  string name = "MC_res_Ups_lowxi.root";
  double dataN = 1./1000.;
  double histoN = 2000.;
  
  double y_min[n_y] = {0, 2, 2.5, 3, 3.5, 4};
  double y_max[n_y] = {1.2, 2.5, 3, 3.5, 4, 4.5};
  
  double xi_max[n_y] = {40, 11, 11, 11, 11, 11};
  double xi_min[n_y] = {-1, -1, -1, -1, -1, -1};

  double yax_min[n_y] = {1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5};
  double yax_max[n_y] = {1e4, 1e4, 1e4, 1e4, 1e4, 1e4};

  double bin_width = 0.1;
  int n_bins[n_y];

  for(int i = 0; i < n_y; i++) {
    n_bins[i] = (int)((xi_max[i] - xi_min[i]) / bin_width);
  }

  int base_MC = 1;
  int dataPsi = 0;
  double mass = dataPsi == 1 ? 3.097 : 9.46;
  
  // part 1: get the MC information and fill histograms for each y bin
  
  // making the histograms  
  TH1F **xi_h = new TH1F*[n_y];
  for(int j = 0; j < n_y; j++)
    xi_h[j] = new TH1F(Form("xi_y%d", j), Form("%s %.1f < |y| < %.1f", j < 1 ? "CMS" : "LHCb", y_min[j], y_max[j]), n_bins[j], xi_min[j], xi_max[j]); 
  

  // opening the files
  
  Double_t xi;
  Double_t w_qg;
  Double_t w_gg;
  Double_t w_cos;
  Double_t w_sstar;
  
  Double_t w_xi;
  Double_t w_sig;
  Double_t y;

  TFile *fin = new TFile(Form("%s", name.c_str()));
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
  cout << nentries << " events in tree " << endl;
  int chk = nentries / 100;
  for( Int_t i = 0; i < nentries; i++)
    {
      tree->GetEntry(i);

      for (int k = 0; k < n_y; k++)
	if(abs(y) < y_max[k] && abs(y) > y_min[k] && xi > 1 && xi < 50) {
	  if(base_MC == 0) xi_h[k]->Fill(xi, w_xi * w_sig);
	  else if(base_MC == 1) xi_h[k]->Fill(xi, w_gg * w_cos);
	}
      
      if((i+1)%chk == 0) {
	cout << (i+1)/chk << "% | " << flush;
      }
    }
  cout << endl;

  // part 2: read the data to plot
  TGraphAsymmErrors **g = new TGraphAsymmErrors*[n_y];
  
  if (dataPsi == 1) {
    ifstream file;
    string data;
    int n_pts;
    double nums[12];
    
    string filelist[n_y] = {"data/CMS_jpsi_7_cs_y1.txt",
			    "data/LHCb_jpsi_7_cs_y1.txt",
			    "data/LHCb_jpsi_7_cs_y2.txt",
			    "data/LHCb_jpsi_7_cs_y3.txt",
			    "data/LHCb_jpsi_7_cs_y4.txt",
			    "data/LHCb_jpsi_7_cs_y5.txt"};
    for(int k = 0; k < n_y; k++) {
      file.open(filelist[k].c_str());
      n_pts = 0;
      while(getline(file,data)) n_pts+=1;
      
      double pt[n_pts], dpt_lo[n_pts], dpt_hi[n_pts], sig[n_pts], dsig[n_pts];
      
      file.clear();
      file.seekg(0);
      for(int j = 0; j < n_pts; j++) {
	for(int i = 0; i < 12; i++)
	  file >> nums[i];
	
	pt[j] = nums[2]/mass;
	dpt_lo[j] = (nums[2]-nums[3])/mass;
	dpt_hi[j] = (nums[4]-nums[2])/mass;
	
	double signorm = k < 1 ? 5.961e-2*1e3 : 1;
	sig[j] = nums[5]*mass/(signorm);
	dsig[j] = sqrt(nums[6]*nums[6]+nums[8]*nums[8])*mass/(signorm);
      }
      file.close();
      
      // normalize to first point
      double aux;
      if(k < 1) aux = sig[0]*dataN;
      for(int i = 0; i < n_pts; i++) {
	sig[i] /= aux;
	dsig[i] /= aux;
      }
      
      g[k] = new TGraphAsymmErrors(n_pts, pt, sig, dpt_lo, dpt_hi, dsig, dsig);
      
    }

  }
  else if (dataPsi == 0) {
    ifstream file;
    string data;
    int n_pts;
    double nums[12];
    
    string filelist[n_y] = {"data/CMS_ups1_7_cs_y1.txt",
			    "data/LHCb_ups1_7_cs_y1.txt",
			    "data/LHCb_ups1_7_cs_y2.txt",
			    "data/LHCb_ups1_7_cs_y3.txt",
			    "data/LHCb_ups1_7_cs_y4.txt",
			    "data/LHCb_ups1_7_cs_y5.txt"};

    for(int k = 0; k < n_y; k++) {
      file.open(filelist[k].c_str());
      n_pts = 0;
      while(getline(file,data)) n_pts+=1;
      
      double pt[n_pts], dpt_lo[n_pts], dpt_hi[n_pts], sig[n_pts], dsig[n_pts];
      
      file.clear();
      file.seekg(0);
      for(int j = 0; j < n_pts; j++) {
	for(int i = 0; i < 12; i++)
	  file >> nums[i];
	
	pt[j] = nums[2]/mass;
	dpt_lo[j] = (nums[2]-nums[3])/mass;
	dpt_hi[j] = (nums[4]-nums[2])/mass;
	
	double signorm = k < 1 ? 2.4*2.48e-2*1e6 : 2.48e-2*1e3;
	sig[j] = nums[5]*mass/(signorm);
	dsig[j] = sqrt(nums[6]*nums[6]+nums[8]*nums[8])*mass/(signorm);
      }
      file.close();
      
      // normalize to first point
      double aux;
      if(k < 1) aux = sig[0]*dataN;
      for(int i = 0; i < n_pts; i++) {
	sig[i] /= aux;
	dsig[i] /= aux;
      }
      
      g[k] = new TGraphAsymmErrors(n_pts, pt, sig, dpt_lo, dpt_hi, dsig, dsig);
      
    }
    
  }
  
  // part 3: draw the results (1 canvas for each y bin)
  
  TCanvas *can = new TCanvas("", "", 700, 700);
  can->SetLogy();
  
  // scaling the histos to the CMS integral
  double n = xi_h[0]->Integral("width")/histoN;
  for(int j = 0; j < n_y; j++) {
    xi_h[j]->Scale(j < 1 ? 1./(2.4*n) : 1./n);
  }
  
  // double cycle for all 6 plots of 3 histos each  
  for(int j_y = 0; j_y < n_y; j_y++) {    
    xi_h[j_y]->SetStats(0);
    xi_h[j_y]->SetLineColor(1);
    xi_h[j_y]->GetXaxis()->SetTitle("#xi");
    xi_h[j_y]->GetYaxis()->SetTitle(Form("Events/%.1f", bin_width));
    xi_h[j_y]->GetYaxis()->SetTitleOffset(1.3);
    xi_h[j_y]->GetYaxis()->SetRangeUser(yax_min[j_y], yax_max[j_y]);
    xi_h[j_y]->Draw("hist");
    
    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->SetTextSize(0.03);
    leg->AddEntry(xi_h[j_y], Form("MC"), "l");
    
    g[j_y]->SetMarkerStyle(20);
    g[j_y]->SetMarkerSize(.75);
    g[j_y]->Draw("p");
    leg->AddEntry(g[j_y], "data", "pl");
    
    leg->Draw();
    can->SaveAs(Form("plots/xi_new_%s.pdf", j_y < 1? "CMS" : Form("LHCb_y%d", j_y)));
    can->Clear();
  }
}
