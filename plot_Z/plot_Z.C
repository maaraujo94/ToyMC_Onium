/* code to plot toy MC histos vs data graphs
   - make TGraphs (external code)
   - make THistos
   - plot both together
   process is repeated for both xi plots in y bins and y plots in xi bins
*/

#import "plot_codes.C"

void plot_Z()
{
  /////////////////////////////////////////
  // part 0 : auxiliary variables

  string loc = "/home/mariana/Documents/2020_PhD_work/Phenom/ToyMC_Onium/MC/";
  
  string dataName = "Z"; // "state" we're plotting
  const double mass = 91.1876;
  const int n_ybins[2][2] = {{7, 6}, {1, 6}};

  double yminA13[] = {0.0};
  double ymaxA13[] = {2.5};
  double yminC13[] = {0.0, 0.0, 0.4, 0.8, 1.2, 1.6};
  double ymaxC13[] = {2.4, 0.4, 0.8, 1.2, 1.6, 2.4};

  double yminA8[] = {0.0, 0.0, 0.4, 0.8, 1.2, 1.6, 2.0};
  double ymaxA8[] = {2.4, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};
  double yminC8[] = {0.0, 0.0, 0.4, 0.8, 1.2, 1.6};
  double ymaxC8[] = {2.0, 0.4, 0.8, 1.2, 1.6, 2.0};
  
  double xi_min_global = 0.1;
  double xi_max_global = 30;
  
  /////////////////////////////////////////
  // part 1 : variables to be changed on each run

  // type of MC files
  string typein = "ZMC_symm_pT/";
  string typeout = "ZMC_noSc_pT/";

  // choose sqrt(s) we're plotting
  const int nsqs = 2;
  string sqsNames[nsqs] = {"8", "13"}; 
  int nentries[nsqs];
  double jacV[nsqs];
  
  double xi_bin_width = 0.1;   // constant for easier normalization
  
  double normfactor[nsqs] = {40., 40.}; // same for all y and sqrt(s)
  
  // data TGraphs (xi and y dists)
  TGraphAsymmErrors ***xi_g_A = new TGraphAsymmErrors**[nsqs];
  TH1F ***xi_h_A = new TH1F**[nsqs]; // xi histogram, defined before the data cycle
  TGraphAsymmErrors ***xi_g_C = new TGraphAsymmErrors**[nsqs];
  TH1F ***xi_h_C = new TH1F**[nsqs]; // xi histogram, defined before the data cycle
  for(int i = 0; i < nsqs; i++){
    xi_g_A[i] = new TGraphAsymmErrors*[n_ybins[i][0]];
    xi_h_A[i] = new TH1F*[n_ybins[i][0]];
    xi_g_C[i] = new TGraphAsymmErrors*[n_ybins[i][1]];
    xi_h_C[i] = new TH1F*[n_ybins[i][1]];
  }

  double xi_min = 1e-1, xi_max = 30; // limits for plotting the xi dists (x-axis)
  double xi_yax_min = 1e-7, xi_yax_max = 1e3; // limits for plotting the xi dists (y-axis)

  /////////////////////////////////////////
  // part 2 : create ROOT file and save results

  TFile *tf = new TFile(Form("%sMC_vs_Data_%s.root", typeout.c_str(), dataName.c_str()), "RECREATE", "MC histograms and Data graphs");
  TTree *norm_factor = new TTree("norm", "norm factor");

  norm_factor->Branch("normfactor_8", &normfactor[0]);
  norm_factor->Branch("normfactor_13", &normfactor[1]);
  norm_factor->Fill();
  
  tf->Write();
  tf->Close();    
  
  // cycle in sqrt(s)
  for(int i = 0; i < nsqs; i++) {
    
    string sqsName = sqsNames[i];
    string fName = Form("%s%sMC_res_%s_%s", loc.c_str(), typein.c_str(), dataName.c_str(), sqsName.c_str());   // file to open

    /////////////////////////////////////////
    // part 3 : making the data TGraphs

    if(sqsName == "8") {
      xi_Z_8_A(xi_g_A, xi_h_A);
      xi_Z_8_C(xi_g_C, xi_h_C);
    }
    else if(sqsName == "13") {
      xi_Z_13_A(xi_g_A, xi_h_A);
      xi_Z_13_C(xi_g_C, xi_h_C);
    }
    else {
      cout << "incorrect sqsName formatting!" << endl;
      break;
    }
    
    /////////////////////////////////////////
    // part 4 : making the MC histograms
    
    // variables to open the MC
    Double_t xi, y;
    Double_t jac, w_ggHX_tr;
    Double_t w_ug, w_dg, w_uub, w_ddb, w_ubg, w_dbg;
    TLorentzVector* muonP = 0;
    TLorentzVector* muonN = 0;
    
    TFile *fin = new TFile(Form("%s.root", fName.c_str()));
    TTree *tree = (TTree*)fin->Get("ZMC");
    
    tree->SetBranchAddress("xi",        &xi);
    tree->SetBranchAddress("y",         &y);
    tree->SetBranchAddress("jac",       &jac);
    tree->SetBranchAddress("w_ggHX_tr", &w_ggHX_tr);
    tree->SetBranchAddress("w_ug",      &w_ug);
    tree->SetBranchAddress("w_dg",      &w_dg);
    tree->SetBranchAddress("w_uub",     &w_uub);
    tree->SetBranchAddress("w_ddb",     &w_ddb);
    tree->SetBranchAddress("w_ubg",     &w_ubg);
    tree->SetBranchAddress("w_dbg",     &w_dbg);
    tree->SetBranchAddress("muonP",     &muonP);
    tree->SetBranchAddress("muonN",     &muonN);
    
    nentries[i] = tree->GetEntries();
    cout << nentries[i] << " events in " << sqsName << " TeV tree " << endl;
    int chk = nentries[i] / 100;
    jacV[i] = 0;

    if(sqsName == "13") {
      for( int in = 0; in < nentries[i]; in++)
	{
	  tree->GetEntry(in);
	  jacV[i] += jac;
  
	  // ATLAS part
	  for (int k = 0; k < n_ybins[i][0]; k++) 
	    if(abs(y) > yminA13[k] && abs(y) < ymaxA13[k] &&
	       xi > xi_min_global && xi < xi_max_global &&
	       muonP->Pt() > 27 && muonN->Pt() > 27 &&
	       abs(muonP->Eta()) < 2.5 && abs(muonN->Eta()) < 2.5) {
	      double w_sum = w_ug+w_dg+w_uub+w_ddb+w_ubg+w_dbg;
	      xi_h_A[i][k]->Fill(xi, jac*w_ggHX_tr*w_sum);
	    }
	  
	  // CMS part
	  for (int k = 0; k < n_ybins[i][1]; k++) 
	    if(abs(y) > yminC13[k] && abs(y) < ymaxC13[k] &&
	       xi > xi_min_global && xi < xi_max_global &&
	       muonP->Pt() > 25 && muonN->Pt() > 25 &&
	       abs(muonP->Eta()) < 2.4 && abs(muonN->Eta()) < 2.4) {
	      double w_sum = w_ug+w_dg+w_uub+w_ddb+w_ubg+w_dbg;
	      xi_h_C[i][k]->Fill(xi, jac*w_ggHX_tr*w_sum);
	    }
	  
	  if((in+1)%chk == 0) {
	    cout << (in+1)/chk << "% | " << flush;
	  }
	}
    }
    else if(sqsName == "8") {
      for( int in = 0; in < nentries[i]; in++)
	{
	  tree->GetEntry(in);
	  jacV[i] += jac;
	  
	  // ATLAS part
	  for (int k = 0; k < n_ybins[i][0]; k++) 
	    if(abs(y) > yminA8[k] && abs(y) < ymaxA8[k] &&
	       xi > xi_min_global && xi < xi_max_global &&
	       muonP->Pt() > 20 && muonN->Pt() > 20 &&
	       abs(muonP->Eta()) < 2.4 && abs(muonN->Eta()) < 2.4) {
	      double w_sum = w_ug+w_dg+w_uub+w_ddb+w_ubg+w_dbg;
	      xi_h_A[i][k]->Fill(xi, jac*w_ggHX_tr*w_sum);
	    }
	  
	  // CMS part
	  for (int k = 0; k < n_ybins[i][1]; k++) 
	    if(abs(y) > yminC8[k] && abs(y) < ymaxC8[k] &&
	       xi > xi_min_global && xi < xi_max_global &&

	       ((muonP->Pt() > 25 && muonN->Pt() > 10 &&
		 abs(muonP->Eta()) < 2.1 && abs(muonN->Eta()) < 2.4) ||
		(muonP->Pt() > 10 && muonN->Pt() > 25 &&
		 abs(muonP->Eta()) < 2.4 && abs(muonN->Eta()) < 2.1))  ){
	      double w_sum = w_ug+w_dg+w_uub+w_ddb+w_ubg+w_dbg;
	      xi_h_C[i][k]->Fill(xi, jac*w_ggHX_tr*w_sum);
	    }
	  
	  if((in+1)%chk == 0) {
	    cout << (in+1)/chk << "% | " << flush;
	  }
	}
    }
    else {
      cout << "sqsName formatted incorrectly!" << endl;
      break;
    }
    
    cout << endl;
    
    fin->Close();

    // scaling the histos to nr events + y bin width
    double n = 1.;//jacV[i];
    for(int j = 0; j < n_ybins[i][0]; j++) { // ATLAS cycle
      double yM = (sqsName == "13" ? ymaxA13[j] : ymaxA8[j]);
      double ym = (sqsName == "13" ? yminA13[j] : yminA8[j]);
      xi_h_A[i][j]->Scale(normfactor[i]/(n*2.*(yM-ym)));
      // scaling each xi bin by its width
      int nb = xi_h_A[i][j]->GetNbinsX();
      for(int bin = 0; bin < nb; bin++) {
	double wd = xi_h_A[i][j]->GetBinWidth(bin+1) * mass;
	double vl = xi_h_A[i][j]->GetBinContent(bin+1);
	double ve = xi_h_A[i][j]->GetBinError(bin+1);
	xi_h_A[i][j]->SetBinContent(bin+1, vl/wd);
	xi_h_A[i][j]->SetBinError(bin+1, ve/wd);
	if(xi_h_A[i][j]->GetBinLowEdge(bin+1) < xi_min_global)
	  xi_h_A[i][j]->SetBinContent(bin+1, 0);
      }
    }
    for(int j = 0; j < n_ybins[i][1]; j++) { // CMS cycle
      double yM = (sqsName == "13" ? ymaxC13[j] : ymaxC8[j]);
      double ym = (sqsName == "13" ? yminC13[j] : yminC8[j]);
      xi_h_C[i][j]->Scale(normfactor[i]/(n*2.*(yM-ym)));
      // scaling each xi bin by its width
      int nb = xi_h_C[i][j]->GetNbinsX();
      for(int bin = 0; bin < nb; bin++) {
	double wd = xi_h_C[i][j]->GetBinWidth(bin+1) * mass;
	double vl = xi_h_C[i][j]->GetBinContent(bin+1);
	double ve = xi_h_C[i][j]->GetBinError(bin+1);
	xi_h_C[i][j]->SetBinContent(bin+1, vl/wd);
	xi_h_C[i][j]->SetBinError(bin+1, ve/wd);
	if(xi_h_C[i][j]->GetBinLowEdge(bin+1) < xi_min_global)
	  xi_h_C[i][j]->SetBinContent(bin+1, 0);
      }
    }
    
  }
  
  /////////////////////////////////////////
  // part 5 : plotting
  TFile *tf_h = new TFile(Form("%sMC_vs_Data_%s.root", typeout.c_str(), dataName.c_str()), "UPDATE", "MC histograms and Data graphs");

  TCanvas *can = new TCanvas("", "", 700, 700);
  can->SetLogy();
  can->SetLogx();
  
  // xi plots
  
  // cycle for all the sqrt(s) and y bins (histo + graph)
  // repeat for CMS and ATLAS
  for(int is = 0; is < nsqs; is++) {
    // ATLAS cycle
    for(int j_y = 0; j_y < n_ybins[is][0]; j_y++) {
      
      double yM = (sqsNames[is] == "13" ? ymaxA13[j_y] : ymaxA8[j_y]);
      double ym = (sqsNames[is] == "13" ? yminA13[j_y] : yminA8[j_y]);


      // plot frame + title etc
      TH1F *fc = can->DrawFrame(xi_min, xi_yax_min, xi_max, xi_yax_max);
      fc->SetXTitle("#xi");
      fc->SetYTitle("d#sigma / dp_{T}dy (pb/GeV)");
      fc->GetYaxis()->SetTitleOffset(1.3);
      fc->SetTitle(Form("%s TeV ATLAS %.1f < | y | < %.1f", sqsNames[is].c_str(), ym, yM));
      can->Modified();
      can->SetTitle("");
      
      // draw legend
      TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
      leg->SetTextSize(0.03);
      
      // draw histogram
      xi_h_A[is][j_y]->SetStats(0);
      xi_h_A[is][j_y]->SetLineColor(1);
      xi_h_A[is][j_y]->Draw("hist same");
      //xi_h_A[is][j_y]->Draw("error same");
      
      //save in ROOT fle
      xi_g_A[is][j_y]->SetName(Form("%s_g", xi_h_A[is][j_y]->GetTitle()));
      xi_h_A[is][j_y]->SetTitle(Form("%s TeV ATLAS %.1f < | y | < %.1f", sqsNames[is].c_str(), ym, yM));
      xi_h_A[is][j_y]->Write();
    
      leg->AddEntry(xi_h_A[is][j_y], Form("MC"), "l");
	
      // draw data graph
      xi_g_A[is][j_y]->SetMarkerStyle(20);
      xi_g_A[is][j_y]->SetMarkerSize(.75);
      xi_g_A[is][j_y]->Draw("p same");
      // save in ROOT file
      xi_g_A[is][j_y]->SetTitle(Form("%s TeV ATLAS %.1f < y < %.1f", sqsNames[is].c_str(), ym, yM));
      xi_g_A[is][j_y]->Write();
      
      leg->AddEntry(xi_g_A[is][j_y], "data", "pl");
	
      // save plot
      leg->Draw();
      can->SaveAs(Form("%splots/xi_%s_A_y%d.pdf", typeout.c_str(), sqsNames[is].c_str(), j_y));
      can->Clear();
    }
    // CMS cycle
    for(int j_y = 0; j_y < n_ybins[is][1]; j_y++) {
      
      double yM = (sqsNames[is] == "13" ? ymaxC13[j_y] : ymaxC8[j_y]);
      double ym = (sqsNames[is] == "13" ? yminC13[j_y] : yminC8[j_y]);


      // plot frame + title etc
      TH1F *fc = can->DrawFrame(xi_min, xi_yax_min, xi_max, xi_yax_max);
      fc->SetXTitle("#xi");
      fc->SetYTitle("d#sigma / dp_{T}dy (pb/GeV)");
      fc->GetYaxis()->SetTitleOffset(1.3);
      fc->SetTitle(Form("%s TeV CMS %.1f < | y | < %.1f", sqsNames[is].c_str(), ym, yM));
      can->Modified();
      can->SetTitle("");
      
      // draw legend
      TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
      leg->SetTextSize(0.03);
      
      // draw histogram
      xi_h_C[is][j_y]->SetStats(0);
      xi_h_C[is][j_y]->SetLineColor(1);
      xi_h_C[is][j_y]->Draw("hist same");
      //xi_h_C[is][j_y]->Draw("error same");
      
      //save in ROOT fle
      xi_g_C[is][j_y]->SetName(Form("%s_g", xi_h_C[is][j_y]->GetTitle()));
      xi_h_C[is][j_y]->SetTitle(Form("%s TeV CMS %.1f < | y | < %.1f", sqsNames[is].c_str(), ym, yM));
      xi_h_C[is][j_y]->Write();
    
      leg->AddEntry(xi_h_C[is][j_y], Form("MC"), "l");
	
      // draw data graph
      xi_g_C[is][j_y]->SetMarkerStyle(20);
      xi_g_C[is][j_y]->SetMarkerSize(.75);
      xi_g_C[is][j_y]->Draw("p same");
      // save in ROOT file
      xi_g_C[is][j_y]->SetTitle(Form("%s TeV CMS %.1f < y < %.1f", sqsNames[is].c_str(), ym, yM));
      xi_g_C[is][j_y]->Write();
      
      leg->AddEntry(xi_g_C[is][j_y], "data", "pl");
	
      // save plot
      leg->Draw();
      can->SaveAs(Form("%splots/xi_%s_C_y%d.pdf", typeout.c_str(), sqsNames[is].c_str(), j_y));
      can->Clear();
    }
  }
      
  can->Destructor();

  tf_h->Write();
  tf_h->Close();

}
