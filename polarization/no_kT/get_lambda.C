// code to determine the lambda parameters from the TProfiles of angular quantities
// reads TProfiles from ROOT file

double lambda_theta(double cos2th)
{
  return (1.-3.*cos2th)/(cos2th-3./5.);
  
}
double e_lth(double cos2th, double ec)
{
  return 20./pow(5.*cos2th-3, 2)*ec;
}

double lambda_phi(double cos2phi, double lth)
{
  return (3.+lth)*cos2phi;
}
double e_lp(double cos2phi, double ec, double lth, double elth)
{
  double p1 = cos2phi*elth;
  double p2 = (3.+lth)*ec;

  return sqrt(p1*p1+p2*p2);
}

double lambda_thetaphi(double sc, double lth)
{
  return 5./4.*(3.+lth)*sc;
}
double e_ltp(double sc, double esc, double lth, double elth)
{
  double p1 = 5./4.*sc*elth;
  double p2 = 5./4.*(3.+lth)*esc;

  return sqrt(p1*p1+p2*p2);
}

double lambda_tilde(double lth, double lph)
{
  return (lth+3.*lph)/(1.-lph);
}

void get_lambda()
{ 
  string dataName = "jpsi", sqsName = "7";

  TFile *infile = new TFile(Form("ang_profs/plots_%s.root", dataName.c_str()));

  int y_bins = 4;
  double ylims[y_bins+1];
  for(int i = 0; i <= y_bins; i++)
    ylims[i] = i;
  
  TProfile ***costh = new TProfile**[4];
  TProfile ***phi = new TProfile**[4];
  TProfile ***phith = new TProfile**[4];
  for(int i = 0; i < 4; i++) {
    costh[i] = new TProfile*[y_bins];
    phi[i] = new TProfile*[y_bins];
    phith[i] = new TProfile*[y_bins];
    for(int j = 0; j < y_bins; j++) {
      costh[i][j] = (TProfile*)infile->Get(Form("h_ct%d_y%d", i, j));
      phi[i][j] = (TProfile*)infile->Get(Form("h_phi%d_y%d", i, j));
      phith[i][j] = (TProfile*)infile->Get(Form("h_pt%d_y%d", i, j));
    }
  }

  int pt_bins = costh[0][0]->GetNbinsX();

  TGraphErrors ***g_lth = new TGraphErrors**[4];
  TGraphErrors ***g_lp = new TGraphErrors**[4];
  TGraphErrors ***g_ltp = new TGraphErrors**[4];
  TGraphErrors ***g_ltil = new TGraphErrors**[4];

  for(int i_pol = 0; i_pol < 4; i_pol++) {
    g_lth[i_pol] = new TGraphErrors*[y_bins];
    g_lp[i_pol] = new TGraphErrors*[y_bins];
    g_ltp[i_pol] = new TGraphErrors*[y_bins];
    g_ltil[i_pol] = new TGraphErrors*[y_bins];
    
    for(int i_y = 0; i_y < y_bins; i_y++) {
      double lth[pt_bins], elth[pt_bins], lp[pt_bins], elp[pt_bins], ltp[pt_bins], eltp[pt_bins], ltil[pt_bins], eltil[pt_bins];
      double pt[pt_bins], ept[pt_bins];
      
      for(int i_pt = 0; i_pt < pt_bins; i_pt++) {
	double pt_min = costh[i_pol][i_y]->GetXaxis()->GetBinLowEdge(i_pt+1);
	double pt_max = costh[i_pol][i_y]->GetXaxis()->GetBinUpEdge(i_pt+1);
	pt[i_pt] = 0.5*(pt_max+pt_min);
	ept[i_pt] = 0.5*(pt_max-pt_min);

	double cos2th = costh[i_pol][i_y]->GetBinContent(i_pt+1);
	double ecth = costh[i_pol][i_y]->GetBinError(i_pt+1);
	lth[i_pt] = lambda_theta(cos2th);
	elth[i_pt] = e_lth(cos2th, ecth);

	double cos2phi = phi[i_pol][i_y]->GetBinContent(i_pt+1);
	double ecp = phi[i_pol][i_y]->GetBinError(i_pt+1);
	lp[i_pt] = lambda_phi(cos2phi, lth[i_pt]);
	elp[i_pt] = e_lp(cos2phi, ecp, lth[i_pt], elth[i_pt]);

      	double sc = phith[i_pol][i_y]->GetBinContent(i_pt+1);
	double esc = phith[i_pol][i_y]->GetBinError(i_pt+1);
	ltp[i_pt] = lambda_thetaphi(sc, lth[i_pt]);
	eltp[i_pt] = e_ltp(sc, esc, lth[i_pt], elth[i_pt]);

	ltil[i_pt] = lambda_tilde(lth[i_pt], lp[i_pt]);
	eltil[i_pt] = 0.;
	
      }

      g_lth[i_pol][i_y] = new TGraphErrors(pt_bins, pt, lth, ept, elth);
      g_lth[i_pol][i_y]->SetLineColor(i_y+1);
      g_lth[i_pol][i_y]->SetMarkerColor(i_y+1);
      if(i_pol%2 == 1) g_lth[i_pol][i_y]->SetLineStyle(kDashed);
      g_lth[i_pol][i_y]->SetMarkerStyle(i_pol%2 == 0 ? 20 : 24);
      //if(i_pol%2 == 0)
      g_lth[i_pol][i_y]->SetMarkerSize(.75);
      
      g_lp[i_pol][i_y] = new TGraphErrors(pt_bins, pt, lp, ept, elp);
      g_lp[i_pol][i_y]->SetLineColor(i_y+1);
      g_lp[i_pol][i_y]->SetMarkerColor(i_y+1);
      if(i_pol%2 == 1) g_lp[i_pol][i_y]->SetLineStyle(kDashed);
      g_lp[i_pol][i_y]->SetMarkerStyle(i_pol%2 == 0 ? 20 : 24);
      //if(i_pol%2 == 0)
      g_lp[i_pol][i_y]->SetMarkerSize(.75);

      g_ltp[i_pol][i_y] = new TGraphErrors(pt_bins, pt, ltp, ept, eltp);
      g_ltp[i_pol][i_y]->SetLineColor(i_y+1);
      g_ltp[i_pol][i_y]->SetMarkerColor(i_y+1);
      if(i_pol%2 == 1) g_ltp[i_pol][i_y]->SetLineStyle(kDashed);
      g_ltp[i_pol][i_y]->SetMarkerStyle(i_pol%2 == 0 ? 20 : 24);
      //if(i_pol%2 == 0)
      g_ltp[i_pol][i_y]->SetMarkerSize(.75);

      g_ltil[i_pol][i_y] = new TGraphErrors(pt_bins, pt, ltil, ept, eltil);
      g_ltil[i_pol][i_y]->SetLineColor(i_y+1);
      g_ltil[i_pol][i_y]->SetMarkerColor(i_y+1);
      if(i_pol%2 == 1) g_ltil[i_pol][i_y]->SetLineStyle(kDashed);
      g_ltil[i_pol][i_y]->SetMarkerStyle(i_pol%2 == 0 ? 20 : 24);
      //if(i_pol%2 == 0)
      g_ltil[i_pol][i_y]->SetMarkerSize(.75);

    }
  }

  TCanvas *c = new TCanvas("","",700,700);
  string polN[2] = {"trans", "long"};
  string fullN[2] = {"transverse", "longitudinal"};
  
  TLegend *leg = new TLegend(0.6, 0.6, 0.9, 0.9);
  leg->SetTextSize(0.03);
  for(int i = 0; i < y_bins; i++) {
    leg->AddEntry(g_lth[0][i], Form("%.0f < |y| < %.0f ppHX", ylims[i], ylims[i+1]), "pl");
    leg->AddEntry(g_lth[1][i], Form("%.0f < |y| < %.0f ggHX", ylims[i], ylims[i+1]), "pl");
  }
  
  double thLim[4] = {0, 2.0, -1.5, 0};
  for(int i_pol = 0; i_pol < 4; i_pol++) {
    if(i_pol%2 == 0) {
      TH1F *fc = c->DrawFrame(0, thLim[i_pol], 40, thLim[i_pol+1]);
      fc->SetXTitle("#xi");
      fc->SetYTitle("#lambda_{#theta}");
      fc->GetYaxis()->SetTitleOffset(1.3);
      fc->SetTitle(Form("%s TeV #lambda_{#theta} %s", sqsName.c_str(), fullN[i_pol/2].c_str()));
    }
    
    for(int i_y = 0; i_y < y_bins; i_y++) {
      g_lth[i_pol][i_y]->Draw("psame");
    }

    if(i_pol%2 == 1) {
      leg->Draw();
      c->SaveAs(Form("lambda_plots/%s_l_theta_%s_%s.pdf", dataName.c_str(), polN[(i_pol-1)/2].c_str(), sqsName.c_str()));
      c->Clear();
    }
  }

  for(int i_pol = 0; i_pol < 4; i_pol++) {
    if(i_pol%2 == 0) {
      TH1F *fc = c->DrawFrame(0, -1, 40, 1);
      fc->SetXTitle("#xi");
      fc->SetYTitle("#lambda_{#phi}");
      fc->GetYaxis()->SetTitleOffset(1.3);
      fc->SetTitle(Form("%s TeV #lambda_{#phi} %s", sqsName.c_str(), fullN[i_pol/2].c_str()));
    }
    
    for(int i_y = 0; i_y < y_bins; i_y++) {
      g_lp[i_pol][i_y]->Draw("psame");
    }

    if(i_pol%2 == 1) {
      leg->Draw();
      c->SaveAs(Form("lambda_plots/%s_l_phi_%s_%s.pdf", dataName.c_str(), polN[(i_pol-1)/2].c_str(), sqsName.c_str()));
      c->Clear();
    }
  }

  for(int i_pol = 0; i_pol < 4; i_pol++) {
    if(i_pol%2 == 0) {
      TH1F *fc = c->DrawFrame(0, -1, 40, 1);
      fc->SetXTitle("#xi");
      fc->SetYTitle("#lambda_{#theta#phi}");
      fc->GetYaxis()->SetTitleOffset(1.3);
      fc->SetTitle(Form("%s TeV #lambda_{#theta#phi} %s", sqsName.c_str(), fullN[i_pol/2].c_str()));
    }
    
    for(int i_y = 0; i_y < y_bins; i_y++) {
      g_ltp[i_pol][i_y]->Draw("psame");
    }

    if(i_pol%2 == 1) {
      leg->Draw();
      c->SaveAs(Form("lambda_plots/%s_l_thetaphi_%s_%s.pdf", dataName.c_str(), polN[(i_pol-1)/2].c_str(), sqsName.c_str()));
      c->Clear();
    }
  }
    
  for(int i_pol = 0; i_pol < 4; i_pol++) {
    if(i_pol%2 == 0) {
      TH1F *fc = c->DrawFrame(0, thLim[i_pol], 40, thLim[i_pol+1]);
      //TH1F *fc = c->DrawFrame(0, -1, 40, 1);
      fc->SetXTitle("#xi");
      fc->SetYTitle("#tilde#lambda");
      fc->GetYaxis()->SetTitleOffset(1.3);
      fc->SetTitle(Form("%s TeV #tilde{#lambda} %s", sqsName.c_str(), fullN[i_pol/2].c_str()));
    }
    
    for(int i_y = 0; i_y < y_bins; i_y++) {
      g_ltil[i_pol][i_y]->Draw("psame");
    }

    if(i_pol%2 == 1) {
      leg->Draw();
      c->SaveAs(Form("lambda_plots/%s_l_tilde_%s_%s.pdf", dataName.c_str(), polN[(i_pol-1)/2].c_str(), sqsName.c_str()));
      c->Clear();
    }
  }

  
  infile->Close();
}
