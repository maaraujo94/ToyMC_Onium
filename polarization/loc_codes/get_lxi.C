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

void get_lxi()
{ 
  string dataName = "jpsi", sqsName = "7";

  TFile *infile = new TFile(Form("ang_profs/plots_%s.root", dataName.c_str()));

  int y_bins = 8;
  double ylims[y_bins+1];
  for(int i = 0; i <= y_bins; i++)
    ylims[i] = 0.5*i;

  int xi_bins = 5;
  
  TProfile ***costh = new TProfile**[4];
  TProfile ***phi = new TProfile**[4];
  TProfile ***phith = new TProfile**[4];
  for(int i = 0; i < 4; i++) {
    costh[i] = new TProfile*[y_bins];
    phi[i] = new TProfile*[y_bins];
    phith[i] = new TProfile*[y_bins];
    for(int j = 0; j < y_bins; j++) {
      costh[i][j] = (TProfile*)infile->Get(Form("%s%s_ct%d_y%d", dataName.c_str(), sqsName.c_str(), i, j));
      phi[i][j] = (TProfile*)infile->Get(Form("%s%s_phi%d_y%d", dataName.c_str(), sqsName.c_str(), i, j));
      phith[i][j] = (TProfile*)infile->Get(Form("%s%s_pt%d_y%d", dataName.c_str(), sqsName.c_str(), i, j));
    }
  }


  TGraphErrors ***g_lth = new TGraphErrors**[4];
  TGraphErrors ***g_lp = new TGraphErrors**[4];
  TGraphErrors ***g_ltp = new TGraphErrors**[4];
  TGraphErrors ***g_ltil = new TGraphErrors**[4];

  for(int i_pol = 0; i_pol < 4; i_pol++) {
    g_lth[i_pol] = new TGraphErrors*[xi_bins];
    g_lp[i_pol] = new TGraphErrors*[xi_bins];
    g_ltp[i_pol] = new TGraphErrors*[xi_bins];
    g_ltil[i_pol] = new TGraphErrors*[xi_bins];
    
    for(int i_xi = 0; i_xi < xi_bins; i_xi++) {
      double lth[y_bins], elth[y_bins], lp[y_bins], elp[y_bins], ltp[y_bins], eltp[y_bins], ltil[y_bins], eltil[y_bins];
      double pt[y_bins], ept[y_bins];
      
      for(int i_y = 0; i_y < y_bins; i_y++) {
	pt[i_y] = 0.5*(ylims[i_y+1]+ylims[i_y]);
	ept[i_y] = 0.5*(ylims[i_y+1]-ylims[i_y]);

	double cos2th = costh[i_pol][i_y]->GetBinContent(i_xi+1);
	double ecth = costh[i_pol][i_y]->GetBinError(i_xi+1);
	lth[i_y] = lambda_theta(cos2th);
	elth[i_y] = e_lth(cos2th, ecth);

	double cos2phi = phi[i_pol][i_y]->GetBinContent(i_xi+1);
	double ecp = phi[i_pol][i_y]->GetBinError(i_xi+1);
	lp[i_y] = lambda_phi(cos2phi, lth[i_y]);
	elp[i_y] = e_lp(cos2phi, ecp, lth[i_y], elth[i_y]);

      	double sc = phith[i_pol][i_y]->GetBinContent(i_xi+1);
	double esc = phith[i_pol][i_y]->GetBinError(i_xi+1);
	ltp[i_y] = lambda_thetaphi(sc, lth[i_y]);
	eltp[i_y] = e_ltp(sc, esc, lth[i_y], elth[i_y]);

	ltil[i_y] = lambda_tilde(lth[i_y], lp[i_y]);
	eltil[i_y] = 0.;
	
      }
      
      g_lth[i_pol][i_xi] = new TGraphErrors(y_bins, pt, lth, ept, elth);
      g_lth[i_pol][i_xi]->SetLineColor(i_xi+1);
      g_lth[i_pol][i_xi]->SetMarkerColor(i_xi+1);
      if(i_pol%2 == 1) g_lth[i_pol][i_xi]->SetLineStyle(kDashed);
      g_lth[i_pol][i_xi]->SetMarkerStyle(i_pol%2 == 0 ? 20 : 24);
      //if(i_pol%2 == 0)
      g_lth[i_pol][i_xi]->SetMarkerSize(.75);
      
      g_lp[i_pol][i_xi] = new TGraphErrors(y_bins, pt, lp, ept, elp);
      g_lp[i_pol][i_xi]->SetLineColor(i_xi+1);
      g_lp[i_pol][i_xi]->SetMarkerColor(i_xi+1);
      if(i_pol%2 == 1) g_lp[i_pol][i_xi]->SetLineStyle(kDashed);
      g_lp[i_pol][i_xi]->SetMarkerStyle(i_pol%2 == 0 ? 20 : 24);
      //if(i_pol%2 == 0)
      g_lp[i_pol][i_xi]->SetMarkerSize(.75);

      g_ltp[i_pol][i_xi] = new TGraphErrors(y_bins, pt, ltp, ept, eltp);
      g_ltp[i_pol][i_xi]->SetLineColor(i_xi+1);
      g_ltp[i_pol][i_xi]->SetMarkerColor(i_xi+1);
      if(i_pol%2 == 1) g_ltp[i_pol][i_xi]->SetLineStyle(kDashed);
      g_ltp[i_pol][i_xi]->SetMarkerStyle(i_pol%2 == 0 ? 20 : 24);
      //if(i_pol%2 == 0)
      g_ltp[i_pol][i_xi]->SetMarkerSize(.75);

      g_ltil[i_pol][i_xi] = new TGraphErrors(y_bins, pt, ltil, ept, eltil);
      g_ltil[i_pol][i_xi]->SetLineColor(i_xi+1);
      g_ltil[i_pol][i_xi]->SetMarkerColor(i_xi+1);
      if(i_pol%2 == 1) g_ltil[i_pol][i_xi]->SetLineStyle(kDashed);
      g_ltil[i_pol][i_xi]->SetMarkerStyle(i_pol%2 == 0 ? 20 : 24);
      //if(i_pol%2 == 0)
      g_ltil[i_pol][i_xi]->SetMarkerSize(.75);

    }
  }

  TCanvas *c = new TCanvas("","",700,700);
  string polN[2] = {"trans", "long"};
  string fullN[2] = {"transverse", "longitudinal"};

  double xilims[xi_bins+1];
  for(int i = 0; i <= xi_bins; i++)
    xilims[i] = i+1.;
  
  TLegend *leg = new TLegend(0.6, 0.75, 0.9, 0.9);
  leg->SetTextSize(0.03);
  for(int i = 0; i < xi_bins; i++) {
    leg->AddEntry(g_lth[1][i], Form("%.0f < #xi < %.0f ggHX", xilims[i], xilims[i+1]), "pl");
  }
  
  double thLim[4] = {0, 1.5, -1.5, 0};
  for(int i_pol = 0; i_pol < 4; i_pol++) {
    if(i_pol%2 == 0) {
      TH1F *fc = c->DrawFrame(0, thLim[i_pol], 4, thLim[i_pol+1]);
      fc->SetXTitle("|y|");
      fc->SetYTitle("#lambda_{#theta}");
      fc->GetYaxis()->SetTitleOffset(1.3);
      fc->SetTitle(Form("%s TeV #lambda_{#theta} %s", sqsName.c_str(), fullN[i_pol/2].c_str()));
    }


    if(i_pol%2 == 1) {
    for(int i_xi = 0; i_xi < xi_bins; i_xi++) {
      g_lth[i_pol][i_xi]->Draw("psame");
    }
      leg->Draw();
      c->SaveAs(Form("lambda_plots/%s_l_theta_y_%s_%s.pdf", dataName.c_str(), polN[(i_pol-1)/2].c_str(), sqsName.c_str()));
      c->Clear();
    }
  }

  for(int i_pol = 0; i_pol < 4; i_pol++) {
    if(i_pol%2 == 0) {
      TH1F *fc = c->DrawFrame(0, -1, 4, 1);
      fc->SetXTitle("|y|");
      fc->SetYTitle("#lambda_{#phi}");
      fc->GetYaxis()->SetTitleOffset(1.3);
      fc->SetTitle(Form("%s TeV #lambda_{#phi} %s", sqsName.c_str(), fullN[i_pol/2].c_str()));
    }
    
    for(int i_xi = 0; i_xi < xi_bins; i_xi++) {
      g_lp[i_pol][i_xi]->Draw("psame");
    }

    if(i_pol%2 == 1) {
      leg->Draw();
      c->SaveAs(Form("lambda_plots/%s_l_phi_y_%s_%s.pdf", dataName.c_str(), polN[(i_pol-1)/2].c_str(), sqsName.c_str()));
      c->Clear();
    }
  }

  for(int i_pol = 0; i_pol < 4; i_pol++) {
    if(i_pol%2 == 0) {
      TH1F *fc = c->DrawFrame(0, -1, 4, 1);
      fc->SetXTitle("#xi");
      fc->SetYTitle("#lambda_{#theta#phi}");
      fc->GetYaxis()->SetTitleOffset(1.3);
      fc->SetTitle(Form("%s TeV #lambda_{#theta#phi} %s", sqsName.c_str(), fullN[i_pol/2].c_str()));
    }
    
    for(int i_xi = 0; i_xi < xi_bins; i_xi++) {
      g_ltp[i_pol][i_xi]->Draw("psame");
    }

    if(i_pol%2 == 1) {
      leg->Draw();
      c->SaveAs(Form("lambda_plots/%s_l_thetaphi_y_%s_%s.pdf", dataName.c_str(), polN[(i_pol-1)/2].c_str(), sqsName.c_str()));
      c->Clear();
    }
  }
    
  for(int i_pol = 0; i_pol < 4; i_pol++) {
    if(i_pol%2 == 0) {
      TH1F *fc = c->DrawFrame(0, thLim[i_pol], 4, thLim[i_pol+1]);
      //TH1F *fc = c->DrawFrame(0, -1, 40, 1);
      fc->SetXTitle("#xi");
      fc->SetYTitle("#tilde#lambda");
      fc->GetYaxis()->SetTitleOffset(1.3);
      fc->SetTitle(Form("%s TeV #tilde{#lambda} %s", sqsName.c_str(), fullN[i_pol/2].c_str()));
    }
    
    for(int i_xi = 0; i_xi < xi_bins; i_xi++) {
      g_ltil[i_pol][i_xi]->Draw("psame");
    }

    if(i_pol%2 == 1) {
      leg->Draw();
      c->SaveAs(Form("lambda_plots/%s_l_tilde_y_%s_%s.pdf", dataName.c_str(), polN[(i_pol-1)/2].c_str(), sqsName.c_str()));
      c->Clear();
    }
  }

  
  infile->Close();
}
