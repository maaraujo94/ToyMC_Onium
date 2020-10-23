void get_axang()
{ 
  string dataName = "jpsi", sqsName = "7";

  TFile *infile = new TFile(Form("axAng_plots/plots_%s.root", dataName.c_str()));

  int y_bins = 8;
  double ylims[y_bins+1];
  for(int i = 0; i <= y_bins; i++)
    ylims[i] = 0.5*i;

  int xi_bins = 5;
  
  TProfile ***costh = new TProfile**[4];
  for(int i = 0; i < 4; i++) {
    costh[i] = new TProfile*[xi_bins];
    for(int j = 0; j < xi_bins; j++) {
      costh[i][j] = (TProfile*)infile->Get(Form("%s%s_ang_prof%d_y%d", dataName.c_str(), sqsName.c_str(), i, j));
    }
  }

  TGraphErrors ***g_lth = new TGraphErrors**[4];

  for(int i_pol = 0; i_pol < 4; i_pol++) {
    g_lth[i_pol] = new TGraphErrors*[xi_bins];
    
    for(int i_xi = 0; i_xi < xi_bins; i_xi++) {
      double lth[y_bins], elth[y_bins];
      double pt[y_bins], ept[y_bins];
      
      for(int i_y = 0; i_y < y_bins; i_y++) {
	pt[i_y] = 0.5*(ylims[i_y+1]+ylims[i_y]);
	ept[i_y] = 0.5*(ylims[i_y+1]-ylims[i_y]);

	double cos2th = costh[i_pol][i_xi]->GetBinContent(i_y+1);
	double ecth = costh[i_pol][i_xi]->GetBinError(i_y+1);
	lth[i_y] = cos2th;
	elth[i_y] = ecth;	
      }
      
      g_lth[i_pol][i_xi] = new TGraphErrors(y_bins, pt, lth, ept, elth);
      g_lth[i_pol][i_xi]->SetLineColor(i_xi+1);
      g_lth[i_pol][i_xi]->SetMarkerColor(i_xi+1);
      if(i_pol%2 == 1) g_lth[i_pol][i_xi]->SetLineStyle(kDashed);
      g_lth[i_pol][i_xi]->SetMarkerStyle(i_pol%2 == 0 ? 20 : 24);
      //if(i_pol%2 == 0)
      g_lth[i_pol][i_xi]->SetMarkerSize(.75);
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
  
  for(int i_pol = 0; i_pol < 4; i_pol++) {
    if(i_pol%2 == 0) {
      TH1F *fc = c->DrawFrame(0, 0,  4, 0.6);
      fc->SetXTitle("|y|");
      fc->SetYTitle("angle");
      fc->GetYaxis()->SetTitleOffset(1.3);
      fc->SetTitle(Form("%s TeV HX-ggHX angle %s", sqsName.c_str(), fullN[i_pol/2].c_str()));
    }


    if(i_pol%2 == 1) {
    for(int i_xi = 0; i_xi < xi_bins; i_xi++) {
      g_lth[i_pol][i_xi]->Draw("psame");
    }
      leg->Draw();
      c->SaveAs(Form("axAng_plots/%s_angle_%s_%s.pdf", dataName.c_str(), polN[(i_pol-1)/2].c_str(), sqsName.c_str()));
      c->Clear();
    }
  }

  
  infile->Close();
}
