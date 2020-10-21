// code to plot cleaned up versions of y profiles
// to avoid repeating lengthy tree reading process
// reads TProfiles off of ROOT file

void plot_yprof()
{
  string type = "rho2delta0/";
  string dataName = "jpsi";
  string sqsName = "7";
  TFile *infile = new TFile(Form("y_profs/plots_%s.root", dataName.c_str()));

  int xi_bins = 4;
  double xilims[xi_bins+1];
  for(int i = 0; i <= xi_bins; i++)
    xilims[i] = i*10.;
  xilims[0] = 1.;

  TProfile ***y_prof = new TProfile**[5];
  for(int i = 0; i < 5; i++) {
    y_prof[i] = new TProfile*[xi_bins];
    for(int j = 0; j < xi_bins; j++) {
      y_prof[i][j] = (TProfile*)infile->Get(Form("y_prof%d_xi%d", i, j));
      //     y_prof[i][j]->Sumw2();
    }
  }
  
  TCanvas *can = new TCanvas("","",700,700);
  for(int j = 0; j < 5; j++) {
    for(int i = 0; i < xi_bins; i++) {
      y_prof[j][i]->GetYaxis()->SetRangeUser(0, 0.3);
      y_prof[j][i]->SetStats(0);
    }
  }

  string fr[5] = {"HX", "ggHX", "HX", "ggHX", ""};
  string pol[5] = {"tr", "tr", "lg", "lg", "unp"};
  
  // draw the plots (one canvas for each pol)
  // same legend for all
  TLegend *legt = new TLegend(0.65, 0.7, 0.9, 0.9);
  legt->SetTextSize(0.03);
  for(int i = 0; i < xi_bins; i++) {
    legt->AddEntry(y_prof[0][i], Form("%.0f < #xi < %.0f", xilims[i], xilims[i+1]), "pl");
  }

  for(int i = 0; i < 5; i++) {
    y_prof[i][0]->SetTitle(Form("%s TeV #hat{y}/y %s %s", sqsName.c_str(), fr[i].c_str(), pol[i].c_str()));
    y_prof[i][0]->GetYaxis()->SetTitle("#hat{y}/y");
    y_prof[i][0]->GetYaxis()->SetTitleOffset(1.2);
    y_prof[i][0]->GetXaxis()->SetTitle("y");
    y_prof[i][0]->Draw("error");
    for(int j = 1; j < xi_bins; j++)
      y_prof[i][j]->Draw("same");
    legt->Draw();

    can->SaveAs(Form("y_profs/%s_y_%s_%s_%s.pdf", dataName.c_str(), fr[i].c_str(), pol[i].c_str(), sqsName.c_str()));
    can->Clear();

  }

  
  infile->Close();
}
