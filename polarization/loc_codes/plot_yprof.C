// code to plot cleaned up versions of y profiles
// to avoid repeating lengthy tree reading process
// reads TProfiles off of ROOT file

void plot_yprof()
{
  string type = "rho2delta0/";
  string dataName = "ups1";
  string sqsName = "7";
  TFile *infile = new TFile(Form("y_profs/plots_%s.root", dataName.c_str()));

  int y_bins = 4;
  double ylims[y_bins+1];
  for(int i = 0; i <= y_bins; i++)
    ylims[i] = i;
  
  TProfile ***y_prof = new TProfile**[4];
  for(int i = 0; i < 4; i++) {
    y_prof[i] = new TProfile*[y_bins];
    for(int j = 0; j < y_bins; j++) {
      y_prof[i][j] = (TProfile*)infile->Get(Form("%s%s_y_prof%d_y%d", dataName.c_str(), sqsName.c_str(), i, j));
    }
  }
  
  TCanvas *can = new TCanvas("","",700,700);
  for(int j = 0; j < 4; j++) {
    for(int i = 0; i < y_bins; i++) {
      y_prof[j][i]->GetYaxis()->SetRangeUser(0.0, 0.6);
      y_prof[j][i]->GetXaxis()->SetRangeUser(-2, 40);
      y_prof[j][i]->SetStats(0);

      y_prof[j][i]->GetYaxis()->SetTitle("#hat{y}/y");
      y_prof[j][i]->GetYaxis()->SetTitleOffset(1.2);
      y_prof[j][i]->GetXaxis()->SetTitle("#xi");
      y_prof[j][i]->Draw("error");
   
    }
  }

  string fr[4] = {"HX", "ggHX", "HX", "ggHX"};
  string pol[4] = {"tr", "tr", "lg", "lg"};
  
  // draw the plots
  y_prof[0][0]->SetTitle(Form("%s TeV #hat{y}/y transverse", sqsName.c_str()));
  y_prof[0][0]->Draw("error");
  for(int i = 0; i < y_bins; i++) {
    if(i!=0) y_prof[0][i]->Draw("same");
    y_prof[1][i]->Draw("same");
  }

  TLegend *legb = new TLegend(0.6, 0.1, 0.9, 0.4);
  legb->SetTextSize(0.03);
  for(int i = 0; i < y_bins; i++) {
    legb->AddEntry(y_prof[0][i], Form("%.0f < |y| < %.0f ppHX", ylims[i], ylims[i+1]), "pl");
    legb->AddEntry(y_prof[1][i], Form("%.0f < |y| < %.0f ggHX", ylims[i], ylims[i+1]), "pl");
  }
  legb->Draw();
  
  can->SaveAs(Form("y_profs/%s_y_trans_%s.pdf", dataName.c_str(), sqsName.c_str()));
  can->Clear();
	
  y_prof[2][0]->SetTitle(Form("%s TeV #hat{y}/y longitudinal", sqsName.c_str()));
  y_prof[2][0]->Draw("error");
  for(int i = 0; i < y_bins; i++) {
    if(i!=0) y_prof[2][i]->Draw("same");
    y_prof[3][i]->Draw("same");
  }
	
  legb->Draw();

  can->SaveAs(Form("y_profs/%s_y_long_%s.pdf", dataName.c_str(), sqsName.c_str()));
  can->Clear();
  
  infile->Close();
}
