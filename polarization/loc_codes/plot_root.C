// code to plot cleaned up versions of profiles
// to avoid repeating lengthy tree reading process
// reads TProfiles off of ROOT file

void plot_root()
{
  string type = "rho2delta0/";
  string dataName = "jpsi";
  string sqsName = "7";
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
      costh[i][j] = (TProfile*)infile->Get(Form("%s%s_ct%d_y%d", dataName.c_str(), sqsName.c_str(), i, j));
      phi[i][j] = (TProfile*)infile->Get(Form("%s%s_phi%d_y%d", dataName.c_str(), sqsName.c_str(), i, j));
      phith[i][j] = (TProfile*)infile->Get(Form("%s%s_pt%d_y%d", dataName.c_str(), sqsName.c_str(), i, j));
    }
  }

  TCanvas *can = new TCanvas("","",700,700);
  for(int j = 0; j < 4; j++) {
    for(int i = 0; i < y_bins; i++) {
      costh[j][i]->GetYaxis()->SetRangeUser(0.1,0.5);
      phi[j][i]->GetYaxis()->SetRangeUser(-0.6,0.6);
      phith[j][i]->GetYaxis()->SetRangeUser(-0.5,0.5);

      costh[j][i]->SetStats(0);
      phi[j][i]->SetStats(0);
      phith[j][i]->SetStats(0);
    }
  }
  
  // draw the costh plots
  costh[0][0]->SetTitle(Form("%s TeV cos^{2}#theta_{HX} transverse", sqsName.c_str()));
  costh[0][0]->Draw("error");
  for(int i = 0; i < y_bins; i++) {
    if(i!=0) costh[0][i]->Draw("same");
    costh[1][i]->Draw("same");
  }
  TLine *ct1 = new TLine(1, 0.2, 40, 0.2);
  ct1->SetLineStyle(kDotted);
  ct1->Draw();
  TLine *ct2 = new TLine(1, 0.4, 40, 0.4);
  ct2->SetLineStyle(kDotted);
  ct2->Draw();

  TLegend *legt = new TLegend(0.6, 0.6, 0.9, 0.9);
  legt->SetTextSize(0.03);
  TLegend *legb = new TLegend(0.6, 0.1, 0.9, 0.4);
  legb->SetTextSize(0.03);
  for(int i = 0; i < y_bins; i++) {
    legt->AddEntry(costh[0][i], Form("%.0f < |y| < %.0f ppHX", ylims[i], ylims[i+1]), "pl");
    legt->AddEntry(costh[1][i], Form("%.0f < |y| < %.0f ggHX", ylims[i], ylims[i+1]), "pl");
    legb->AddEntry(costh[0][i], Form("%.0f < |y| < %.0f ppHX", ylims[i], ylims[i+1]), "pl");
    legb->AddEntry(costh[1][i], Form("%.0f < |y| < %.0f ggHX", ylims[i], ylims[i+1]), "pl");
  }
  legb->Draw();
  
  can->SaveAs(Form("ang_profs/%s_cos2th_trans_%s.pdf", dataName.c_str(), sqsName.c_str()));
  can->Clear();
	
  costh[2][0]->SetTitle(Form("%s TeV cos^{2}#theta_{HX} longitudinal", sqsName.c_str()));
  costh[2][0]->Draw("error");
  for(int i = 0; i < y_bins; i++) {
    if(i!=0) costh[2][i]->Draw("same");
    costh[3][i]->Draw("same");
  }
  ct1->Draw();
  ct2->Draw();
	
  legt->Draw();

  can->SaveAs(Form("ang_profs/%s_cos2th_long_%s.pdf", dataName.c_str(), sqsName.c_str()));
  can->Clear();
	
  // draw the phi plots
  phi[0][0]->SetTitle(Form("%s TeV cos2#phi_{HX} transverse", sqsName.c_str()));
  phi[0][0]->Draw("error");
  for(int i = 0; i < y_bins; i++) {
    if(i!=0) phi[0][i]->Draw("same");
    phi[1][i]->Draw("same");
  }
	
  TLine *ph1 = new TLine(1, -0.5, 40, -0.5);
  ph1->SetLineStyle(kDotted);
  ph1->Draw();
  TLine *ph2 = new TLine(1, 0.5, 40, 0.5);
  ph2->SetLineStyle(kDotted);
  ph2->Draw();

  legt->Draw();
  
  can->SaveAs(Form("ang_profs/%s_cos2phi_trans_%s.pdf", dataName.c_str(), sqsName.c_str()));
  can->Clear();
	
  phi[2][0]->SetTitle(Form("%s TeV cos2#phi_{HX} longitudinal", sqsName.c_str()));
  phi[2][0]->Draw("error");
  for(int i = 0; i < y_bins; i++) {
    if(i!=0) phi[2][i]->Draw("same");
    phi[3][i]->Draw("same");
  }
	
  ph1->Draw();
  ph2->Draw();

  legt->Draw();
  
  can->SaveAs(Form("ang_profs/%s_cos2phi_long_%s.pdf", dataName.c_str(), sqsName.c_str()));
  can->Clear();
	
  // draw the phith plots
  phith[0][0]->SetTitle(Form("%s TeV sin2#theta_{HX}cos#phi_{HX} transverse", sqsName.c_str()));
  phith[0][0]->Draw("error");
  for(int i = 0; i < y_bins; i++) {
    if(i!=0) phith[0][i]->Draw("same");
    phith[1][i]->Draw("same");
  }

  TLine *pt1 = new TLine(1, -0.4, 40, -0.4);
  pt1->SetLineStyle(kDotted);
  pt1->Draw();
  TLine *pt2 = new TLine(1, 0.4, 40, 0.4);
  pt2->SetLineStyle(kDotted);
  pt2->Draw();

  legt->Draw();

  can->SaveAs(Form("ang_profs/%s_phitheta_trans_%s.pdf", dataName.c_str(), sqsName.c_str()));
  can->Clear();
	
  phith[2][0]->SetTitle(Form("%s TeV sin2#theta_{HX}cos#phi_{HX} longitudinal", sqsName.c_str()));
  phith[2][0]->Draw("error");
  for(int i = 0; i < y_bins; i++) {
    if(i!=0) phith[2][i]->Draw("same");
    phith[3][i]->Draw("same");
  }

  pt1->Draw();
  pt2->Draw();


  legt->Draw();

  can->SaveAs(Form("ang_profs/%s_phitheta_long_%s.pdf", dataName.c_str(), sqsName.c_str()));
  can->Clear();
  
infile->Close();
}
