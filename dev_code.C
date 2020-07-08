void xi_Devs(TGraphAsymmErrors** graph, TH1F** histo, TGraphAsymmErrors** devs, int nvals)
{
  for(int i = 0; i < nvals; i++)
    {
      int nb = histo[i]->GetNbinsX();
      if(nb == 1)
	continue;
      
      // get bin at which the comparison starts (pT/M > 2)
      int fb;
      for(int j = 0; j < nb; j++)
	if(histo[i]->GetBinContent(j+1) > 0) {
	  fb = j;
	  break;
	}
      
      int npts = nb-fb;
      double xvals[nb-fb], yvals[nb-fb], exlo[nb-fb], exhi[nb-fb], eylo[nb-fb], eyhi[nb-fb];
      double *g_x = graph[i]->GetX();
      double *g_y = graph[i]->GetY();
      double *g_eyl = graph[i]->GetEYlow();
      double *g_eyh = graph[i]->GetEYhigh();
      double *g_exl = graph[i]->GetEXlow();
      double *g_exh = graph[i]->GetEXhigh();
      
      for(int j = 0; j < nb-fb; j++)
	{
	  double mc_val = histo[i]->GetBinContent(j+fb+1);
	  double data_val = g_y[j+fb];

	  xvals[j] = g_x[j+fb];
	  exlo[j] = g_exl[j+fb];
	  exhi[j] = g_exh[j+fb];
	  yvals[j] = (data_val-mc_val)/mc_val*100.;
	  //if (yvals[j] > 60) cout << xvals[j] << " " << data_val << " " << mc_val << endl;
	  eylo[j] = g_eyl[j+fb]/mc_val*100.;
	  eyhi[j] = g_eyh[j+fb]/mc_val*100.;
	}
      devs[i] = new TGraphAsymmErrors(npts, xvals, yvals, exlo, exhi, eylo, eyhi); 
    }

}
