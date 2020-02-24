void xi_Pulls(TGraphAsymmErrors** graph, TH1F** histo, TGraphAsymmErrors** pull, int nvals)
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
      double xvals[nb-fb], yvals[nb-fb], exlo[nb-fb], exhi[nb-fb], ey[nb-fb];
      double *g_x = graph[i]->GetX();
      double *g_y = graph[i]->GetY();
      double *g_e = graph[i]->GetEYhigh();
      double *g_exl = graph[i]->GetEXlow();
      double *g_exh = graph[i]->GetEXhigh();
      
      for(int j = 0; j < nb-fb; j++)
	{
	  double mc_val = histo[i]->GetBinContent(j+fb+1);
	  double data_val = g_y[j+fb];
	  double data_err = g_e[j+fb];

	  xvals[j] = g_x[j+fb];
	  exlo[j] = g_exl[j+fb];
	  exhi[j] = g_exh[j+fb];
	  yvals[j] = (data_val-mc_val)/data_err;
	  ey[j] = 0;
	}
      pull[i] = new TGraphAsymmErrors(npts, xvals, yvals, exlo, exhi, ey, ey); 
    }
}


void y_Pulls(TGraphAsymmErrors** graph, TH1F** histo, TGraphAsymmErrors** pull, int nvals)
{
  for(int i = 0; i < nvals; i++)
    {
      int npts = graph[i]->GetN();
      double xvals[npts], yvals[npts], exlo[npts], exhi[npts], ey[npts];

      double *g_x = graph[i]->GetX();
      double *g_y = graph[i]->GetY();
      double *g_e = graph[i]->GetEYhigh();
      double *g_exl = graph[i]->GetEXlow();
      double *g_exh = graph[i]->GetEXhigh();

      TAxis *xaxis = histo[i]->GetXaxis();
      for(int j = 0; j < npts; j++)
	{
	  int bin = xaxis->FindBin(g_x[j]);
	  
	  double mc_val = histo[i]->GetBinContent(bin);
	  double data_val = g_y[j];
	  double data_err = g_e[j];

	  xvals[j] = g_x[j];
	  exlo[j] = g_exl[j];
	  exhi[j] = g_exh[j];
	  yvals[j] = (data_val-mc_val)/data_err;
	  ey[j] = 0;

	}
      pull[i] = new TGraphAsymmErrors(npts, xvals, yvals, exlo, exhi, ey, ey); 

    }
}
