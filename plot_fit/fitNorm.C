// code that reads stored TGraph data and TH1 MC plots
// and fits normalizations, lumi NP and f_beta1

#import "pull_codes.C"
#import "dev_code.C"

// struct to store the pairs of MC vs data results
struct bin_coll {
  TH1F* MC_beta1;
  TH1F* MC_beta2;
  TGraphAsymmErrors* Data_graph;
  string exp;
  string state;
  int sqrts;
  double ylim[2];
};

vector <bin_coll> bin_c;

// get relative position on an axis (pi, pf)
double getPos(double pi, double pf, double mult, bool isLog) {
  if(isLog) return pow(10, log10(pi)+mult*(log10(pf)-log10(pi)));
  else return pi + mult*(pf-pi);
}

// function to parse a string into components separated by "deli"
vector< string > parseString( string line, string deli) {
  vector< string > out;
  string aux = line;
  size_t pos = 0;
  
  while( (pos = aux.find(deli)) != string::npos)
    {
      out.push_back( aux.substr( 0, pos ) );
      aux.erase( 0, pos + deli.length() );	    
    }
  out.push_back( aux );
  
  return out;
}

// chisquare function
double myFunction(double *par)
{
  const int n_norm = 2;
  double norm[n_norm];
  for(int i = 0; i < n_norm; i++)
    norm[i] = par[i];
  double LC7 = par[n_norm];
  double LC13 = par[n_norm+1];
  double LL7_j = par[n_norm+2];
  double LL7 = par[n_norm+3];
  double LL13 = par[n_norm+4];
  double f_beta1 = par[n_norm+5];
  
  double chisquare = 0;
  
  // cycle over all structs
  for(int i = 0; i <(int)bin_c.size(); i++) {
    // only fitting the xi plots for now
    TH1F* h1 = bin_c[i].MC_beta1;
    TH1F* h2 = bin_c[i].MC_beta2;
    TGraphAsymmErrors* g = bin_c[i].Data_graph;
    
    int nb = h1->GetNbinsX();
    
    double *g_y = g->GetY();
    double *g_e = g->GetEYhigh();
    
    // run over each bin in the histogram (same binning in both)
    for(int j = 0; j < nb; j++) 
      // if pT/M > 2, get chisquare contrib
      if(h1->GetBinLowEdge(j+1) > 2) {
	
	double mc_val = f_beta1*h1->GetBinContent(j+1)+(1.-f_beta1)*h2->GetBinContent(j+1);
	double data_val = g_y[j];
	double data_err = g_e[j];
	
	// state defines overall norm
	if(bin_c[i].state == "jpsi") mc_val *= norm[0];
	if(bin_c[i].state == "psi2") mc_val *= norm[1];
	//if(bin_c[i].state == "ups1") mc_val *= norm[2];
	//if(bin_c[i].state == "ups2") mc_val *= norm[3];
	//if(bin_c[i].state == "ups3") mc_val *= norm[4];
	
	// exp, sqrt(s) define lumi unc NP (in general)
	if(bin_c[i].exp == "CMS" && bin_c[i].sqrts == 7) {
	  data_val *= LC7;   data_err *= LC7; }
	else if(bin_c[i].exp == "CMS" && bin_c[i].sqrts == 13) {
	  data_val *= LC13;   data_err *= LC13; }
	else if(bin_c[i].exp == "LHCb" && bin_c[i].sqrts == 7 && bin_c[i].state == "jpsi") {
	  data_val *= LL7_j;   data_err *= LL7_j; }
	else if(bin_c[i].exp == "LHCb" && bin_c[i].sqrts == 7 && bin_c[i].state != "jpsi") {
	  data_val *= LL7;   data_err *= LL7; }
	else if(bin_c[i].exp == "LHCb" && bin_c[i].sqrts == 13) {
	  data_val *= LL13;   data_err *= LL13; }
	
	double dist = (data_val-mc_val)/data_err;
	chisquare += dist*dist;
      }
    
  
  }
  // constraints for the NP
  chisquare += pow((LC7-1.)/0.022,2);
  chisquare += pow((LC13-1.)/0.023,2);
  //chisquare += pow((LL7_j-1.)/0.035,2);
  //chisquare += pow((LL7-1.)/0.017,2);
  //chisquare += pow((LL13-1.)/0.039,2);
  
  return chisquare;
}

void minuitFunction(int& nDim, double* gout, double& result, double par[], int flg)
{
  result = myFunction(par);
}


// main fitting code
void fitNorm()
{
  const int n_state = 2;
  const int n_beta = 2;
  string sList[] = {"jpsi", "psi2"};
  int bList[] = {1, 2};
  double normFactor[n_beta][n_state];
  int n_plots[n_state];

  // open file first time to get normalizations and nr histos/graphs
  for(int is = 0; is < n_state; is++) {
    string state = sList[is];
    TFile *f1 = new TFile(Form("MC_rho3delta0/MC_vs_Data_%s.root", state.c_str()));
    
    TIter next(f1->GetListOfKeys());
    TKey *key;
    int n_tot = 0;

    while ((key = (TKey*)next())) {
      TClass *cl = gROOT->GetClass(key->GetClassName());
      if(cl->InheritsFrom("TTree")) {
	TTree *t = (TTree*)key->ReadObj();
	TLeaf *lf_nf1 = t->GetLeaf("normFactor_b1");
	lf_nf1->GetBranch()->GetEntry(0);
	normFactor[0][is] = lf_nf1->GetValue();
	TLeaf *lf_nf2 = t->GetLeaf("normFactor_b2");
	lf_nf2->GetBranch()->GetEntry(0);
	normFactor[1][is] = lf_nf2->GetValue();
      }
      else if (cl->InheritsFrom("TGraphAsymmErrors")) {
	TH1 *h = (TH1*)key->ReadObj();
	string name = h->GetName();
	if(name.find("p_{T}") != string::npos) continue;
	n_tot++;
      }
    }
    f1->Close();
    n_plots[is] = n_tot;
  }
  
  // create the array of struct, store data there
  int n_tot = 0;
  for(int is = 0; is < n_state; is++) n_tot += n_plots[is];
  bin_c.resize(n_tot);

  int ct = 0;
  for(int is = 0; is < n_state; is++) {
    string state = sList[is];
    TFile *f2 = new TFile(Form("MC_rho3delta0/MC_vs_Data_%s.root", state.c_str()));

    TIter n_ent(f2->GetListOfKeys());
    TKey *key;
    int test_v = 0;

    // structure is beta1 - beta2 - graph
    while ((key = (TKey*)n_ent())) {
      TClass *cl = gROOT->GetClass(key->GetClassName());
      // start by reading the two MC histos
      if (cl->InheritsFrom("TH1") && test_v < 2) {
	TH1F *h = (TH1F*)key->ReadObj();
	string name = h->GetName();
	if(name.find("p_{T}") != string::npos) {
	  continue;
	}
	// on the first (xi) histo, store the extra struct stuff
	if(test_v == 0) {
	  bin_c[ct].MC_beta1 = (TH1F*)key->ReadObj();
	  bin_c[ct].MC_beta1->SetDirectory(0);
	  bin_c[ct].MC_beta1->Scale(1./normFactor[test_v][is]);

	  vector <string> aux = parseString(name, " ");
	  bin_c[ct].ylim[0] = stof(aux[2]);
	  bin_c[ct].ylim[1] = stof(aux[6]);
	  if(stof(aux[2]) < 2) bin_c[ct].exp = "CMS";
	  else bin_c[ct].exp = "LHCb";
	  bin_c[ct].sqrts = stoi(aux[0]);
	  bin_c[ct].state = state;
	}
	// second histo just stores the histo itself
	else if(test_v == 1) {
	  bin_c[ct].MC_beta2 = (TH1F*)key->ReadObj();
	  bin_c[ct].MC_beta2->SetDirectory(0);
	  bin_c[ct].MC_beta2->Scale(1./normFactor[test_v][is]);
	}
 	
	test_v ++;
      }

      // just store the associated data graph
      else if(cl->InheritsFrom("TGraph") && test_v == 2) {
	bin_c[ct].Data_graph = (TGraphAsymmErrors*)key->ReadObj();
	ct++;
	test_v = 0;
      }
    }

    f2->Close();
  }

  // fitting process
  int n_pars = n_state + 5 + 1;
  
  TFitter *fit = new TFitter(n_pars);
  fit->SetFCN(minuitFunction);

  for(int is = 0; is < n_state; is++)
    fit->SetParameter(is, Form("%s_norm", sList[is].c_str()), normFactor[1][is], normFactor[1][is]/10., 0, 0);
  
  fit->SetParameter(n_state, "L_CMS_7", 1., 0.1, 0, 0);
  fit->SetParameter(n_state+1, "L_CMS_13", 1., 0.1, 0, 0);
  fit->SetParameter(n_state+2, "L_LHCb_7_j", 1., 0.1, 0, 0);
  fit->SetParameter(n_state+3, "L_LHCb_7", 1., 0.1, 0, 0);
  fit->SetParameter(n_state+4, "L_LHCb_13", 1., 0.1, 0, 0);
  //for(int inp = 0; inp < 5; inp++) fit->FixParameter(n_state+inp);

  fit->SetParameter(n_pars-1, "f_beta1", 0.1, 0.1, 0, 0);
  
  fit->ExecuteCommand("MIGRAD", 0, 0);

  //store results
  double pars[n_pars];
  for(int i = 0; i < n_pars; i++) pars[i] = fit->GetParameter(i);

  double minimum = myFunction(pars);

  int npts = 0;
  // get nr pts in fit
  for(int i = 0; i <(int)bin_c.size(); i++) {
    TH1F* h1 = bin_c[i].MC_beta1;
    int nb = h1->GetNbinsX();
    for(int j = 0; j < nb; j++)
      if(h1->GetBinLowEdge(j+1) > 2) {
	npts++;
      }
  }
  
  // print results
  cout << npts << " data pts" << endl;
  npts = npts-fit->GetNumberFreeParameters()+2; // ndf = npts - nr(non-NP free params)
  cout << npts << " ndf" << endl;

  cout << "chi^2 norm = " << minimum/(float)npts << endl;
  cout << "chi^2 prob = " << TMath::Prob(minimum, npts) << endl;

  // print fit results in a TeX file
  ofstream tex;
  tex.open("fitPar.tex");

  string stateTex[] = {"J/\\psi", "\\psi(2S)"};
  string NPTex[5] = {"\\text{CMS,7}", "\\text{CMS,13}", "\\text{LHCb,7},J/\\psi", "\\text{LHCb,7}", "\\text{LHCb,13}"};
  float lumi_unc[5] = {0.022, 0.023, 0.035, 0.017, 0.039};
  
  tex << "\\begin{table}[h!]" << endl;
  tex << "\\centering" << endl;
  tex << "\\begin{tabular}{c | c || c | c | c}" << endl;
  tex << "Parameter & Value $(\\times10^3)$ & Parameter & Value & Dev ($\\sigma$) \\\\" << endl;
  tex << "\\hline" << endl;
  for(int i_tex = 0; i_tex < 5; i_tex++) {
    if(i_tex < n_state) {
      int p_norm = ceil(-log10(fit->GetParError(i_tex)/1000.))+1;
      int p_NP = ceil(-log10(fit->GetParError(n_state+i_tex)))+1;

      float sig_dev = (pars[i_tex+n_state]-1.)/lumi_unc[i_tex];
      int p_sig = 1;
      if(abs(sig_dev) < 1) p_sig = ceil(-log10(abs(sig_dev)))+1;
    
      tex << "$N_{" << stateTex[i_tex] << "}$ & $" << setprecision(p_norm) << fixed << pars[i_tex]/1000. << "\\pm" << fit->GetParError(i_tex)/1000. << "$ & $L_{" << NPTex[i_tex] << "}$ & $" << setprecision(p_NP) << fixed << pars[i_tex+n_state] << "\\pm" << fit->GetParError(n_state+i_tex) << "$ & " << setprecision(p_sig) << fixed << sig_dev << " \\\\" << endl;
    }
    else {
      int p_NP = ceil(-log10(fit->GetParError(n_state+i_tex)))+1;

      float sig_dev = (pars[i_tex+n_state]-1.)/lumi_unc[i_tex];
      int p_sig = 1;
      if(abs(sig_dev) < 1) p_sig = ceil(-log10(abs(sig_dev)))+1;
      if(i_tex == n_state) tex << "\\hline " << endl;
      tex << " &  & $L_{" << NPTex[i_tex] << "}$ & $" << setprecision(p_NP) << fixed << pars[i_tex+n_state] << "\\pm" << fit->GetParError(n_state+i_tex) << "$ & " << setprecision(p_sig) << fixed << sig_dev << " \\\\" << endl;
    }
  }
  tex << "\\hline" << endl;
  int p_fbeta = ceil(-log10(fit->GetParError(n_pars-1)*100.))+1;
  tex << "& & $f_{\\beta1}$ & $(" << setprecision(p_fbeta) << fixed << pars[n_pars-1]*100. << "\\pm" << fit->GetParError(n_pars-1)*100 << ")$\\% & " << endl;
  tex << "\\end{tabular}" << endl;
  tex << "\\caption{Results of the normalization fit. " << Form("$\\chi^2/$ndf = %.0f / %d.}", minimum, npts) << endl;
  tex << "\\label{t:fit}" << endl;
  tex << "\\end{table}" << endl;
  tex.close();
    
  // store all results and plots in same file
  TFile *fout = new TFile("Fit_results.root", "recreate"); // will store updated plots
  
  TTree *fit_norms = new TTree("fit_norms", "fit norm factors");
  
  for(int is = 0; is < n_state; is++)
    TBranch *xi_nf = fit_norms->Branch(Form("%s_norm", sList[is].c_str()), &pars[is]);
  TBranch *C_7_unc = fit_norms->Branch("CMS_7_lumi_unc", &pars[n_state]);
  TBranch *C_13_unc = fit_norms->Branch("CMS_13_lumi_unc", &pars[n_state+1]);
  TBranch *L_7_j_unc = fit_norms->Branch("LHCb_7_jpsi_lumi_unc", &pars[n_state+2]);
  TBranch *L_7_unc = fit_norms->Branch("LHCb_7_lumi_unc", &pars[n_state+3]);
  TBranch *L_13_unc = fit_norms->Branch("LHCb_13_lumi_unc", &pars[n_state+4]);
  TBranch *f_b1 = fit_norms->Branch("f_beta1", &pars[n_pars-1]);
  
  fit_norms->Fill();

  // auxiliary funcs for plotting
  int xi_n = 1, y_n = 1;
  int sqsVals[] = {7,13};
  double y_min_ups[7] = {0, 0.6, 2, 2.5, 3, 3.5, 4};
  double y_max_ups[7] = {0.6, 1.2, 2.5, 3, 3.5, 4, 4.5};
  double y_min_psi[9] = {0, 0.3, 0.6, 0.9, 2, 2.5, 3, 3.5, 4};
  double y_max_psi[9] = {0.3, 0.6, 0.9, 1.2, 2.5, 3, 3.5, 4, 4.5};
  double y_yax_min_log, y_yax_max_log;

  TCanvas *can = new TCanvas("", "", 700, 700);
  can->SetLogy();

  // plotting: cycle over all n_state states to plot xi dists
  for(int is = 0; is < n_state; is++) {
    string state = sList[is];

    //for plotting need to define the limits etc
    if(state == "jpsi" || state == "psi2")
      xi_n = 9;
    else if(state == "ups1" || state == "ups2" || state == "ups3")
      xi_n = 7;

    double y_min[xi_n], y_max[xi_n]; // y bins of each xi dist
    double xi_min[xi_n], xi_max[xi_n]; // limits for plotting the xi dists (x-axis)
    double xi_yax_min[xi_n], xi_yax_max[xi_n]; // limits for plotting the xi dists (y-axis)

    for(int i_sqs = 0; i_sqs < 2; i_sqs++) {
      int sqsVal = sqsVals[i_sqs];

      // y limits are the same for both psi states in both sqrt(s)
      if(state == "jpsi" || state == "psi2") {
	for(int i = 0; i < xi_n; i++) {
	  y_min[i] = y_min_psi[i];
	  y_max[i] = y_max_psi[i]; 
	}
      }
      // y limits are the same for the 3 ups states in both sqrt(s)
      else if (state == "ups1" || state == "ups2" || state == "ups3") {
	for(int i = 0; i < xi_n; i++) {
	  y_min[i] = y_min_ups[i];
	  y_max[i] = y_max_ups[i]; 
	}
      }

      
      // one if condition for each state and sqrt(s)
      if(sqsVal == 7) {
	// filling in the j/psi 7 TeV values
	if (state == "jpsi") {  
	  for(int i = 0; i < xi_n; i++) {
	    xi_min[i] = -1;
	    xi_max[i] = i < 4 ? 35 : 6;
	    xi_yax_min[i] = i < 4 ? 1e-5 : 1;
	    xi_yax_max[i] = i < 4 ? 1e2 : 1e4;
	  }
	  y_yax_min_log = 1;
	  y_yax_max_log = 1e4;
	}
	
	// filling in the psi(2S) 7 TeV values
	else if (state == "psi2") {
	  for(int i = 0; i < xi_n; i++) {
	    xi_min[i] = 0;
	    xi_max[i] = i < 4 ? 25 : 5;
	    xi_yax_min[i] = i < 4 ? 1e-4 : 1e-1;
	    xi_yax_max[i] = i < 4 ? 1e2 : 1e4;
	  }
	  y_yax_min_log = 1e-1;
	  y_yax_max_log = 1e3;
	}
	
	else if (state == "ups1") {
	  for(int i = 0; i < xi_n; i++) {
	    xi_min[i] = -1;
	    xi_max[i] = i < 2 ? 14 : 5;
	    xi_yax_min[i] = i < 2 ? 1e-4 : 1e-2;
	    xi_yax_max[i] = i < 2 ? 1e2 : 1e3;
	  }
	  y_yax_min_log = 1e-1;
	  y_yax_max_log = 1e2;
	}

	else if (state == "ups2") {
	  for(int i = 0; i < xi_n; i++) {
	    xi_min[i] = -1;
	    xi_max[i] = i < 2 ? 14 : 5;
	    xi_yax_min[i] = i < 2 ? 1e-4 : 1e-2;
	    xi_yax_max[i] = i < 2 ? 1e2 : 1e3;
	  }
	  y_yax_min_log = 1e-1;
	  y_yax_max_log = 1e2;
	}
	
	else if (state == "ups3") {
	  for(int i = 0; i < xi_n; i++) {
	    xi_min[i] = -1;
	    xi_max[i] = i < 2 ? 14 : 5;
	    xi_yax_min[i] = i < 2 ? 1e-4 : 1e-2;
	    xi_yax_max[i] = i < 2 ? 1e1 : 1e2;
	  }
	  y_yax_min_log = 1e-1;
	  y_yax_max_log = 1e2;
	}
      }
      
      // same for 13 TeV
      else if(sqsVal == 13) {
	if (state == "jpsi") {
	  for(int i = 0; i < xi_n; i++) {
	    xi_min[i] = -1;
	    xi_max[i] = i < 4 ? 40 : 6;	
	    xi_yax_min[i] = i < 4 ? 1e-5 : 1;
	    xi_yax_max[i] = i < 4 ? 1e2 : 1e4;
	  }
	  y_yax_min_log = 1;
	  y_yax_max_log = 1e4;	
	}

	else if (state == "psi2") {
	  for(int i = 0; i < xi_n; i++) {
	    xi_min[i] = 0;
	    xi_max[i] = i < 4 ? 30 : 6;
	    xi_yax_min[i] = i < 4 ? 1e-4 : 1e-1;
	    xi_yax_max[i] = i < 4 ? 1e2 : 1e4;
	  }
	  y_yax_min_log = 5e-1;
	  y_yax_max_log = 5e3;
	}

	else if (state == "ups1") {
	  for(int i = 0; i < xi_n; i++) {
	    xi_min[i] = -1;
	    xi_max[i] = i < 2 ? 14 : 5;
	    xi_yax_min[i] = i < 2 ? 1e-4 : 1e-2;
	    xi_yax_max[i] = i < 2 ? 1e2 : 1e3;
	  }
	  y_yax_min_log = 1e-1;
	  y_yax_max_log = 1e3;
	}

	else if (state == "ups2") {
	  for(int i = 0; i < xi_n; i++) {
	    xi_min[i] = -1;
	    xi_max[i] = i < 2 ? 14 : 5;
	    xi_yax_min[i] = i < 2 ? 1e-4 : 1e-2;
	    xi_yax_max[i] = i < 2 ? 1e2 : 1e3;
	  }
	  y_yax_min_log = 1e-1;
	  y_yax_max_log = 1e3;
	}

	else if (state == "ups3") {
	  for(int i = 0; i < xi_n; i++) {
	    xi_min[i] = -1;
	    xi_max[i] = i < 2 ? 14 : 5;
	    xi_yax_min[i] = i < 2 ? 1e-4 : 1e-2;
	    xi_yax_max[i] = i < 2 ? 1e1 : 1e2;
	  }
	  y_yax_min_log = 1e-1;
	  y_yax_max_log = 1e3;
	}
      }

      ct = 1;
      // cycle over all structs
      int j_y = 0;
      for(int i = 0; i <(int)bin_c.size(); i++) {
	// only fitting (+plotting) the xi plots for now
	if(bin_c[i].state == state && bin_c[i].sqrts == sqsVal) {
	  
	  TH1F *fc = can->DrawFrame(xi_min[j_y], xi_yax_min[j_y], xi_max[j_y], xi_yax_max[j_y]);
	  fc->SetXTitle("#xi");
	  fc->SetYTitle("d#sigma / d#xidy (nb/GeV)");
	  fc->GetYaxis()->SetTitleOffset(1.3);
	  fc->SetTitle(Form("%d TeV %.1f < y < %.1f", sqsVal, y_min[j_y], y_max[j_y]));
	  can->Modified();
	  can->SetTitle("");
	  
	  //TGraphAsymmErrors* g = bin_c[i].Data_graph;
	
	  // apply norms and corrections to struct histos and graph
	  bin_c[i].MC_beta1->Scale(pars[is]);
	  bin_c[i].MC_beta2->Scale(pars[is]);

	  for(int i_bin = 0; i_bin < bin_c[i].Data_graph->GetN(); i_bin++) {
	    if(bin_c[i].exp == "CMS" && sqsVal == 7) 
	      bin_c[i].Data_graph->GetY()[i_bin] *= pars[n_state];
	    else if(bin_c[i].exp == "CMS" && sqsVal == 13)
	      bin_c[i].Data_graph->GetY()[i_bin] *= pars[n_state+1];
	    else if(bin_c[i].exp == "LHCb" && sqsVal == 7 && state == "jpsi")
	      bin_c[i].Data_graph->GetY()[i_bin] *= pars[n_state+2];
	    else if(bin_c[i].exp == "LHCb" && sqsVal == 7 && state != "jpsi")
	      bin_c[i].Data_graph->GetY()[i_bin] *= pars[n_state+3];
	    else if(bin_c[i].exp == "LHCb" && sqsVal == 13)
	      bin_c[i].Data_graph->GetY()[i_bin] *= pars[n_state+4];
	  }
	  
	  TH1F *h_comb = (TH1F*)bin_c[i].MC_beta1->Clone();
	  h_comb->Scale(pars[n_pars-1]);
	  h_comb->Add(bin_c[i].MC_beta2, 1.-pars[n_pars-1]);

	  h_comb->Draw("hist same");
	  bin_c[i].Data_graph->Draw("p same");

	  vector <string> aux = parseString(h_comb->GetName(), " ");
	  string newname = "";
	  for(int i_elm = 0; i_elm < aux.size()-3; i_elm++) {
	    newname += aux[i_elm];
	    if(i_elm < aux.size() - 4)
	      newname += " ";
	  }
	  h_comb->SetName(Form("Fitted %s %s", state.c_str(), newname.c_str()));
	  bin_c[i].MC_beta1->SetName(Form("%s %s", state.c_str(), bin_c[i].MC_beta1->GetName()));
	  bin_c[i].MC_beta2->SetName(Form("%s %s", state.c_str(), bin_c[i].MC_beta2->GetName()));
	  
	  bin_c[i].MC_beta1->Write();
	  bin_c[i].MC_beta2->Write();
	  //h_comb->Write();

	  TLine* xi_cut = new TLine(2, xi_yax_min[j_y], 2, xi_yax_max[j_y]);
	  xi_cut->SetLineStyle(kDashed);
	  xi_cut->SetLineColor(kBlack);
	  xi_cut->Draw("same");
	  
	  TLatex lc;
	  double xpos, ypos;
	  lc.SetTextSize(0.03);
	  xpos = getPos(xi_min[j_y], xi_max[j_y], 0.7, 0);
	  ypos = getPos(xi_yax_min[j_y], xi_yax_max[j_y], 0.8, 1);
	  lc.DrawLatex(xpos, ypos, Form("f_{#beta=1} = %.1f%%", pars[n_pars-1]*100.));
	  xpos = getPos(xi_min[j_y], xi_max[j_y], 0.65, 0);
	  ypos = getPos(xi_yax_min[j_y], xi_yax_max[j_y], 0.9, 1);
	  lc.DrawLatex(xpos, ypos, Form("#chi^{2}/ndf = %.0f / %d", minimum, npts));
	  
	  
	  can->SaveAs(Form("plots_%s/xi_%d_y%d.pdf", state.c_str(), bin_c[i].sqrts, ct));
	  can->Clear();

	  j_y++;
	  ct++;
	}
      }
    }
  }
      
  fout->Write();
  fout->Close();
  
  // running the pulls in summarized mode
  can->SetLogy(0);
  can->SetLogx();
  string sexp[] = {"CMS", "LHCb"};
  
  // frame for pulls
  TH1F *fc_pull = can->DrawFrame(1, -10, 40, 10);
  fc_pull->SetXTitle("#xi");
  fc_pull->SetYTitle("pulls");
  fc_pull->GetYaxis()->SetTitleOffset(1.3);
  fc_pull->SetTitle(Form("#xi pulls"));
  can->Modified();
  can->SetTitle("");
  
  // legend for the pulls
  double diff = 0.05*6.;
  //    TLegend *leg_xip = new TLegend(0.1, 0.9-diff, 0.3, 0.9);
  TLegend *leg_xip = new TLegend(0.65, 0.15, 0.9, 0.15+diff);
  leg_xip->SetTextSize(0.03);

  int i_plot = 0;  

  // the three cycles for all the pulls
  for(int i_sqs = 0; i_sqs < 2; i_sqs++) {
    int sqsVal = sqsVals[i_sqs];
    for(int is = 0; is < n_state; is++) {
      string state = sList[is];      
      for(int ie = 0; ie < 2; ie++) {
	string exp = sexp[ie];
	
	int nsqs = 0;
	for(int i = 0; i <(int)bin_c.size(); i++) {
	  // only fitting the xi plots for now
	  if(bin_c[i].sqrts == sqsVal && bin_c[i].exp == exp && bin_c[i].state == state) nsqs++;
	}
	
	if(nsqs == 0) continue;
	
	// get plots in format that imported code can use
	TGraphAsymmErrors** xi_g = new TGraphAsymmErrors*[nsqs];
	TH1F** xi_h = new TH1F*[nsqs];
	TGraphAsymmErrors** xi_p = new TGraphAsymmErrors*[nsqs];
	
	ct = 0;
	for(int i = 0; i <(int)bin_c.size(); i++) {
	  // only fitting the xi plots for now
	  if(bin_c[i].sqrts == sqsVal && bin_c[i].exp == exp && bin_c[i].state == state) {
	    xi_g[ct] = bin_c[i].Data_graph;
	    TH1F *h_comb = (TH1F*)bin_c[i].MC_beta1->Clone();
	    h_comb->Scale(pars[n_pars-1]);
	    h_comb->Add(bin_c[i].MC_beta2, 1.-pars[n_pars-1]);

	    xi_h[ct] = h_comb;
	    ct++;
	  }
	}
	
	xi_Pulls(xi_g, xi_h, xi_p, nsqs);

	int i_pull;
	if(exp == "LHCb") i_pull = 2;
	else if(state == "jpsi" || state == "psi2") i_pull = 2;
	else i_pull = 1;
	
	xi_p[i_pull]->SetMarkerColor(is+1);
	xi_p[i_pull]->SetMarkerStyle(i_sqs == 0 ? 25 : 20);
	xi_p[i_pull]->SetMarkerSize(i_sqs == 0 ? 1 : .75);
	xi_p[i_pull]->SetLineColor(is+1);
	xi_p[i_pull]->SetLineStyle(i_sqs == 0 ? kDashed : kSolid);
	if(is+1 == 5)
	  {
	    xi_p[i_pull]->SetLineColor(6);
	    xi_p[i_pull]->SetMarkerColor(6);
	  }
	xi_p[i_pull]->Draw("pl same");
	//xi_p[i_pull]->Draw("p same");
	if(i_plot%2 == 0 && i_sqs == 1) leg_xip->AddEntry(xi_p[i_pull], Form("%s", state.c_str()), "l");

	i_plot++;
      }
    }
  }
  
  TF1* line = new TF1("line", "0", -1, 40);
  line->SetLineStyle(kDashed);
  line->SetLineColor(kBlack);
  line->Draw("same");
  
  TLine* vert = new TLine(2, -10, 2, 10);
  vert->SetLineStyle(kDashed);
  vert->SetLineColor(kBlack);
  vert->Draw("same");
  
  leg_xip->Draw();
  
  can->SaveAs(Form("xi_pulls_full.pdf"));
  can->Clear();
  
  // running the deviation btw data, model (each and global)
  for(int is = 0; is < n_state; is++) {
    string state = sList[is];
    for(int i_sqs = 0; i_sqs < 2; i_sqs++) {
      int sqsVal = sqsVals[i_sqs];
      
      // frame for devs
      TH1F *fc_devs = can->DrawFrame(1, -100, 40, 100);
      fc_devs->SetXTitle("#xi");
      fc_devs->SetYTitle("(data-model)/model (%)");
      fc_devs->GetYaxis()->SetTitleOffset(1.3);
      fc_devs->SetTitle(Form("%d TeV %s deviation data/model", sqsVal, state.c_str()));
      can->Modified();
      can->SetTitle("");
    
      // legend for the pulls
      //    TLegend *leg_xid = new TLegend(0.1, 0.9-diff, 0.3, 0.9);
      TLegend *leg_xid = new TLegend(0.55, 0.15, 0.9, 0.15+diff);
      leg_xid->SetTextSize(0.03);
      
      i_plot = 0;  
      
      // the three cycles for all the devs
      for(int ie = 0; ie < 2; ie++) {
	string exp = sexp[ie];
	
	int nsqs = 0;
	for(int i = 0; i <(int)bin_c.size(); i++) {
	  // only fitting the xi plots for now
	  if(bin_c[i].sqrts == sqsVal && bin_c[i].exp == exp && bin_c[i].state == state) nsqs++;
	}
	
	if(nsqs == 0) continue;
	
	// get plots in format that imported code can use
	TGraphAsymmErrors** xi_g = new TGraphAsymmErrors*[nsqs];
	TH1F** xi_h = new TH1F*[nsqs];
	TGraphAsymmErrors** xi_d = new TGraphAsymmErrors*[nsqs];
	double ylim[nsqs][2];
	
	ct = 0;
	for(int i = 0; i <(int)bin_c.size(); i++) {
	  // only fitting the xi plots for now
	  if(bin_c[i].sqrts == sqsVal && bin_c[i].exp == exp && bin_c[i].state == state) {
	    xi_g[ct] = bin_c[i].Data_graph;
	    TH1F *h_comb = (TH1F*)bin_c[i].MC_beta1->Clone();
	    h_comb->Scale(pars[n_pars-1]);
	    h_comb->Add(bin_c[i].MC_beta2, 1.-pars[n_pars-1]);
	    xi_h[ct] = h_comb;
	    
	    ylim[ct][0] = bin_c[i].ylim[0];
	    ylim[ct][1] = bin_c[i].ylim[1];
	    ct++;
	  }
	}
	
	xi_Devs(xi_g, xi_h, xi_d, nsqs);
	
	for(int i_dev = 0; i_dev < nsqs; i_dev++) {	  
	  xi_d[i_dev]->SetMarkerColor(i_dev+1);
	  xi_d[i_dev]->SetMarkerStyle(20);
	  xi_d[i_dev]->SetMarkerSize(.75);
	  xi_d[i_dev]->SetLineColor(i_dev+1);
	  if (ie == 1) xi_d[i_dev]->SetLineStyle(kDashed);
	  if(i_dev+1 == 5)
	    {
	      xi_d[i_dev]->SetLineColor(6);
	      xi_d[i_dev]->SetMarkerColor(6);
	    }
	  xi_d[i_dev]->Draw("pl same");
	  //xi_d[i_dev]->Draw("p same");
	  leg_xid->AddEntry(xi_d[i_dev], Form("%s %.1f < y < %.1f", exp.c_str(), ylim[i_dev][0], ylim[i_dev][1]), "l");
	  
	  i_plot++;
	}
      }
    
      
      TF1* line_d = new TF1("line_d", "0", -1, 40);
      line_d->SetLineStyle(kDashed);
      line_d->SetLineColor(kBlack);
      line_d->Draw("same");
      
      TLine* vert_d = new TLine(2, -100, 2, 100);
      vert_d->SetLineStyle(kDashed);
      vert_d->SetLineColor(kBlack);
      vert_d->Draw("same");
      
      leg_xid->Draw();
      
      can->SaveAs(Form("xi_devs_%s_%d.pdf", state.c_str(), sqsVal));
      can->Clear();
    }
  }
  
  // now all states together
  // frame for devs
  TH1F *fc_devs = can->DrawFrame(1, -100, 40, 100);
  fc_devs->SetXTitle("#xi");
  fc_devs->SetYTitle("(data-model)/model (%)");
  fc_devs->GetYaxis()->SetTitleOffset(1.3);
  fc_devs->SetTitle(Form("deviation data/model"));
  can->Modified();
  can->SetTitle("");
  
  // legend for the devs
  //    TLegend *leg_xid = new TLegend(0.1, 0.9-diff, 0.3, 0.9);
  TLegend *leg_xid = new TLegend(0.65, 0.15, 0.9, 0.15+diff);
  leg_xid->SetTextSize(0.03);

  i_plot = 0;  

  // the three cycles for all the devs
  for(int i_sqs = 0; i_sqs < 2; i_sqs++) {
    int sqsVal = sqsVals[i_sqs];
    for(int is = 0; is < n_state; is++) {
      string state = sList[is];      
      for(int ie = 0; ie < 2; ie++) {
	string exp = sexp[ie];
	
	int nsqs = 0;
	for(int i = 0; i <(int)bin_c.size(); i++) {
	  // only fitting the xi plots for now
	  if(bin_c[i].sqrts == sqsVal && bin_c[i].exp == exp && bin_c[i].state == state) nsqs++;
	}
	
	if(nsqs == 0) continue;
	
	// get plots in format that imported code can use
	TGraphAsymmErrors** xi_g = new TGraphAsymmErrors*[nsqs];
	TH1F** xi_h = new TH1F*[nsqs];
	TGraphAsymmErrors** xi_d = new TGraphAsymmErrors*[nsqs];
	
	ct = 0;
	for(int i = 0; i <(int)bin_c.size(); i++) {
	  // only fitting the xi plots for now
	  if(bin_c[i].sqrts == sqsVal && bin_c[i].exp == exp && bin_c[i].state == state) {
	    xi_g[ct] = bin_c[i].Data_graph;
	    TH1F *h_comb = (TH1F*)bin_c[i].MC_beta1->Clone();
	    h_comb->Scale(pars[n_pars-1]);
	    h_comb->Add(bin_c[i].MC_beta2, 1.-pars[n_pars-1]);
	    xi_h[ct] = h_comb;
	    	    
	    ct++;
	  }
	}
	
	xi_Devs(xi_g, xi_h, xi_d, nsqs);

	int i_dev;
	if(exp == "LHCb") i_dev = 2;
	else if(state == "jpsi" || state == "psi2") i_dev = 2;
	else i_dev = 1;
	
	xi_d[i_dev]->SetMarkerColor(is+1);
	xi_d[i_dev]->SetMarkerStyle(i_sqs == 0 ? 25 : 20);
	xi_d[i_dev]->SetMarkerSize(i_sqs == 0 ? 1 : .75);
	xi_d[i_dev]->SetLineColor(is+1);
	xi_d[i_dev]->SetLineStyle(i_sqs == 0 ? kDashed : kSolid);
	if(is+1 == 5)
	  {
	    xi_d[i_dev]->SetLineColor(6);
	    xi_d[i_dev]->SetMarkerColor(6);
	  }
	xi_d[i_dev]->Draw("pl same");
	//xi_d[i_dev]->Draw("p same");
	if(i_plot%2 == 0 && i_sqs == 1) leg_xid->AddEntry(xi_d[i_dev], Form("%s", state.c_str()), "l");

	i_plot++;
      }
    }
  }
    
  TF1* line_d = new TF1("line_d", "0", -1, 40);
  line_d->SetLineStyle(kDashed);
  line_d->SetLineColor(kBlack);
  line_d->Draw("same");
  
  TLine* vert_d = new TLine(2, -100, 2, 100);
  vert_d->SetLineStyle(kDashed);
  vert_d->SetLineColor(kBlack);
  vert_d->Draw("same");
  
  leg_xid->Draw();
  
  can->SaveAs(Form("xi_devs_full.pdf"));
  can->Clear();
  
}


