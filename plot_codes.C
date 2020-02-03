/* auxiliary file with all the data-reading codes
   - xi-ordered for J/psi, Psi(2S), Ups(1S) @ 7 TeV
   - y-ordered for same
*/

// xi codes
void xi_jpsi_7(TGraphAsymmErrors **g, double norm)
{
  double mass = 3.097;

  ifstream file;
  string data;
  int n_pts;
  double nums[12];
  const int n_y = 6; 

  string filelist[n_y] = {"data/CMS_jpsi_7_cs_y1.txt",
			  "data/LHCb_jpsi_7_cs_y1.txt",
			  "data/LHCb_jpsi_7_cs_y2.txt",
			  "data/LHCb_jpsi_7_cs_y3.txt",
			  "data/LHCb_jpsi_7_cs_y4.txt",
			  "data/LHCb_jpsi_7_cs_y5.txt"};
  for(int k = 0; k < n_y; k++) {
    file.open(filelist[k].c_str());
    n_pts = 0;
    while(getline(file,data)) n_pts+=1;
      
    double pt[n_pts], dpt_lo[n_pts], dpt_hi[n_pts], sig[n_pts], dsig[n_pts];
      
    file.clear();
    file.seekg(0);
    for(int j = 0; j < n_pts; j++) {
      for(int i = 0; i < 12; i++)
	file >> nums[i];
      
      pt[j] = nums[2]/mass;
      dpt_lo[j] = (nums[2]-nums[3])/mass;
      dpt_hi[j] = (nums[4]-nums[2])/mass;
	
      double signorm = k < 1 ? 5.961e-2*1e3 : 1;
      sig[j] = nums[5]*mass/(signorm);
      dsig[j] = sqrt(nums[6]*nums[6]+nums[8]*nums[8])*mass/(signorm);
    }
    file.close();
    
    // normalize to first point
    double aux;
    if(k < 1) aux = sig[0]/norm;
    for(int i = 0; i < n_pts; i++) {
      sig[i] /= aux;
      dsig[i] /= aux;
    }
      
    g[k] = new TGraphAsymmErrors(n_pts, pt, sig, dpt_lo, dpt_hi, dsig, dsig);
    
  }
  
}

// here there will be a TGraph that is just empty:
// check how to do this!
void xi_psi2_7(TGraphAsymmErrors **g, double norm)
{
  double mass = 3.686;
  
  ifstream file;
  string data;
  int n_pts;
  double nums[12];
  const int n_y = 6;

  string filelist[n_y] = {"",
			  "data/LHCb_psi2_7_cs_y1.txt",
			  "data/LHCb_psi2_7_cs_y2.txt",
			  "data/LHCb_psi2_7_cs_y3.txt",
			  "data/LHCb_psi2_7_cs_y4.txt",
			  "data/LHCb_psi2_7_cs_y5.txt"};
  for(int k = 0; k < n_y; k++) {
    file.open(filelist[k].c_str());
    n_pts = 0;
    while(getline(file,data)) n_pts+=1;
      
    double pt[n_pts], dpt_lo[n_pts], dpt_hi[n_pts], sig[n_pts], dsig[n_pts];
    
    file.clear();
    file.seekg(0);
    for(int j = 0; j < n_pts; j++) {
      for(int i = 0; i < 12; i++)
	file >> nums[i];
      
      pt[j] = nums[2]/mass;
      dpt_lo[j] = (nums[2]-nums[3])/mass;
      dpt_hi[j] = (nums[4]-nums[2])/mass;
      
      double signorm = 1.;
      sig[j] = nums[5]*mass/(signorm);
      dsig[j] = sqrt(nums[6]*nums[6]+nums[8]*nums[8])*mass/(signorm);
    }
    file.close();
    
    // normalize to first point
    double aux;
    if(k == 1) aux = sig[0]/norm;
    for(int i = 0; i < n_pts; i++) {
      sig[i] /= aux;
      dsig[i] /= aux;
    }

    g[k] = new TGraphAsymmErrors(n_pts, pt, sig, dpt_lo, dpt_hi, dsig, dsig);
  }
  
}

void xi_ups1_7(TGraphAsymmErrors **g, double norm)
{
  double mass = 9.46;
  
  ifstream file;
  string data;
  int n_pts;
  double nums[12];
  const int n_y = 6; 

  
  string filelist[n_y] = {"data/CMS_ups1_7_cs_y1.txt",
			  "data/LHCb_ups1_7_cs_y1.txt",
			  "data/LHCb_ups1_7_cs_y2.txt",
			  "data/LHCb_ups1_7_cs_y3.txt",
			  "data/LHCb_ups1_7_cs_y4.txt",
			  "data/LHCb_ups1_7_cs_y5.txt"};

  for(int k = 0; k < n_y; k++) {
    file.open(filelist[k].c_str());
    n_pts = 0;
    while(getline(file,data)) n_pts+=1;
    
    double pt[n_pts], dpt_lo[n_pts], dpt_hi[n_pts], sig[n_pts], dsig[n_pts];
    
    file.clear();
    file.seekg(0);
    for(int j = 0; j < n_pts; j++) {
      for(int i = 0; i < 12; i++)
	file >> nums[i];
      
      pt[j] = nums[2]/mass;
      dpt_lo[j] = (nums[2]-nums[3])/mass;
      dpt_hi[j] = (nums[4]-nums[2])/mass;
      
      double signorm = k < 1 ? 2.4*2.48e-2*1e6 : 2.48e-2*1e3;
      sig[j] = nums[5]*mass/(signorm);
      dsig[j] = sqrt(nums[6]*nums[6]+nums[8]*nums[8])*mass/(signorm);
    }
    file.close();
    
    // normalize to first point
    double aux;
    if(k < 1) aux = sig[0]/norm;
    for(int i = 0; i < n_pts; i++) {
      sig[i] /= aux;
      dsig[i] /= aux;
    }
    
    g[k] = new TGraphAsymmErrors(n_pts, pt, sig, dpt_lo, dpt_hi, dsig, dsig);
    
  }
  
}

// y codes
void y_jpsi_7(TGraphAsymmErrors **g, double norm, int y_n)
{
  double mass = 3.097;
  
  int n_val[y_n];
  for (int i = 0; i < y_n; i++)
    n_val[i] = i == 6 ? 4 : i == 3 ? 6 : 5;

  ifstream file;
  string data;
  int n_pts[6] = {31, 12, 12, 12, 11, 9};
  double pt_lo[6][n_pts[0]], pt_hi[6][n_pts[0]], y_v[6][n_pts[0]], dy[6][n_pts[0]], sig[6][n_pts[0]], dsig[6][n_pts[0]];
  double nums[12];
  
  string filelist[6] = {"data/CMS_jpsi_7_cs_y1.txt",
			"data/LHCb_jpsi_7_cs_y1.txt",
			"data/LHCb_jpsi_7_cs_y2.txt",
			"data/LHCb_jpsi_7_cs_y3.txt",
			"data/LHCb_jpsi_7_cs_y4.txt",
			"data/LHCb_jpsi_7_cs_y5.txt"};
  for(int k = 0; k < 6; k++) {
    file.open(filelist[k].c_str());
	
    for(int j = 0; j < n_pts[k]; j++) {
      for(int i = 0; i < 12; i++)
	file >> nums[i];
	  
      y_v[k][j] = 0.5 * (nums[1] + nums[0]);
      dy[k][j] = 0.5 * (nums[1] - nums[0]);
	  
      pt_lo[k][j] = nums[3];
      pt_hi[k][j] = nums[4];
	  
      double signorm = k < 1 ? 5.961e-2*1e3 : 1;
      sig[k][j] = nums[5]*mass/(signorm);
      dsig[k][j] = sqrt(nums[6]*nums[6]+nums[8]*nums[8])*mass/(signorm);
    }
    file.close();
	
    // normalize to first point
    double aux;
    if(k < 1) aux = sig[0][0]/norm;
    for(int i = 0; i < n_pts[k]; i++) {
      sig[k][i] /= aux;
      dsig[k][i] /= aux;
    }
    
  }
  
  for(int i = 0; i < y_n; i++) {
    double y_av[n_val[i]], sig_av[n_val[i]], dy_av[n_val[i]], dsig_av[n_val[i]];
    for(int j = 0; j < n_val[i]; j++) {
      y_av[j] = i < 3 ? y_v[j+1][0] : y_v[j][0];
      dy_av[j] = i < 3 ? dy[j+1][0] : dy[j][0];
	  
      sig_av[j] = i < 3 ? sig[j+1][i+5] : sig[j][j < 1 ? i-3 : i+5];
      dsig_av[j] = i < 3 ? dsig[j+1][i+5] : dsig[j][j < 1 ? i-3 : i+5];
	  
    }
    g[i] = new TGraphAsymmErrors(n_val[i], y_av, sig_av, dy_av, dy_av, dsig_av, dsig_av); 
  }
  
  
}

void y_psi2_7(TGraphAsymmErrors **g, double norm, int y_n)
{
  double mass = 3.686;
  
  int n_val[y_n];
  for (int i = 0; i < y_n; i++)
    n_val[i] = 5;

  ifstream file;
  string data;
  int n_pts[5] = {11, 11, 11, 11, 11};
  double pt_lo[5][n_pts[2]], pt_hi[5][n_pts[2]], y_v[5][n_pts[2]], dy[5][n_pts[2]], sig[5][n_pts[2]], dsig[5][n_pts[2]];
  double nums[12];
  
  string filelist[5] = {"data/LHCb_psi2_7_cs_y1.txt",
			"data/LHCb_psi2_7_cs_y2.txt",
			"data/LHCb_psi2_7_cs_y3.txt",
			"data/LHCb_psi2_7_cs_y4.txt",
			"data/LHCb_psi2_7_cs_y5.txt"};
  
  for(int k = 0; k < 5; k++) {
    file.open(filelist[k].c_str());
    
    for(int j = 0; j < n_pts[k]; j++) {
      for(int i = 0; i < 12; i++)
	file >> nums[i];
      
      y_v[k][j] = 0.5 * (nums[1] + nums[0]);
      dy[k][j] = 0.5 * (nums[1] - nums[0]);
      
      pt_lo[k][j] = nums[3];
      pt_hi[k][j] = nums[4];
      
      double signorm = 1.;
      sig[k][j] = nums[5]*mass/(signorm);
      dsig[k][j] = sqrt(nums[6]*nums[6]+nums[8]*nums[8])*mass/(signorm);
    }
    file.close();
    
    // normalize to first point
    double aux;
    if(k < 1) aux = sig[0][0]/norm;
    for(int i = 0; i < n_pts[k]; i++) {
      sig[k][i] /= aux;
      dsig[k][i] /= aux;
    }
    
  }
  
  for(int i = 0; i < y_n; i++) {
    double y_av[n_val[i]], sig_av[n_val[i]], dy_av[n_val[i]], dsig_av[n_val[i]];
    for(int j = 0; j < n_val[i]; j++) {
      y_av[j] = y_v[j][0];
      dy_av[j] = dy[j][0];
      
      sig_av[j] = sig[j][i+5];
      dsig_av[j] = dsig[j][i+5];
      
    }
    g[i] = new TGraphAsymmErrors(n_val[i], y_av, sig_av, dy_av, dy_av, dsig_av, dsig_av);
  }
}

void y_ups1_7(TGraphAsymmErrors **g, double norm, int y_n)
{
  double mass = 9.46;
  
  int n_val[y_n];
  for (int i = 0; i < y_n; i++)
    n_val[i] = 5;

  ifstream file;
  string data;
  int n_pts[6] = {22, 24, 25, 24, 21, 15};
  double pt_lo[6][n_pts[2]], pt_hi[6][n_pts[2]], y_v[6][n_pts[2]], dy[6][n_pts[2]], sig[6][n_pts[2]], dsig[6][n_pts[2]];
  double nums[12];
  
  string filelist[6] = {"data/CMS_ups1_7_cs_y1.txt",
			"data/LHCb_ups1_7_cs_y1.txt",
			"data/LHCb_ups1_7_cs_y2.txt",
			"data/LHCb_ups1_7_cs_y3.txt",
			"data/LHCb_ups1_7_cs_y4.txt",
			"data/LHCb_ups1_7_cs_y5.txt"};
  
  for(int k = 0; k < 6; k++) {
    file.open(filelist[k].c_str());
    
    for(int j = 0; j < n_pts[k]; j++) {
      for(int i = 0; i < 12; i++)
	file >> nums[i];
      
      y_v[k][j] = 0.5 * (nums[1] + nums[0]);
      dy[k][j] = 0.5 * (nums[1] - nums[0]);
      
      pt_lo[k][j] = nums[3];
      pt_hi[k][j] = nums[4];
      
      double signorm = k < 1 ? 2.4*2.48e-2*1e6 : 2.48e-2*1e3;
      sig[k][j] = nums[5]*mass/(signorm);
      dsig[k][j] = sqrt(nums[6]*nums[6]+nums[8]*nums[8])*mass/(signorm);
    }
    file.close();
    
    // normalize to first point
    double aux;
    if(k < 1) aux = sig[0][0]/norm;
    for(int i = 0; i < n_pts[k]; i++) {
      sig[k][i] /= aux;
      dsig[k][i] /= aux;
    }
    
  }
  
  for(int i = 0; i < y_n; i++) {
    double y_av[n_val[i]], sig_av[n_val[i]], dy_av[n_val[i]], dsig_av[n_val[i]];
    for(int j = 0; j < n_val[i]; j++) {
      y_av[j] = y_v[j+1][0];
      dy_av[j] = dy[j+1][0];
      
      sig_av[j] = sig[j+1][i+10];
      dsig_av[j] = dsig[j+1][i+10];
      
    }
    g[i] = new TGraphAsymmErrors(n_val[i], y_av, sig_av, dy_av, dy_av, dsig_av, dsig_av);
  }
}
