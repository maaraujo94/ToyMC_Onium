/* auxiliary file with all the data-reading codes
   - xi-ordered for J/psi, Psi(2S), Ups(1S) @ 7 TeV
   - y-ordered for same
*/

// xi codes - 7 TeV
void xi_jpsi_7(TGraphAsymmErrors **g, TH1F **h)
{
  double mass = 3.097;

  ifstream file;
  string data;
  int n_pts;
  double nums[12];
  const int n_y = 9;

  string filelist[n_y] = {"data/CMS_jpsi_7_cs_y1.txt",
			  "data/CMS_jpsi_7_cs_y2.txt",
			  "data/CMS_jpsi_7_cs_y3.txt",
			  "data/CMS_jpsi_7_cs_y4.txt",
			  "data/LHCb_jpsi_7_cs_y1.txt",
			  "data/LHCb_jpsi_7_cs_y2.txt",
			  "data/LHCb_jpsi_7_cs_y3.txt",
			  "data/LHCb_jpsi_7_cs_y4.txt",
			  "data/LHCb_jpsi_7_cs_y5.txt"};
  for(int k = 0; k < n_y; k++) {
    file.open(filelist[k].c_str());
    n_pts = 0;
    while(getline(file,data)) n_pts+=1;
      
    double pt[n_pts], dpt_lo[n_pts], dpt_hi[n_pts], sig[n_pts], dsig_lo[n_pts], dsig_hi[n_pts], bins[n_pts+1];
      
    file.clear();
    file.seekg(0);
    for(int j = 0; j < n_pts; j++) {
      for(int i = 0; i < 12; i++)
	file >> nums[i];
      
      pt[j] = nums[2]/mass;
      dpt_lo[j] = (nums[2]-nums[3])/mass;
      dpt_hi[j] = (nums[4]-nums[2])/mass;

      bins[j] = nums[3]/mass;
      if(j == n_pts-1) bins[j+1] = nums[4]/mass;

      double signorm = k < 4 ? 5.961e-2*1e3 : 1;
      sig[j] = nums[5]*mass/(signorm);
      dsig_hi[j] = sqrt(nums[6]*nums[6]+nums[8]*nums[8])*mass/(signorm);
      dsig_lo[j] = sqrt(nums[7]*nums[7]+nums[9]*nums[9])*mass/(signorm);
    }
    file.close();
          
    g[k] = new TGraphAsymmErrors(n_pts, pt, sig, dpt_lo, dpt_hi, dsig_lo, dsig_hi);

    h[k] = new TH1F(Form("xi_y%d", k), Form("xi_y%d", k), n_pts, bins); 
  }
  
}

void xi_psi2_7(TGraphAsymmErrors **g, TH1F **h)
{
  double mass = 3.686;
  
  ifstream file;
  string data;
  int n_pts;
  double nums[12];
  const int n_y = 9;

  string filelist[n_y] = {"data/CMS_psi2_7_cs_y1.txt",
			  "data/CMS_psi2_7_cs_y2.txt",
			  "data/CMS_psi2_7_cs_y3.txt",
			  "data/CMS_psi2_7_cs_y4.txt",
			  "data/LHCb_psi2_7_cs_y1.txt",
			  "data/LHCb_psi2_7_cs_y2.txt",
			  "data/LHCb_psi2_7_cs_y3.txt",
			  "data/LHCb_psi2_7_cs_y4.txt",
			  "data/LHCb_psi2_7_cs_y5.txt"};
  for(int k = 0; k < n_y; k++) {
    file.open(filelist[k].c_str());
    n_pts = 0;
    while(getline(file,data)) n_pts+=1;
    
    double pt[n_pts], dpt_lo[n_pts], dpt_hi[n_pts], sig[n_pts], dsig_lo[n_pts], dsig_hi[n_pts], bins[n_pts+1];

    file.clear();
    file.seekg(0);
    for(int j = 0; j < n_pts; j++) {
      for(int i = 0; i < 12; i++)
	file >> nums[i];
      
      pt[j] = nums[2]/mass;
      dpt_lo[j] = (nums[2]-nums[3])/mass;
      dpt_hi[j] = (nums[4]-nums[2])/mass;

      bins[j] = nums[3]/mass;
      if(j == n_pts-1) bins[j+1] = nums[4]/mass;
      
      double signorm = k < 4 ? 8e-3 * 1e3 : 1.;
      sig[j] = nums[5]*mass/(signorm);
      dsig_hi[j] = sqrt(nums[6]*nums[6]+nums[8]*nums[8])*mass/(signorm);
      dsig_lo[j] = sqrt(nums[7]*nums[7]+nums[9]*nums[9])*mass/(signorm);
    }
    file.close();
          
    g[k] = new TGraphAsymmErrors(n_pts, pt, sig, dpt_lo, dpt_hi, dsig_lo, dsig_hi);

    h[k] = new TH1F(Form("xi_y%d", k), Form("xi_y%d", k), n_pts, bins); 
  }
  
}

void xi_ups1_7(TGraphAsymmErrors **g, TH1F **h)
{
  double mass = 9.46;
  
  ifstream file;
  string data;
  int n_pts;
  double nums[12];
  const int n_y = 7; 

  string filelist[n_y] = {"data/CMS_ups1_7_cs_y1.txt",
			  "data/CMS_ups1_7_cs_y2.txt",
			  "data/LHCb_ups1_7_cs_y1.txt",
			  "data/LHCb_ups1_7_cs_y2.txt",
			  "data/LHCb_ups1_7_cs_y3.txt",
			  "data/LHCb_ups1_7_cs_y4.txt",
			  "data/LHCb_ups1_7_cs_y5.txt"};

  for(int k = 0; k < n_y; k++) {
    file.open(filelist[k].c_str());
    n_pts = 0;
    while(getline(file,data)) n_pts+=1;
    
    double pt[n_pts], dpt_lo[n_pts], dpt_hi[n_pts], sig[n_pts], dsig_lo[n_pts], dsig_hi[n_pts], bins[n_pts+1];
    
    file.clear();
    file.seekg(0);
    for(int j = 0; j < n_pts; j++) {
      for(int i = 0; i < 12; i++)
	file >> nums[i];
      
      pt[j] = nums[2]/mass;
      dpt_lo[j] = (nums[2]-nums[3])/mass;
      dpt_hi[j] = (nums[4]-nums[2])/mass;

      bins[j] = nums[3]/mass;
      if(j == n_pts-1) bins[j+1] = nums[4]/mass;
      
      double signorm = k < 2 ? 1.2*2.48e-2*1e6 : 2.48e-2*1e3;
      sig[j] = nums[5]*mass/(signorm);
      dsig_hi[j] = sqrt(nums[6]*nums[6]+nums[8]*nums[8])*mass/(signorm);
      dsig_lo[j] = sqrt(nums[7]*nums[7]+nums[9]*nums[9])*mass/(signorm);
    }
    file.close();
          
    g[k] = new TGraphAsymmErrors(n_pts, pt, sig, dpt_lo, dpt_hi, dsig_lo, dsig_hi);

    h[k] = new TH1F(Form("xi_y%d", k), Form("xi_y%d", k), n_pts, bins); 
  }
  
}

void xi_ups2_7(TGraphAsymmErrors **g, TH1F **h)
{
  double mass = 10.023;
  
  ifstream file;
  string data;
  int n_pts;
  double nums[12];
  const int n_y = 7; 

  string filelist[n_y] = {"data/CMS_ups2_7_cs_y1.txt",
			  "data/CMS_ups2_7_cs_y2.txt",
			  "data/LHCb_ups2_7_cs_y1.txt",
			  "data/LHCb_ups2_7_cs_y2.txt",
			  "data/LHCb_ups2_7_cs_y3.txt",
			  "data/LHCb_ups2_7_cs_y4.txt",
			  "data/LHCb_ups2_7_cs_y5.txt"};

  for(int k = 0; k < n_y; k++) {
    file.open(filelist[k].c_str());
    n_pts = 0;
    while(getline(file,data)) n_pts+=1;
    
    double pt[n_pts], dpt_lo[n_pts], dpt_hi[n_pts], sig[n_pts], dsig_lo[n_pts], dsig_hi[n_pts], bins[n_pts+1];
    
    file.clear();
    file.seekg(0);
    for(int j = 0; j < n_pts; j++) {
      for(int i = 0; i < 12; i++)
	file >> nums[i];
      
      pt[j] = nums[2]/mass;
      dpt_lo[j] = (nums[2]-nums[3])/mass;
      dpt_hi[j] = (nums[4]-nums[2])/mass;

      bins[j] = nums[3]/mass;
      if(j == n_pts-1) bins[j+1] = nums[4]/mass;
      
      double signorm = k < 2 ? 1.2*1.93e-2*1e6 : 1.93e-2*1e3;
      sig[j] = nums[5]*mass/(signorm);
      dsig_hi[j] = sqrt(nums[6]*nums[6]+nums[8]*nums[8])*mass/(signorm);
      dsig_lo[j] = sqrt(nums[7]*nums[7]+nums[9]*nums[9])*mass/(signorm);
    }
    file.close();
          
    g[k] = new TGraphAsymmErrors(n_pts, pt, sig, dpt_lo, dpt_hi, dsig_lo, dsig_hi);

    h[k] = new TH1F(Form("xi_y%d", k), Form("xi_y%d", k), n_pts, bins); 
  }
 
}

void xi_ups3_7(TGraphAsymmErrors **g, TH1F **h)
{ 
  double mass = 10.355;
  
  ifstream file;
  string data;
  int n_pts;
  double nums[12];
  const int n_y = 7; 

  string filelist[n_y] = {"data/CMS_ups3_7_cs_y1.txt",
			  "data/CMS_ups3_7_cs_y2.txt",
			  "data/LHCb_ups3_7_cs_y1.txt",
			  "data/LHCb_ups3_7_cs_y2.txt",
			  "data/LHCb_ups3_7_cs_y3.txt",
			  "data/LHCb_ups3_7_cs_y4.txt",
			  "data/LHCb_ups3_7_cs_y5.txt"};

  for(int k = 0; k < n_y; k++) {
    file.open(filelist[k].c_str());
    n_pts = 0;
    while(getline(file,data)) n_pts+=1;
    
    double pt[n_pts], dpt_lo[n_pts], dpt_hi[n_pts], sig[n_pts], dsig_lo[n_pts], dsig_hi[n_pts], bins[n_pts+1];
    
    file.clear();
    file.seekg(0);
    for(int j = 0; j < n_pts; j++) {
      for(int i = 0; i < 12; i++)
	file >> nums[i];
      
      pt[j] = nums[2]/mass;
      dpt_lo[j] = (nums[2]-nums[3])/mass;
      dpt_hi[j] = (nums[4]-nums[2])/mass;

      bins[j] = nums[3]/mass;
      if(j == n_pts-1) bins[j+1] = nums[4]/mass;
      
      double signorm = k < 2 ? 1.2*2.18e-2*1e6 : 2.18e-2*1e3;
      sig[j] = nums[5]*mass/(signorm);
      dsig_hi[j] = sqrt(nums[6]*nums[6]+nums[8]*nums[8])*mass/(signorm);
      dsig_lo[j] = sqrt(nums[7]*nums[7]+nums[9]*nums[9])*mass/(signorm);
    }
    file.close();
          
    g[k] = new TGraphAsymmErrors(n_pts, pt, sig, dpt_lo, dpt_hi, dsig_lo, dsig_hi);

    h[k] = new TH1F(Form("xi_y%d", k), Form("xi_y%d", k), n_pts, bins); 
  }
 
 
}

// xi codes - 13 TeV
void xi_jpsi_13(TGraphAsymmErrors **g, TH1F **h)
{
  double mass = 3.097;

  ifstream file;
  string data;
  int n_pts;
  double nums[12];
  const int n_y = 9;

  string filelist[n_y] = {"data/CMS_jpsi_13_cs_y1.txt",
			  "data/CMS_jpsi_13_cs_y2.txt",
			  "data/CMS_jpsi_13_cs_y3.txt",
			  "data/CMS_jpsi_13_cs_y4.txt",
			  "data/LHCb_jpsi_13_cs_y1.txt",
			  "data/LHCb_jpsi_13_cs_y2.txt",
			  "data/LHCb_jpsi_13_cs_y3.txt",
			  "data/LHCb_jpsi_13_cs_y4.txt",
			  "data/LHCb_jpsi_13_cs_y5.txt"};
  for(int k = 0; k < n_y; k++) {
    file.open(filelist[k].c_str());
    n_pts = 0;
    while(getline(file,data)) n_pts+=1;
      
    double pt[n_pts], dpt_lo[n_pts], dpt_hi[n_pts], sig[n_pts], dsig_lo[n_pts], dsig_hi[n_pts], bins[n_pts+1];
      
    file.clear();
    file.seekg(0);
    for(int j = 0; j < n_pts; j++) {
      for(int i = 0; i < 12; i++)
	file >> nums[i];
      
      pt[j] = nums[2]/mass;
      dpt_lo[j] = (nums[2]-nums[3])/mass;
      dpt_hi[j] = (nums[4]-nums[2])/mass;

      bins[j] = nums[3]/mass;
      if(j == n_pts-1) bins[j+1] = nums[4]/mass;

      double signorm = k < 4 ? 5.961e-2*1e3 : 1.;
      sig[j] = nums[5]*mass/(signorm);
      dsig_hi[j] = sqrt(nums[6]*nums[6]+nums[8]*nums[8])*mass/(signorm);
      dsig_lo[j] = sqrt(nums[7]*nums[7]+nums[9]*nums[9])*mass/(signorm);
    }
    file.close();
          
    g[k] = new TGraphAsymmErrors(n_pts, pt, sig, dpt_lo, dpt_hi, dsig_lo, dsig_hi);

    h[k] = new TH1F(Form("xi_y%d", k), Form("xi_y%d", k), n_pts, bins); 
  }
  

}

void xi_psi2_13(TGraphAsymmErrors **g, TH1F **h)
{
  double mass = 3.686;
  
  ifstream file;
  string data;
  int n_pts;
  double nums[12];
  const int n_y = 9;

  string filelist[n_y] = {"data/CMS_psi2_13_cs_y1.txt",
			  "data/CMS_psi2_13_cs_y2.txt",
			  "data/CMS_psi2_13_cs_y3.txt",
			  "data/CMS_psi2_13_cs_y4.txt",
			  "data/LHCb_psi2_13_cs_y1.txt",
			  "data/LHCb_psi2_13_cs_y2.txt",
			  "data/LHCb_psi2_13_cs_y3.txt",
			  "data/LHCb_psi2_13_cs_y4.txt",
			  "data/LHCb_psi2_13_cs_y5.txt"};

  for(int k = 0; k < n_y; k++) {
    file.open(filelist[k].c_str());
    n_pts = 0;
    while(getline(file,data)) n_pts+=1;
      
    double pt[n_pts], dpt_lo[n_pts], dpt_hi[n_pts], sig[n_pts], dsig_lo[n_pts], dsig_hi[n_pts], bins[n_pts+1];
    
    file.clear();
    file.seekg(0);
    for(int j = 0; j < n_pts; j++) {
      for(int i = 0; i < 12; i++)
	file >> nums[i];
      
      pt[j] = nums[2]/mass;
      dpt_lo[j] = (nums[2]-nums[3])/mass;
      dpt_hi[j] = (nums[4]-nums[2])/mass;

      bins[j] = nums[3]/mass;
      if(j == n_pts-1) bins[j+1] = nums[4]/mass;
      
      double signorm = k < 4 ? 8e-3 * 1e3 : 1.;
      sig[j] = nums[5]*mass/(signorm);
      dsig_hi[j] = sqrt(nums[6]*nums[6]+nums[8]*nums[8])*mass/(signorm);
      dsig_lo[j] = sqrt(nums[7]*nums[7]+nums[9]*nums[9])*mass/(signorm);
    }
    file.close();
          
    g[k] = new TGraphAsymmErrors(n_pts, pt, sig, dpt_lo, dpt_hi, dsig_lo, dsig_hi);
    
    h[k] = new TH1F(Form("xi_y%d", k), Form("xi_y%d", k), n_pts, bins); 

  }
  
}

void xi_ups1_13(TGraphAsymmErrors **g, TH1F **h)
{
  double mass = 9.46;
  
  ifstream file;
  string data;
  int n_pts;
  double nums[12];
  const int n_y = 7; 

  string filelist[n_y] = {"data/CMS_ups1_13_cs_y1.txt",
			  "data/CMS_ups1_13_cs_y2.txt",
			  "data/LHCb_ups1_13_cs_y1.txt",
			  "data/LHCb_ups1_13_cs_y2.txt",
			  "data/LHCb_ups1_13_cs_y3.txt",
			  "data/LHCb_ups1_13_cs_y4.txt",
			  "data/LHCb_ups1_13_cs_y5.txt"};

  for(int k = 0; k < n_y; k++) {
    file.open(filelist[k].c_str());
    n_pts = 0;
    while(getline(file,data)) n_pts+=1;
    
    double pt[n_pts], dpt_lo[n_pts], dpt_hi[n_pts], sig[n_pts], dsig_lo[n_pts], dsig_hi[n_pts], bins[n_pts+1];
    
    file.clear();
    file.seekg(0);
    for(int j = 0; j < n_pts; j++) {
      for(int i = 0; i < 12; i++)
	file >> nums[i];
      
      pt[j] = nums[2]/mass;
      dpt_lo[j] = (nums[2]-nums[3])/mass;
      dpt_hi[j] = (nums[4]-nums[2])/mass;

      bins[j] = nums[3]/mass;
      if(j == n_pts-1) bins[j+1] = nums[4]/mass;
      
      double signorm = 2.48e-2*1e3;
      sig[j] = nums[5]*mass/(signorm);
      dsig_hi[j] = sqrt(nums[6]*nums[6]+nums[8]*nums[8])*mass/(signorm);
      dsig_lo[j] = sqrt(nums[7]*nums[7]+nums[9]*nums[9])*mass/(signorm);
    }
    file.close();
          
    g[k] = new TGraphAsymmErrors(n_pts, pt, sig, dpt_lo, dpt_hi, dsig_lo, dsig_hi);

    h[k] = new TH1F(Form("xi_y%d", k), Form("xi_y%d", k), n_pts, bins); 
  }

}

void xi_ups2_13(TGraphAsymmErrors **g, TH1F **h)
{
  double mass = 10.023;
  
  ifstream file;
  string data;
  int n_pts;
  double nums[12];
  const int n_y = 7; 

  string filelist[n_y] = {"data/CMS_ups2_13_cs_y1.txt",
			  "data/CMS_ups2_13_cs_y2.txt",
			  "data/LHCb_ups2_13_cs_y1.txt",
			  "data/LHCb_ups2_13_cs_y2.txt",
			  "data/LHCb_ups2_13_cs_y3.txt",
			  "data/LHCb_ups2_13_cs_y4.txt",
			  "data/LHCb_ups2_13_cs_y5.txt"};

  for(int k = 0; k < n_y; k++) {
    file.open(filelist[k].c_str());
    n_pts = 0;
    while(getline(file,data)) n_pts+=1;
    
    double pt[n_pts], dpt_lo[n_pts], dpt_hi[n_pts], sig[n_pts], dsig_lo[n_pts], dsig_hi[n_pts], bins[n_pts+1];
    
    file.clear();
    file.seekg(0);
    for(int j = 0; j < n_pts; j++) {
      for(int i = 0; i < 12; i++)
	file >> nums[i];
      
      pt[j] = nums[2]/mass;
      dpt_lo[j] = (nums[2]-nums[3])/mass;
      dpt_hi[j] = (nums[4]-nums[2])/mass;

      bins[j] = nums[3]/mass;
      if(j == n_pts-1) bins[j+1] = nums[4]/mass;
      
      double signorm = 1.93e-2*1e3;
      sig[j] = nums[5]*mass/(signorm);
      dsig_hi[j] = sqrt(nums[6]*nums[6]+nums[8]*nums[8])*mass/(signorm);
      dsig_lo[j] = sqrt(nums[7]*nums[7]+nums[9]*nums[9])*mass/(signorm);
    }
    file.close();
          
    g[k] = new TGraphAsymmErrors(n_pts, pt, sig, dpt_lo, dpt_hi, dsig_lo, dsig_hi);

    h[k] = new TH1F(Form("xi_y%d", k), Form("xi_y%d", k), n_pts, bins); 
  }

  
}

void xi_ups3_13(TGraphAsymmErrors **g, TH1F **h)
{
  double mass = 10.355;
  
  ifstream file;
  string data;
  int n_pts;
  double nums[12];
  const int n_y = 7; 

  string filelist[n_y] = {"data/CMS_ups3_13_cs_y1.txt",
			  "data/CMS_ups3_13_cs_y2.txt",
			  "data/LHCb_ups3_13_cs_y1.txt",
			  "data/LHCb_ups3_13_cs_y2.txt",
			  "data/LHCb_ups3_13_cs_y3.txt",
			  "data/LHCb_ups3_13_cs_y4.txt",
			  "data/LHCb_ups3_13_cs_y5.txt"};

  for(int k = 0; k < n_y; k++) {
    file.open(filelist[k].c_str());
    n_pts = 0;
    while(getline(file,data)) n_pts+=1;
    
    double pt[n_pts], dpt_lo[n_pts], dpt_hi[n_pts], sig[n_pts], dsig_lo[n_pts], dsig_hi[n_pts], bins[n_pts+1];
    
    file.clear();
    file.seekg(0);
    for(int j = 0; j < n_pts; j++) {
      for(int i = 0; i < 12; i++)
	file >> nums[i];
      
      pt[j] = nums[2]/mass;
      dpt_lo[j] = (nums[2]-nums[3])/mass;
      dpt_hi[j] = (nums[4]-nums[2])/mass;

      bins[j] = nums[3]/mass;
      if(j == n_pts-1) bins[j+1] = nums[4]/mass;
      
      double signorm = 2.18e-2*1e3;
      sig[j] = nums[5]*mass/(signorm);
      dsig_hi[j] = sqrt(nums[6]*nums[6]+nums[8]*nums[8])*mass/(signorm);
      dsig_lo[j] = sqrt(nums[7]*nums[7]+nums[9]*nums[9])*mass/(signorm);
    }
    file.close();
          
    g[k] = new TGraphAsymmErrors(n_pts, pt, sig, dpt_lo, dpt_hi, dsig_lo, dsig_hi);

    h[k] = new TH1F(Form("xi_y%d", k), Form("xi_y%d", k), n_pts, bins); 
  }

  
}

// y codes - 7 TeV
void y_jpsi_7(TGraphAsymmErrors **g, int y_n)
{
  double mass = 3.097;
  const int nfiles = 9;
  
  int n_val[y_n];
  for (int i = 0; i < y_n; i++){
    if(i < 3) n_val[i] = 5;
    else if(i == 3) n_val[i] = 9;
    else if(i > 3 && i < 6) n_val[i] = 8;
    else n_val[i] = 7;
  }

  ifstream file;
  string data;
  int n_pts[nfiles] = {30, 30, 30, 30, 12, 12, 12, 11, 9};
  double pt_lo[nfiles][n_pts[0]], pt_hi[nfiles][n_pts[0]], y_v[nfiles][n_pts[0]], dy[nfiles][n_pts[0]], sig[nfiles][n_pts[0]], dsig_lo[nfiles][n_pts[0]], dsig_hi[nfiles][n_pts[0]];
  double nums[12];
  
  string filelist[nfiles] = {"data/CMS_jpsi_7_cs_y1.txt",
			     "data/CMS_jpsi_7_cs_y2.txt",
			     "data/CMS_jpsi_7_cs_y3.txt",
			     "data/CMS_jpsi_7_cs_y4.txt",
			     "data/LHCb_jpsi_7_cs_y1.txt",
			     "data/LHCb_jpsi_7_cs_y2.txt",
			     "data/LHCb_jpsi_7_cs_y3.txt",
			     "data/LHCb_jpsi_7_cs_y4.txt",
			     "data/LHCb_jpsi_7_cs_y5.txt"};
  for(int k = 0; k < nfiles; k++) {
    file.open(filelist[k].c_str());
	
    for(int j = 0; j < n_pts[k]; j++) {
      for(int i = 0; i < 12; i++)
	file >> nums[i];
	  
      y_v[k][j] = 0.5 * (nums[1] + nums[0]);
      dy[k][j] = 0.5 * (nums[1] - nums[0]);
	  
      pt_lo[k][j] = nums[3];
      pt_hi[k][j] = nums[4];
	  
      double signorm = k < 4 ? 5.961e-2*1e3 : 1;
      sig[k][j] = nums[5]*mass/(signorm);
      dsig_hi[k][j] = sqrt(nums[6]*nums[6]+nums[8]*nums[8])*mass/(signorm);
      dsig_lo[k][j] = sqrt(nums[7]*nums[7]+nums[9]*nums[9])*mass/(signorm);
    }
    file.close();
    
  }
  
  for(int i = 0; i < y_n; i++) { // i reps pt bin
    double y_av[n_val[i]], sig_av[n_val[i]], dy_av[n_val[i]], dsig_lo_av[n_val[i]], dsig_hi_av[n_val[i]];
    for(int j = 0; j < n_val[i]; j++) { // j repts y bin
      if(i < 3) {
	y_av[j] = y_v[j+4][0];
	dy_av[j] = dy[j+4][0];
	
	sig_av[j] = sig[j+4][i+5];
	dsig_lo_av[j] = dsig_lo[j+4][i+5];
	dsig_hi_av[j] = dsig_hi[j+4][i+5];
      }
      else {
	y_av[j] = y_v[j][0];
	dy_av[j] = dy[j][0];
	
	sig_av[j] = sig[j][j<4 ? i-3 : i+5];
	dsig_lo_av[j] = dsig_lo[j][j<4 ? i-3 : i+5];
	dsig_hi_av[j] = dsig_hi[j][j<4 ? i-3 : i+5];
      }
    }
    g[i] = new TGraphAsymmErrors(n_val[i], y_av, sig_av, dy_av, dy_av, dsig_lo_av, dsig_hi_av); 
  }
}

void y_psi2_7(TGraphAsymmErrors **g, int y_n)
{
  double mass = 3.686;
  const int nfiles = 9;
  
  int n_val[y_n];
  for (int i = 0; i < y_n; i++)
    n_val[i] = i < 2 ? 5 : 9;

  ifstream file;
  string data;
  int n_pts[nfiles] = {18, 18, 18, 18, 11, 11, 11, 11, 11};
  double pt_lo[nfiles][n_pts[2]], pt_hi[nfiles][n_pts[2]], y_v[nfiles][n_pts[2]], dy[nfiles][n_pts[2]], sig[nfiles][n_pts[2]], dsig_lo[nfiles][n_pts[0]], dsig_hi[nfiles][n_pts[0]];
  double nums[12];
  
  string filelist[nfiles] = {"data/CMS_psi2_7_cs_y1.txt",
			     "data/CMS_psi2_7_cs_y2.txt",
			     "data/CMS_psi2_7_cs_y3.txt",
			     "data/CMS_psi2_7_cs_y4.txt",
			     "data/LHCb_psi2_7_cs_y1.txt",
			     "data/LHCb_psi2_7_cs_y2.txt",
			     "data/LHCb_psi2_7_cs_y3.txt",
			     "data/LHCb_psi2_7_cs_y4.txt",
			     "data/LHCb_psi2_7_cs_y5.txt"};
  
  for(int k = 0; k < nfiles; k++) {
    file.open(filelist[k].c_str());
    
    for(int j = 0; j < n_pts[k]; j++) {
      for(int i = 0; i < 12; i++)
	file >> nums[i];
      
      y_v[k][j] = 0.5 * (nums[1] + nums[0]);
      dy[k][j] = 0.5 * (nums[1] - nums[0]);
      
      pt_lo[k][j] = nums[3];
      pt_hi[k][j] = nums[4];
      
      double signorm = k < 4 ? 8e-3 * 1e3 : 1.;
      sig[k][j] = nums[5]*mass/(signorm);
      dsig_hi[k][j] = sqrt(nums[6]*nums[6]+nums[8]*nums[8])*mass/(signorm);
      dsig_lo[k][j] = sqrt(nums[7]*nums[7]+nums[9]*nums[9])*mass/(signorm);
    }
    file.close();
    
  }
  
  for(int i = 0; i < y_n; i++) { // i reps pt bin
    double y_av[n_val[i]], sig_av[n_val[i]], dy_av[n_val[i]], dsig_lo_av[n_val[i]], dsig_hi_av[n_val[i]];
    for(int j = 0; j < n_val[i]; j++) { // j reps y bin
      y_av[j] = i < 2 ? y_v[j+4][0] : y_v[j][0];
      dy_av[j] = i < 2 ? dy[j+4][0] : dy[j][0];
      
      sig_av[j] = i < 2 ? sig[j+4][i+5] : sig[j][j < 4 ? i-2 : i+5];
      dsig_lo_av[j] = i < 2 ? dsig_lo[j+4][i+5] : dsig_lo[j][j < 4 ? i-2 : i+5];
      dsig_hi_av[j] = i < 2 ? dsig_hi[j+4][i+5] : dsig_hi[j][j < 4 ? i-2 : i+5];

    }
    g[i] = new TGraphAsymmErrors(n_val[i], y_av, sig_av, dy_av, dy_av, dsig_lo_av, dsig_hi_av);
  }
}

void y_ups1_7(TGraphAsymmErrors **g, int y_n)
{
  double mass = 9.46;
  const int nfiles = 7;
  
  int n_val[y_n];
  for (int i = 0; i < y_n; i++)
    n_val[i] = i < 1 ? 7 : 6;

  ifstream file;
  string data;
  int n_pts[nfiles] = {22, 22, 24, 25, 24, 21, 15};
  double pt_lo[nfiles][n_pts[3]], pt_hi[nfiles][n_pts[3]], y_v[nfiles][n_pts[3]], dy[nfiles][n_pts[3]], sig[nfiles][n_pts[3]], dsig_lo[nfiles][n_pts[3]], dsig_hi[nfiles][n_pts[3]];
  double nums[12];
  
  string filelist[nfiles] = {"data/CMS_ups1_7_cs_y1.txt",
			     "data/CMS_ups1_7_cs_y2.txt",
			     "data/LHCb_ups1_7_cs_y1.txt",
			     "data/LHCb_ups1_7_cs_y2.txt",
			     "data/LHCb_ups1_7_cs_y3.txt",
			     "data/LHCb_ups1_7_cs_y4.txt",
			     "data/LHCb_ups1_7_cs_y5.txt"};
  
  for(int k = 0; k < nfiles; k++) {
    file.open(filelist[k].c_str());
    
    for(int j = 0; j < n_pts[k]; j++) {
      for(int i = 0; i < 12; i++)
	file >> nums[i];
      
      y_v[k][j] = 0.5 * (nums[1] + nums[0]);
      dy[k][j] = 0.5 * (nums[1] - nums[0]);
      
      pt_lo[k][j] = nums[3];
      pt_hi[k][j] = nums[4];
      
      double signorm = k < 2 ? 1.2*2.48e-2*1e6 : 2.48e-2*1e3;
      sig[k][j] = nums[5]*mass/(signorm);
      dsig_hi[k][j] = sqrt(nums[6]*nums[6]+nums[8]*nums[8])*mass/(signorm);
      dsig_lo[k][j] = sqrt(nums[7]*nums[7]+nums[9]*nums[9])*mass/(signorm);
    }
    file.close();
        
  }
  
  for(int i = 0; i < y_n; i++) {
    double y_av[n_val[i]], sig_av[n_val[i]], dy_av[n_val[i]], dsig_lo_av[n_val[i]], dsig_hi_av[n_val[i]];
    for(int j = 0; j < n_val[i]; j++) {
      y_av[j] = y_v[j][0];
      dy_av[j] = dy[j][0];
      
      if(j<2) { // CMS bins, 2 GeV
	sig_av[j] = sig[j][i];
	dsig_lo_av[j] = dsig_lo[j][i];
	dsig_hi_av[j] = dsig_hi[j][i];
      }
      else { // LHCb bins, 1 GeV
	sig_av[j] = 0.5*(sig[j][2*i+10]+sig[j][2*i+11]);
	dsig_lo_av[j] = 0.5*(dsig_lo[j][2*i+10]+dsig_lo[j][2*i+11]);
	dsig_hi_av[j] = 0.5*(dsig_hi[j][2*i+10]+dsig_hi[j][2*i+11]);
      }
      
    }
    g[i] = new TGraphAsymmErrors(n_val[i], y_av, sig_av, dy_av, dy_av, dsig_lo_av, dsig_hi_av);
  }
}

void y_ups2_7(TGraphAsymmErrors **g, int y_n)
{
  double mass = 10.023;
  const int nfiles = 7;
  
  int n_val[y_n];
  for (int i = 0; i < y_n; i++)
    n_val[i] = i < 1 ? 7 : 6;

  ifstream file;
  string data;
  int n_pts[nfiles] = {22, 22, 24, 25, 24, 21, 15};
  double pt_lo[nfiles][n_pts[3]], pt_hi[nfiles][n_pts[3]], y_v[nfiles][n_pts[3]], dy[nfiles][n_pts[3]], sig[nfiles][n_pts[3]], dsig_lo[nfiles][n_pts[3]], dsig_hi[nfiles][n_pts[3]];
  double nums[12];
  
  string filelist[nfiles] = {"data/CMS_ups2_7_cs_y1.txt",
			     "data/CMS_ups2_7_cs_y2.txt",
			     "data/LHCb_ups2_7_cs_y1.txt",
			     "data/LHCb_ups2_7_cs_y2.txt",
			     "data/LHCb_ups2_7_cs_y3.txt",
			     "data/LHCb_ups2_7_cs_y4.txt",
			     "data/LHCb_ups2_7_cs_y5.txt"};
  
  for(int k = 0; k < nfiles; k++) {
    file.open(filelist[k].c_str());
    
    for(int j = 0; j < n_pts[k]; j++) {
      for(int i = 0; i < 12; i++)
	file >> nums[i];
      
      y_v[k][j] = 0.5 * (nums[1] + nums[0]);
      dy[k][j] = 0.5 * (nums[1] - nums[0]);
      
      pt_lo[k][j] = nums[3];
      pt_hi[k][j] = nums[4];
      
      double signorm = k < 2 ? 1.2*1.93e-2*1e6 : 1.93e-2*1e3;
      sig[k][j] = nums[5]*mass/(signorm);
      dsig_hi[k][j] = sqrt(nums[6]*nums[6]+nums[8]*nums[8])*mass/(signorm);
      dsig_lo[k][j] = sqrt(nums[7]*nums[7]+nums[9]*nums[9])*mass/(signorm);
    }
    file.close();
        
  }
  
  for(int i = 0; i < y_n; i++) {
    double y_av[n_val[i]], sig_av[n_val[i]], dy_av[n_val[i]], dsig_lo_av[n_val[i]], dsig_hi_av[n_val[i]];
    for(int j = 0; j < n_val[i]; j++) {
      y_av[j] = y_v[j][0];
      dy_av[j] = dy[j][0];

      if(j<2) { // CMS bins, 2 GeV
	sig_av[j] = sig[j][i];
	dsig_hi_av[j] = dsig_hi[j][i];
	dsig_lo_av[j] = dsig_lo[j][i];
      }
      else { // LHCb bins, 1 GeV
	sig_av[j] = 0.5*(sig[j][2*i+10]+sig[j][2*i+11]);
	dsig_lo_av[j] = 0.5*(dsig_lo[j][2*i+10]+dsig_lo[j][2*i+11]);
	dsig_hi_av[j] = 0.5*(dsig_hi[j][2*i+10]+dsig_hi[j][2*i+11]);
      }
      
    }
    g[i] = new TGraphAsymmErrors(n_val[i], y_av, sig_av, dy_av, dy_av, dsig_lo_av, dsig_hi_av);
  }

}

void y_ups3_7(TGraphAsymmErrors **g, int y_n)
{
  double mass = 10.355;
  const int nfiles = 7;
  
  int n_val[y_n];
  for (int i = 0; i < y_n; i++)
    n_val[i] = i < 1 ? 7 : 6;

  ifstream file;
  string data;
  int n_pts[nfiles] = {22, 22, 24, 25, 24, 21, 15};
  double pt_lo[nfiles][n_pts[3]], pt_hi[nfiles][n_pts[3]], y_v[nfiles][n_pts[3]], dy[nfiles][n_pts[3]], sig[nfiles][n_pts[3]], dsig_lo[nfiles][n_pts[3]], dsig_hi[nfiles][n_pts[3]];
  double nums[12];
  
  string filelist[nfiles] = {"data/CMS_ups3_7_cs_y1.txt",
			     "data/CMS_ups3_7_cs_y2.txt",
			     "data/LHCb_ups3_7_cs_y1.txt",
			     "data/LHCb_ups3_7_cs_y2.txt",
			     "data/LHCb_ups3_7_cs_y3.txt",
			     "data/LHCb_ups3_7_cs_y4.txt",
			     "data/LHCb_ups3_7_cs_y5.txt"};
  
  for(int k = 0; k < nfiles; k++) {
    file.open(filelist[k].c_str());
    
    for(int j = 0; j < n_pts[k]; j++) {
      for(int i = 0; i < 12; i++)
	file >> nums[i];
      
      y_v[k][j] = 0.5 * (nums[1] + nums[0]);
      dy[k][j] = 0.5 * (nums[1] - nums[0]);
      
      pt_lo[k][j] = nums[3];
      pt_hi[k][j] = nums[4];
      
      double signorm = k < 2 ? 1.2*2.18e-2*1e6 : 2.18e-2*1e3;
      sig[k][j] = nums[5]*mass/(signorm);
      dsig_hi[k][j] = sqrt(nums[6]*nums[6]+nums[8]*nums[8])*mass/(signorm);
      dsig_lo[k][j] = sqrt(nums[7]*nums[7]+nums[9]*nums[9])*mass/(signorm);
    }
    file.close();
        
  }
  
  for(int i = 0; i < y_n; i++) {
    double y_av[n_val[i]], sig_av[n_val[i]], dy_av[n_val[i]], dsig_lo_av[n_val[i]], dsig_hi_av[n_val[i]];
    for(int j = 0; j < n_val[i]; j++) {
      y_av[j] = y_v[j][0];
      dy_av[j] = dy[j][0];

      if(j<2) { // CMS bins, 2 GeV
	sig_av[j] = sig[j][i];
	dsig_lo_av[j] = dsig_lo[j][i];
	dsig_hi_av[j] = dsig_hi[j][i];
      }
      else { // LHCb bins, 1 GeV
	sig_av[j] = 0.5*(sig[j][2*i+10]+sig[j][2*i+11]);
	dsig_lo_av[j] = 0.5*(dsig_lo[j][2*i+10]+dsig_lo[j][2*i+11]);
	dsig_hi_av[j] = 0.5*(dsig_hi[j][2*i+10]+dsig_hi[j][2*i+11]);
      }
      
    }
    g[i] = new TGraphAsymmErrors(n_val[i], y_av, sig_av, dy_av, dy_av, dsig_lo_av, dsig_hi_av);
  }


}

// y codes - 13 TeV

void y_jpsi_13(TGraphAsymmErrors **g, int y_n)
{
  double mass = 3.097;
  const int nfiles = 9;
  
  int n_val[y_n];
  for (int i = 0; i < y_n; i++)
    n_val[i] = 5;
    
  ifstream file;
  string data;
  int n_pts[nfiles] = {21, 21, 21, 21, 14, 14, 14, 14, 14};
  double pt_lo[nfiles][n_pts[0]], pt_hi[nfiles][n_pts[0]], y_v[nfiles][n_pts[0]], dy[nfiles][n_pts[0]], sig[nfiles][n_pts[0]], dsig_lo[nfiles][n_pts[0]], dsig_hi[nfiles][n_pts[0]];
  double nums[12];
  
  string filelist[nfiles] = {"data/CMS_jpsi_13_cs_y1.txt",
			     "data/CMS_jpsi_13_cs_y2.txt",
			     "data/CMS_jpsi_13_cs_y3.txt",
			     "data/CMS_jpsi_13_cs_y4.txt",
			     "data/LHCb_jpsi_13_cs_y1.txt",
			     "data/LHCb_jpsi_13_cs_y2.txt",
			     "data/LHCb_jpsi_13_cs_y3.txt",
			     "data/LHCb_jpsi_13_cs_y4.txt",
			     "data/LHCb_jpsi_13_cs_y5.txt"};
  for(int k = 0; k < nfiles; k++) {
    file.open(filelist[k].c_str());
	
    for(int j = 0; j < n_pts[k]; j++) {
      for(int i = 0; i < 12; i++)
	file >> nums[i];
	  
      y_v[k][j] = 0.5 * (nums[1] + nums[0]);
      dy[k][j] = 0.5 * (nums[1] - nums[0]);
	  
      pt_lo[k][j] = nums[3];
      pt_hi[k][j] = nums[4];
	  
      double signorm = k < 4 ? 5.961e-2*1e3 : 1;
      sig[k][j] = nums[5]*mass/(signorm);
      dsig_hi[k][j] = sqrt(nums[6]*nums[6]+nums[8]*nums[8])*mass/(signorm);
      dsig_lo[k][j] = sqrt(nums[7]*nums[7]+nums[9]*nums[9])*mass/(signorm);
    }
    file.close();
    
  }
  
  for(int i = 0; i < y_n; i++) { // i reps pt bin
    double y_av[n_val[i]], sig_av[n_val[i]], dy_av[n_val[i]], dsig_lo_av[n_val[i]], dsig_hi_av[n_val[i]];
    for(int j = 0; j < n_val[i]; j++) { // j repts y bin
      y_av[j] = y_v[j+4][0];
      dy_av[j] = dy[j+4][0];
	
      sig_av[j] = sig[j+4][i+7];
      dsig_lo_av[j] = dsig_lo[j+4][i+7];
      dsig_hi_av[j] = dsig_hi[j+4][i+7];
    }
    g[i] = new TGraphAsymmErrors(n_val[i], y_av, sig_av, dy_av, dy_av, dsig_lo_av, dsig_hi_av); 
  }
  

}

void y_psi2_13(TGraphAsymmErrors **g, int y_n)
{
  double mass = 3.686;
  const int nfiles = 9;
  
  int n_val[y_n];
  for (int i = 0; i < y_n; i++)
    n_val[i] = 5;

  ifstream file;
  string data;
  int n_pts[nfiles] = {9, 9, 9, 9, 17, 17, 16, 15, 14};
  double pt_lo[nfiles][n_pts[4]], pt_hi[nfiles][n_pts[4]], y_v[nfiles][n_pts[4]], dy[nfiles][n_pts[4]], sig[nfiles][n_pts[4]], dsig_lo[nfiles][n_pts[4]], dsig_hi[nfiles][n_pts[4]];
  double nums[12];
  
  string filelist[nfiles] = {"data/CMS_psi2_13_cs_y1.txt",
			     "data/CMS_psi2_13_cs_y2.txt",
			     "data/CMS_psi2_13_cs_y3.txt",
			     "data/CMS_psi2_13_cs_y4.txt",
			     "data/LHCb_psi2_13_cs_y1.txt",
			     "data/LHCb_psi2_13_cs_y2.txt",
			     "data/LHCb_psi2_13_cs_y3.txt",
			     "data/LHCb_psi2_13_cs_y4.txt",
			     "data/LHCb_psi2_13_cs_y5.txt"};
  
  for(int k = 0; k < nfiles; k++) {
    file.open(filelist[k].c_str());
    
    for(int j = 0; j < n_pts[k]; j++) {
      for(int i = 0; i < 12; i++)
	file >> nums[i];
      
      y_v[k][j] = 0.5 * (nums[1] + nums[0]);
      dy[k][j] = 0.5 * (nums[1] - nums[0]);
      
      pt_lo[k][j] = nums[3];
      pt_hi[k][j] = nums[4];
      
      double signorm = k < 4 ? 8e-3 * 1e3 : 1.;
      sig[k][j] = nums[5]*mass/(signorm);
      dsig_hi[k][j] = sqrt(nums[6]*nums[6]+nums[8]*nums[8])*mass/(signorm);
      dsig_lo[k][j] = sqrt(nums[7]*nums[7]+nums[9]*nums[9])*mass/(signorm);
    }
    file.close();
    
  }
  
  for(int i = 0; i < y_n; i++) { // i reps pt bin
    double y_av[n_val[i]], sig_av[n_val[i]], dy_av[n_val[i]], dsig_lo_av[n_val[i]], dsig_hi_av[n_val[i]];
    for(int j = 0; j < n_val[i]; j++) { // j reps y bin
      y_av[j] = y_v[j+4][0];
      dy_av[j] = dy[j+4][0];
      
      sig_av[j] = sig[j+4][i+6];
      dsig_lo_av[j] = dsig_lo[j+4][i+6];
      dsig_hi_av[j] = dsig_hi[j+4][i+6];
    }
    
    g[i] = new TGraphAsymmErrors(n_val[i], y_av, sig_av, dy_av, dy_av, dsig_lo_av, dsig_hi_av);

  }
}

void y_ups1_13(TGraphAsymmErrors **g, int y_n)
{
  double mass = 9.46;
  const int nfiles = 7;
  
  int n_val[y_n];
  for (int i = 0; i < y_n; i++)
    n_val[i] = i < 1 ? 5 : 4;

  ifstream file;
  string data;
  int n_pts[nfiles] = {17, 17, 24, 23, 24, 23, 16};
  double pt_lo[nfiles][n_pts[2]], pt_hi[nfiles][n_pts[2]], y_v[nfiles][n_pts[2]], dy[nfiles][n_pts[2]], sig[nfiles][n_pts[2]], dsig_lo[nfiles][n_pts[2]], dsig_hi[nfiles][n_pts[2]];
  double nums[12];
  
  string filelist[nfiles] = {"data/CMS_ups1_13_cs_y1.txt",
			     "data/CMS_ups1_13_cs_y2.txt",
			     "data/LHCb_ups1_13_cs_y1.txt",
			     "data/LHCb_ups1_13_cs_y2.txt",
			     "data/LHCb_ups1_13_cs_y3.txt",
			     "data/LHCb_ups1_13_cs_y4.txt",
			     "data/LHCb_ups1_13_cs_y5.txt"};
  
  for(int k = 0; k < nfiles; k++) {
    file.open(filelist[k].c_str());
    
    for(int j = 0; j < n_pts[k]; j++) {
      for(int i = 0; i < 12; i++)
	file >> nums[i];
      
      y_v[k][j] = 0.5 * (nums[1] + nums[0]);
      dy[k][j] = 0.5 * (nums[1] - nums[0]);
      
      pt_lo[k][j] = nums[3];
      pt_hi[k][j] = nums[4];
      
      double signorm = 2.48e-2*1e3;
      sig[k][j] = nums[5]*mass/(signorm);
      dsig_hi[k][j] = sqrt(nums[6]*nums[6]+nums[8]*nums[8])*mass/(signorm);
      dsig_lo[k][j] = sqrt(nums[7]*nums[7]+nums[9]*nums[9])*mass/(signorm);
    }
    file.close();
        
  }
  
  for(int i = 0; i < y_n; i++) {
    double y_av[n_val[i]], sig_av[n_val[i]], dy_av[n_val[i]], dsig_lo_av[n_val[i]], dsig_hi_av[n_val[i]];
    for(int j = 0; j < n_val[i]; j++) {
      y_av[j] = y_v[j+2][0];
      dy_av[j] = dy[j+2][0];

      sig_av[j] = 0.5*(sig[j+2][2*i+10]+sig[j+2][2*i+11]);
      dsig_lo_av[j] = 0.5*(dsig_lo[j+2][2*i+10]+dsig_lo[j+2][2*i+11]);
      dsig_hi_av[j] = 0.5*(dsig_hi[j+2][2*i+10]+dsig_hi[j+2][2*i+11]);
      
    }
    g[i] = new TGraphAsymmErrors(n_val[i], y_av, sig_av, dy_av, dy_av, dsig_lo_av, dsig_hi_av);
  }

}

void y_ups2_13(TGraphAsymmErrors **g, int y_n)
{
  double mass = 10.023;
  const int nfiles = 7;
  
  int n_val[y_n];
  for (int i = 0; i < y_n; i++)
    n_val[i] = i < 1 ? 5 : 4;

  ifstream file;
  string data;
  int n_pts[nfiles] = {17, 17, 24, 23, 24, 23, 16};
  double pt_lo[nfiles][n_pts[2]], pt_hi[nfiles][n_pts[2]], y_v[nfiles][n_pts[2]], dy[nfiles][n_pts[2]], sig[nfiles][n_pts[2]], dsig_lo[nfiles][n_pts[2]], dsig_hi[nfiles][n_pts[2]];
  double nums[12];
  
  string filelist[nfiles] = {"data/CMS_ups2_13_cs_y1.txt",
			     "data/CMS_ups2_13_cs_y2.txt",
			     "data/LHCb_ups2_13_cs_y1.txt",
			     "data/LHCb_ups2_13_cs_y2.txt",
			     "data/LHCb_ups2_13_cs_y3.txt",
			     "data/LHCb_ups2_13_cs_y4.txt",
			     "data/LHCb_ups2_13_cs_y5.txt"};
  
  for(int k = 0; k < nfiles; k++) {
    file.open(filelist[k].c_str());
    
    for(int j = 0; j < n_pts[k]; j++) {
      for(int i = 0; i < 12; i++)
	file >> nums[i];
      
      y_v[k][j] = 0.5 * (nums[1] + nums[0]);
      dy[k][j] = 0.5 * (nums[1] - nums[0]);
      
      pt_lo[k][j] = nums[3];
      pt_hi[k][j] = nums[4];
      
      double signorm = 1.93e-2*1e3;
      sig[k][j] = nums[5]*mass/(signorm);
      dsig_hi[k][j] = sqrt(nums[6]*nums[6]+nums[8]*nums[8])*mass/(signorm);
      dsig_lo[k][j] = sqrt(nums[7]*nums[7]+nums[9]*nums[9])*mass/(signorm);
    }
    file.close();
        
  }
  
  for(int i = 0; i < y_n; i++) {
    double y_av[n_val[i]], sig_av[n_val[i]], dy_av[n_val[i]], dsig_lo_av[n_val[i]], dsig_hi_av[n_val[i]];
    for(int j = 0; j < n_val[i]; j++) {
      y_av[j] = y_v[j+2][0];
      dy_av[j] = dy[j+2][0];

      sig_av[j] = 0.5*(sig[j+2][2*i+10]+sig[j+2][2*i+11]);
      dsig_lo_av[j] = 0.5*(dsig_lo[j+2][2*i+10]+dsig_lo[j+2][2*i+11]);
      dsig_hi_av[j] = 0.5*(dsig_hi[j+2][2*i+10]+dsig_hi[j+2][2*i+11]);
      
    }
    g[i] = new TGraphAsymmErrors(n_val[i], y_av, sig_av, dy_av, dy_av, dsig_lo_av, dsig_hi_av);
  }

}

void y_ups3_13(TGraphAsymmErrors **g, int y_n)
{
  double mass = 10.355;
  const int nfiles = 7;
  
  int n_val[y_n];
  for (int i = 0; i < y_n; i++)
    n_val[i] = i < 1 ? 5 : 4;

  ifstream file;
  string data;
  int n_pts[nfiles] = {17, 17, 24, 23, 24, 23, 16};
  double pt_lo[nfiles][n_pts[2]], pt_hi[nfiles][n_pts[2]], y_v[nfiles][n_pts[2]], dy[nfiles][n_pts[2]], sig[nfiles][n_pts[2]], dsig_lo[nfiles][n_pts[2]], dsig_hi[nfiles][n_pts[2]];
  double nums[12];
  
  string filelist[nfiles] = {"data/CMS_ups3_13_cs_y1.txt",
			     "data/CMS_ups3_13_cs_y2.txt",
			     "data/LHCb_ups3_13_cs_y1.txt",
			     "data/LHCb_ups3_13_cs_y2.txt",
			     "data/LHCb_ups3_13_cs_y3.txt",
			     "data/LHCb_ups3_13_cs_y4.txt",
			     "data/LHCb_ups3_13_cs_y5.txt"};
  
  for(int k = 0; k < nfiles; k++) {
    file.open(filelist[k].c_str());
    
    for(int j = 0; j < n_pts[k]; j++) {
      for(int i = 0; i < 12; i++)
	file >> nums[i];
      
      y_v[k][j] = 0.5 * (nums[1] + nums[0]);
      dy[k][j] = 0.5 * (nums[1] - nums[0]);
      
      pt_lo[k][j] = nums[3];
      pt_hi[k][j] = nums[4];
      
      double signorm = 2.18e-2*1e3;
      sig[k][j] = nums[5]*mass/(signorm);
      dsig_hi[k][j] = sqrt(nums[6]*nums[6]+nums[8]*nums[8])*mass/(signorm);
      dsig_lo[k][j] = sqrt(nums[7]*nums[7]+nums[9]*nums[9])*mass/(signorm);
    }
    file.close();
        
  }
  
  for(int i = 0; i < y_n; i++) {
    double y_av[n_val[i]], sig_av[n_val[i]], dy_av[n_val[i]], dsig_lo_av[n_val[i]], dsig_hi_av[n_val[i]];
    for(int j = 0; j < n_val[i]; j++) {
      y_av[j] = y_v[j+2][0];
      dy_av[j] = dy[j+2][0];

      sig_av[j] = 0.5*(sig[j+2][2*i+10]+sig[j+2][2*i+11]);
      dsig_lo_av[j] = 0.5*(dsig_lo[j+2][2*i+10]+dsig_lo[j+2][2*i+11]);
      dsig_hi_av[j] = 0.5*(dsig_hi[j+2][2*i+10]+dsig_hi[j+2][2*i+11]);
      
    }
    g[i] = new TGraphAsymmErrors(n_val[i], y_av, sig_av, dy_av, dy_av, dsig_lo_av, dsig_hi_av);
  }

}
