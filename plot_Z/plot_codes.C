/* auxiliary file with all the data-reading codes
   - xi-ordered for Z @ 8 TeV
     also initializes histogram for MC
*/

// xi codes - 8 TeV ATLAS
void xi_Z_8_A(TGraphAsymmErrors ***g, TH1F ***h)
{
  double mass = 91.1876;

  ifstream file;
  string data;
  int n_pts;
  double nums[12];
  const int n_y = 7;

  string filelist[n_y] = {"../data_Z/ATLAS_Z_8_cs_y0.txt",
			  "../data_Z/ATLAS_Z_8_cs_y1.txt",
			  "../data_Z/ATLAS_Z_8_cs_y2.txt",
			  "../data_Z/ATLAS_Z_8_cs_y3.txt",
			  "../data_Z/ATLAS_Z_8_cs_y4.txt",
			  "../data_Z/ATLAS_Z_8_cs_y5.txt",
			  "../data_Z/ATLAS_Z_8_cs_y6.txt"};

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

      double signorm = (nums[1]-nums[0]);
      sig[j] = nums[5]/(signorm);
      dsig_hi[j] = sqrt(nums[6]*nums[6]+nums[8]*nums[8])/(signorm);
      dsig_lo[j] = sqrt(nums[7]*nums[7]+nums[9]*nums[9])/(signorm);
    }
    file.close();
          
    g[0][k] = new TGraphAsymmErrors(n_pts, pt, sig, dpt_lo, dpt_hi, dsig_lo, dsig_hi);

    h[0][k] = new TH1F(Form("xi8A_y%d", k), Form("xi8A_y%d", k), n_pts, bins); 

  }  
}

// xi codes - 8 TeV CMS
void xi_Z_8_C(TGraphAsymmErrors ***g, TH1F ***h)
{
  double mass = 91.1876;

  ifstream file;
  string data;
  int n_pts;
  double nums[12];
  const int n_y = 6;

  string filelist[n_y] = {"../data_Z/CMS_Z_8_cs_y0.txt",
			  "../data_Z/CMS_Z_8_cs_y1.txt",
			  "../data_Z/CMS_Z_8_cs_y2.txt",
			  "../data_Z/CMS_Z_8_cs_y3.txt",
			  "../data_Z/CMS_Z_8_cs_y4.txt",
			  "../data_Z/CMS_Z_8_cs_y5.txt"};

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

      double signorm = 1.;
      sig[j] = nums[5]/(signorm);
      dsig_hi[j] = sqrt(nums[6]*nums[6]+nums[8]*nums[8])/(signorm);
      dsig_lo[j] = sqrt(nums[7]*nums[7]+nums[9]*nums[9])/(signorm);
    }
    file.close();
          
    g[0][k] = new TGraphAsymmErrors(n_pts, pt, sig, dpt_lo, dpt_hi, dsig_lo, dsig_hi);

    h[0][k] = new TH1F(Form("xi8C_y%d", k), Form("xi8C_y%d", k), n_pts, bins);

  }
  
}

// xi codes - 13 TeV ATLAS
void xi_Z_13_A(TGraphAsymmErrors ***g, TH1F ***h)
{
  double mass = 91.1876;

  ifstream file;
  string data;
  int n_pts;
  double nums[12];
  const int n_y = 1;

  string filelist[n_y] = {"../data_Z/ATLAS_Z_13_cs_y0.txt"};

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

      double signorm = 2.*(nums[1]-nums[0]);
      sig[j] = nums[5]/(signorm);
      dsig_hi[j] = sqrt(nums[6]*nums[6]+nums[8]*nums[8])/(signorm);
      dsig_lo[j] = sqrt(nums[7]*nums[7]+nums[9]*nums[9])/(signorm);
    }
    file.close();
          
    g[1][k] = new TGraphAsymmErrors(n_pts, pt, sig, dpt_lo, dpt_hi, dsig_lo, dsig_hi);

    h[1][k] = new TH1F(Form("xi13A_y%d", k), Form("xi13A_y%d", k), n_pts, bins); 

  }
  
}

// xi codes - 13 TeV CMS
void xi_Z_13_C(TGraphAsymmErrors ***g, TH1F ***h)
{
  double mass = 91.1876;

  ifstream file;
  string data;
  int n_pts;
  double nums[12];
  const int n_y = 6;

  string filelist[n_y] = {"../data_Z/CMS_Z_13_cs_y0.txt",
			  "../data_Z/CMS_Z_13_cs_y1.txt",
			  "../data_Z/CMS_Z_13_cs_y2.txt",
			  "../data_Z/CMS_Z_13_cs_y3.txt",
			  "../data_Z/CMS_Z_13_cs_y4.txt",
			  "../data_Z/CMS_Z_13_cs_y5.txt"};

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

      double signorm = 2.*(nums[1]-nums[0]);
      sig[j] = nums[5]/(signorm);
      dsig_hi[j] = sqrt(nums[6]*nums[6]+nums[8]*nums[8])/(signorm);
      dsig_lo[j] = sqrt(nums[7]*nums[7]+nums[9]*nums[9])/(signorm);
    }
    file.close();
          
    g[1][k] = new TGraphAsymmErrors(n_pts, pt, sig, dpt_lo, dpt_hi, dsig_lo, dsig_hi);

    h[1][k] = new TH1F(Form("xi13C_y%d", k), Form("xi13C_y%d", k), n_pts, bins); 

  }
  
}

