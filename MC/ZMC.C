// code to generate a MC sample for Z boson
// generates mutliple sqrt(s) values

#include "Riostream.h"
#include "TSystem.h"
#include "TString.h"
#include "TROOT.h"
#include "TTree.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRotation.h"
#include "TF1.h"
#include "TFile.h"
#include "TError.h"

#include "LHAPDF/LHAPDF.h"
using namespace LHAPDF;
using namespace std;

PDF *pdf_ct = mkPDF("CT14nnlo", 0);

// how many events to generate. The events will be weighted by the partonic cross section, so many of them will count very little.
const int n_events = 1e7;
//const int n_events = 1e4;

const int n_s_val = 3;
double sqrts[n_s_val] = {7000, 8000, 13000};  // write here collision energy in GeV
double s[n_s_val];

const double M = 91.1876;

const double xi_min = 0.;
const double xi_max = 50.;
const double y_max = 5.;

const double sstar_min = pow( xi_min + sqrt(1. + xi_min*xi_min), 2.);
//const double sstar_max = s/(M*M);
double sstar_max = sqrts[n_s_val-1]*sqrts[n_s_val-1]/(M*M);

// defining here all the fit conditions and parameters
const double beta = 2.;
const double rho = 2.;
const double delta = 0.;

const double avgkT2 = 1.;
  
// defs for polarization
double pbeam[n_s_val], Ebeam[n_s_val];
const double Mprot = 0.9382720;
const double Mlepton = 0.10566;  // (muon)
double gPI = TMath::Pi();

TLorentzVector beam1_LAB[n_s_val], beam2_LAB[n_s_val];

// define the (gluon) PDFs:
double func_gluonPDF(double x_gluon, double q2, double sval)
{
  // definition of gluon PDF as a function of x (argument) and Q^2 (parameter 0), with function starting at x_min (parameter 1)

  double x_min = q2 / sval;

  double pdf = 0.;

  if ( x_gluon > x_min && x_gluon < 1. )
    pdf = pdf_ct->xfxQ2(21, x_gluon, q2) / x_gluon;

  return pdf;
}
double func_downquarkPDF(double x_quark, double q2, double sval)
{
  // definition of quark down PDF as a function of x (argument) and Q^2 (parameter 0), with function starting at x_min (parameter 1)

  double x_min = q2 / sval;

  double pdf = 0.;

  if ( x_quark > x_min && x_quark < 1. )
    pdf = pdf_ct->xfxQ2(1, x_quark, q2) / x_quark;

  return pdf;
}
double func_upquarkPDF(double x_quark, double q2, double sval)
{
  // definition of quark up PDF as a function of x (argument) and Q^2 (parameter 0), with function starting at x_min (parameter 1)

  double x_min = q2 / sval;

  double pdf = 0.;

  if ( x_quark > x_min && x_quark < 1. )
    pdf = pdf_ct->xfxQ2(2, x_quark, q2) / x_quark;

  return pdf;
}
double func_downaquarkPDF(double x_quark, double q2, double sval)
{
  // definition of anti-quark down PDF as a function of x (argument) and Q^2 (parameter 0), with function starting at x_min (parameter 1)

  double x_min = q2 / sval;

  double pdf = 0.;

  if ( x_quark > x_min && x_quark < 1. )
    pdf = pdf_ct->xfxQ2(-1, x_quark, q2) / x_quark;

  return pdf;
}
double func_upaquarkPDF(double x_quark, double q2, double sval)
{
  // definition of anti-quark up PDF as a function of x (argument) and Q^2 (parameter 0), with function starting at x_min (parameter 1)

  double x_min = q2 / sval;

  double pdf = 0.;

  if ( x_quark > x_min && x_quark < 1. )
    pdf = pdf_ct->xfxQ2(-2, x_quark, q2) / x_quark;

  return pdf;
}

// define the partonic cross section (quark-antiquark)
double xsect_qq(double s, double t, double u)
{
  return (t*t+u*u+2*s)/(t*u); 
}

// define the partonic cross section (quark-gluon)
double xsect_qg(double s, double t, double u)
{
  return (s*s+u*u+2*t)/(s*u);
}

void ZMC(){

  delete gRandom;
  gRandom = new TRandom3(0);

  //initializing arrays of sqrt(s) values
  for(int i = 0; i < n_s_val; i++) {
    s[i] = sqrts[i] * sqrts[i];
    pbeam[i] = sqrts[i]/2.;
    Ebeam[i] = sqrt(pbeam[i]*pbeam[i] + Mprot*Mprot);
    
    beam1_LAB[i].SetPxPyPzE( 0., 0., pbeam[i], Ebeam[i]);
    beam2_LAB[i].SetPxPyPzE( 0., 0., -pbeam[i], Ebeam[i]);
  }
  
  TFile* hfile = new TFile( "MC_res_Z.root", "RECREATE", "ZMC" );

  //defining the parameters to be used here
  double x1[n_s_val], x2[n_s_val], y[n_s_val], y_gg[n_s_val];
  double sstar, cosalpha, pLstar, xi;
  
  double costh_CS[n_s_val], phi_CS[n_s_val], phith_CS[n_s_val];
  double costh_HX[n_s_val], phi_HX[n_s_val], phith_HX[n_s_val];
  double costh_PX[n_s_val], phi_PX[n_s_val], phith_PX[n_s_val];
  double costh_ggHX[n_s_val], phi_ggHX[n_s_val],  phith_ggHX[n_s_val];

  double cosalpha_inv[n_s_val], axAngle[n_s_val];
  
  // weights to be applied to the events when any distribution is plotted
  double w_uub[n_s_val], w_ddb[n_s_val], w_ug[n_s_val], w_dg[n_s_val], w_ubg[n_s_val], w_dbg[n_s_val];
  double pick;
  double w_CS_tr[n_s_val], w_HX_tr[n_s_val], w_PX_tr[n_s_val], w_ggHX_tr[n_s_val];
  double w_CS_lg[n_s_val], w_HX_lg[n_s_val], w_PX_lg[n_s_val], w_ggHX_lg[n_s_val];
  
  TTree** zbar = new TTree*[n_s_val];
  for(int i = 0; i < n_s_val; i++) {
    zbar[i] = new TTree( Form("zbar_%d", i), Form("zbar (sqrts = %.1f TeV)", sqrts[i]/1000.) );
    
    zbar[i]->Branch( "x1",       &x1[i],       "x1/D" );
    zbar[i]->Branch( "x2",       &x2[i],       "x2/D" );
    zbar[i]->Branch( "sstar",    &sstar,    "sstar/D" );
    zbar[i]->Branch( "cosalpha", &cosalpha, "cosalpha/D" );
    zbar[i]->Branch( "pLstar",   &pLstar,   "pLstar/D" );
    zbar[i]->Branch( "xi",       &xi,       "xi/D" );
    zbar[i]->Branch( "y",        &y[i],        "y/D"  );
    zbar[i]->Branch( "y_gg",     &y_gg[i],     "y_gg/D"  );

    /*    zbar[i]->Branch("costh_CS",   &costh_CS[i],   "costh_CS/D");
    zbar[i]->Branch("phi_CS",     &phi_CS[i],     "phi_CS/D"  );
    zbar[i]->Branch("phith_CS",   &phith_CS[i],   "phith_CS/D");
    zbar[i]->Branch("costh_HX",   &costh_HX[i],   "costh_HX/D");
    zbar[i]->Branch("phi_HX",     &phi_HX[i],     "phi_HX/D"  );
    zbar[i]->Branch("phith_HX",   &phith_HX[i],   "phith_HX/D");
    zbar[i]->Branch("costh_PX",   &costh_PX[i],   "costh_PX/D");
    zbar[i]->Branch("phi_PX",     &phi_PX[i],     "phi_PX/D"  );
    zbar[i]->Branch("phith_PX",   &phith_PX[i],   "phith_PX/D");
    zbar[i]->Branch("costh_ggHX", &costh_ggHX[i], "costh_ggHX/D");
    zbar[i]->Branch("phi_ggHX",   &phi_ggHX[i],   "phi_ggHX/D"  );
    zbar[i]->Branch("phith_ggHX", &phith_ggHX[i], "phith_ggHX/D");

    zbar[i]->Branch("cosalpha_inv", &cosalpha_inv[i], "cosalpha_inv/D");
    zbar[i]->Branch("axAngle", &axAngle[i], "axAngle/D");*/
  
    // weights to be applied to the events when any distribution is plotted (always w_cos, never w_sstar)
    zbar[i]->Branch( "w_uub",      &w_uub[i],      "w_uub/D"  );
    zbar[i]->Branch( "w_ddb",      &w_ddb[i],      "w_ddb/D"  );
    zbar[i]->Branch( "w_ug",      &w_ug[i],      "w_ug/D"  );
    zbar[i]->Branch( "w_dg",      &w_dg[i],      "w_dg/D"  );
    zbar[i]->Branch( "w_ubg",      &w_ubg[i],      "w_ubg/D"  );
    zbar[i]->Branch( "w_dbg",      &w_dbg[i],      "w_dbg/D"  );
    /*    zbar[i]->Branch( "w_CS_tr",   &w_CS_tr[i],   "w_CS_tr/D"  );
    zbar[i]->Branch( "w_HX_tr",   &w_HX_tr[i],   "w_HX_tr/D"  );
    zbar[i]->Branch( "w_PX_tr",   &w_PX_tr[i],   "w_PX_tr/D"  );
    zbar[i]->Branch( "w_ggHX_tr", &w_ggHX_tr[i], "w_ggHX_tr/D"  );
    zbar[i]->Branch( "w_CS_lg",   &w_CS_lg[i],   "w_CS_lg/D"  );
    zbar[i]->Branch( "w_HX_lg",   &w_HX_lg[i],   "w_HX_lg/D"  );
    zbar[i]->Branch( "w_PX_lg",   &w_PX_lg[i],   "w_PX_lg/D"  );
    zbar[i]->Branch( "w_ggHX_lg", &w_ggHX_lg[i], "w_ggHX_lg/D"  );*/
  
    zbar[i]->Branch( "pick",     &pick,     "pick/D"  );
  }
  
  // counter to show progress of the generation
  const int n_step = n_events/50;
  cout << endl;
  cout << "------------------------------------------------------------" << endl;
  cout << "Progress: ";
  
  /////////////////// CYCLE OF EVENTS ////////////////////////
  for( int i_event = 1; i_event <= n_events; i_event++ ){
    
    // generate
    // sstar from uniform
    sstar = gRandom->Uniform(sstar_min, sstar_max);
    double s_hat = M*M * sstar;
    double x_max = 1.; // is it true that x can always reach 1? CHECK

    // random generation of either x1 or x2
    // logx from uniform
    pick = gRandom->Uniform(-1,1);
    for(int i = 0; i < n_s_val; i++) {
      double x_min = s_hat / s[i]; // min of the parton fractional momentum; this expression derives from x1*x2 = shat / s and |x| < 1
      if(pick > 0) {
	double logx1 = gRandom->Uniform(log(x_min), log(x_max));
	x1[i] = exp(logx1);
	x2[i] = s_hat / s[i] / x1[i];
      }
      else {
	double logx2 = gRandom->Uniform(log(x_min), log(x_max));
	x2[i] = exp(logx2);
	x1[i] = s_hat / s[i] / x2[i];
      }
    }
    
    // cosalpha from uniform
    cosalpha = gRandom->Uniform(-1, 1);

    // kinematic variables

    // p, pT, pL, E -> all partonic and /M
    double pstar = (sstar-1)/(2*sqrt(sstar)); 
    pLstar = pstar*cosalpha;
    double Estar = (sstar+1)/(2*sqrt(sstar));
    xi = pstar*sqrt(1-cosalpha*cosalpha);

    // y expressions + pL (pT^hat = pT)
    double exp_y_hat = (Estar + pLstar)/(Estar - pLstar);
    double pL[n_s_val];
    for(int i=0; i< n_s_val; i++) {
      y_gg[i] = 0.5*log(x1[i]/x2[i]);
      y[i] = 0.5*log(exp_y_hat * x1[i]/x2[i]);
      pL[i]=sqrt(1+xi*xi)*sinh(y[i])*M;
    }

    // u and t variables (partonic /M) for the weights
    double tstar = -sqrt(sstar)*pstar*(1-cosalpha);
    double ustar = -sqrt(sstar)*pstar*(1+cosalpha);
    
    // different weight components
    // from the sigma expression we get an s_hat
    double sigma_aux = 1./(s_hat*s_hat);
    for(int i = 0; i < n_s_val; i++) {
      // from the change of vars from dx1 dx2 to dlogx dshat, we get a Jacobian
      double fJacobian =  (s_hat-M*M)/(2.*s[i]) ;
      // the PDFs for the partons are given by x1(u,d), x2(ubar, dbar, g)
      double uubx = 0.5*(func_upquarkPDF(x1[i], s_hat, s[i])*func_upaquarkPDF(x2[i], s_hat, s[i])+func_upaquarkPDF(x1[i], s_hat, s[i])*func_upquarkPDF(x2[i], s_hat, s[i]));
      double ddbx = 0.5*(func_downquarkPDF(x1[i], s_hat, s[i])*func_downaquarkPDF(x2[i], s_hat, s[i])+func_downaquarkPDF(x1[i], s_hat, s[i])*func_downquarkPDF(x2[i], s_hat, s[i]));
      double ugx = 0.5*(func_upquarkPDF(x1[i], s_hat, s[i])*func_gluonPDF(x2[i], s_hat, s[i])+func_gluonPDF(x1[i], s_hat, s[i])*func_upquarkPDF(x2[i], s_hat, s[i]));
      double dgx = 0.5*(func_downquarkPDF(x1[i], s_hat, s[i])*func_gluonPDF(x2[i], s_hat, s[i])+func_gluonPDF(x1[i], s_hat, s[i])*func_downquarkPDF(x2[i], s_hat, s[i]));
      double ubgx = 0.5*(func_upaquarkPDF(x1[i], s_hat, s[i])*func_gluonPDF(x2[i], s_hat, s[i])+func_gluonPDF(x1[i], s_hat, s[i])*func_upaquarkPDF(x2[i], s_hat, s[i]));
      double dbgx = 0.5*(func_downaquarkPDF(x1[i], s_hat, s[i])*func_gluonPDF(x2[i], s_hat, s[i])+func_gluonPDF(x1[i], s_hat, s[i])*func_downaquarkPDF(x2[i], s_hat, s[i]));
      
      // the full weight functions are the product of all this
      w_uub[i] = sigma_aux * fJacobian * uubx * 8./9. * xsect_qq(sstar, tstar, ustar);
      w_ddb[i] = sigma_aux * fJacobian * ddbx * 8./9. * xsect_qq(sstar, tstar, ustar);
      w_ug[i] = sigma_aux * fJacobian * ugx * (-1./3.) * xsect_qg(sstar, tstar, ustar);
      w_dg[i] = sigma_aux * fJacobian * dgx * (-1./3.) * xsect_qg(sstar, tstar, ustar);
      w_ubg[i] = sigma_aux * fJacobian * ubgx * (-1./3.) * xsect_qg(sstar, tstar, ustar);
      w_dbg[i] = sigma_aux * fJacobian * dbgx * (-1./3.) * xsect_qg(sstar, tstar, ustar);
    }

    // other kinematic variables to be stored
    // for the partons
    /*    double Phi1 = 2. * gPI * gRandom->Uniform(1.) - gPI;
    double Phi2 = 2. * gPI * gRandom->Uniform(1.) - gPI;

    double kT1 = abs(gRandom->Gaus(0, avgkT2));
    double kT2 = abs(gRandom->Gaus(0, avgkT2));

    TLorentzVector parton1_LAB[n_s_val], parton2_LAB[n_s_val];

    for(int i = 0; i < n_s_val; i++)
      {
	parton1_LAB[i].SetPxPyPzE( kT1 * cos(Phi1), kT1 * sin(Phi1), pbeam[i]*x1[i], sqrt(kT1*kT1 + pbeam[i]*x1[i]*pbeam[i]*x1[i]) );
	parton2_LAB[i].SetPxPyPzE( kT2 * cos(Phi2), kT2 * sin(Phi2), -pbeam[i]*x2[i], sqrt(kT2*kT2 + pbeam[i]*x2[i]*pbeam[i]*x2[i]) );
	}*/
    
    
    //polarization part
    //generate decay angles
    /*    double costh_gen;
    double sinth_gen;
    double phi_gen;

    costh_gen = -1. + 2. * gRandom->Uniform(1.);
    sinth_gen = sqrt(1-costh_gen*costh_gen);
    phi_gen   = 2. * gPI * gRandom->Uniform(1.);

    //dilepton in lab frame
    double Phi = 2. * gPI * gRandom->Uniform(1.) - gPI;

    //cycle over each sqrt(s) to get angles in each frame
    for(int i = 0; i < n_s_val; i++)
      {
	TLorentzVector dilepton;
	dilepton.SetXYZM( xi*M * cos(Phi) , xi*M * sin(Phi), pL[i], M );

	//lepton in dilepton rest frame
	double p_lepton_DILEP = sqrt( 0.25*M*M - Mlepton*Mlepton );

	TLorentzVector lepton_DILEP;
	
	lepton_DILEP.SetXYZM( p_lepton_DILEP * sinth_gen * cos(phi_gen),
			      p_lepton_DILEP * sinth_gen * sin(phi_gen),
			      p_lepton_DILEP * costh_gen,
			      Mlepton );
	
	// reference directions to calculate angles:
	
	TVector3 lab_to_dilep = -dilepton.BoostVector();
	
	TLorentzVector beam1_DILEP = beam1_LAB[i];
	beam1_DILEP.Boost(lab_to_dilep);         // beam1 in the dilepton rest frame
	TLorentzVector beam2_DILEP = beam2_LAB[i];
	beam2_DILEP.Boost(lab_to_dilep);         // beam2 in the dilepton rest frame

	TVector3 beam1_direction     = beam1_DILEP.Vect().Unit();
	TVector3 beam2_direction     = beam2_DILEP.Vect().Unit();
	TVector3 dilep_direction     = dilepton.Vect().Unit();
	TVector3 beam1_beam2_bisect  = ( beam1_direction - beam2_direction ).Unit();
	
	// all polarization frames have the same Y axis = the normal to the plane formed by
	// the directions of the colliding hadrons
	
	TVector3 Yaxis = ( beam1_direction.Cross( beam2_direction ) ).Unit();

	// flip of y axis with rapidity
	
	if ( y[i] < 0 ) Yaxis = - Yaxis;

	TVector3 perpendicular_to_beam = ( beam1_beam2_bisect.Cross( Yaxis ) ).Unit();

	TLorentzVector lepton_DILEP_xyz = lepton_DILEP;

	// lepton 4-vectors in the LAB frame: CHECK FOR BUGS

	TVector3 dilep_to_lab = dilepton.BoostVector();

	TLorentzVector* lepP = new TLorentzVector(0.,0.,0.,0.);
	TLorentzVector* lepN = new TLorentzVector(0.,0.,0.,0.);

	*lepP = lepton_DILEP_xyz;
	lepP->Boost(dilep_to_lab);
	lepN->SetPxPyPzE(-lepton_DILEP_xyz.Px(),-lepton_DILEP_xyz.Py(),-lepton_DILEP_xyz.Pz(),lepton_DILEP_xyz.E());
	lepN->Boost(dilep_to_lab);

	/////////////////////////////////////////////////////////////////////
	// CS frame
	
	TRotation rotation;
	
	TVector3 newZaxis = beam1_beam2_bisect;
	TVector3 newYaxis = Yaxis;
	TVector3 newXaxis = newYaxis.Cross( newZaxis );
	
	rotation.SetToIdentity();
	rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
	rotation.Invert();  // transforms coordinates from the "xyz" frame to the new frame

	TVector3 lepton_DILEP_rotated = lepton_DILEP_xyz.Vect();

	lepton_DILEP_rotated.Transform(rotation);

	costh_CS[i] = lepton_DILEP_rotated.CosTheta();

	phi_CS[i] = lepton_DILEP_rotated.Phi() * 180. / gPI;

	if ( costh_CS[i] < 0. ) phith_CS[i] = phi_CS[i] - 135.;
	if ( costh_CS[i] > 0. ) phith_CS[i] = phi_CS[i] - 45.;
	
	if ( phith_CS[i] < -180. ) phith_CS[i] = 360. + phith_CS[i];


	/////////////////////////////////////////////////////////////////////
	// HELICITY frame

	newZaxis = dilep_direction;
	TVector3 HXaxis = newZaxis;
	newYaxis = Yaxis;
	newXaxis = newYaxis.Cross( newZaxis );

	rotation.SetToIdentity();
	rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
	rotation.Invert();

	lepton_DILEP_rotated = lepton_DILEP_xyz.Vect();
	
	lepton_DILEP_rotated.Transform(rotation);

	costh_HX[i] = lepton_DILEP_rotated.CosTheta();

	phi_HX[i] = lepton_DILEP_rotated.Phi() * 180. / gPI;

	if ( costh_HX[i] < 0. ) phith_HX[i] = phi_HX[i] - 135.;
	if ( costh_HX[i] > 0. ) phith_HX[i] = phi_HX[i] - 45.;

	if ( phith_HX[i] < -180. ) phith_HX[i] = 360. + phith_HX[i];
    
	
	/////////////////////////////////////////////////////////////////////
	// PERPENDICULAR HELICITY frame

	newZaxis = perpendicular_to_beam;
	newYaxis = Yaxis;
	newXaxis = newYaxis.Cross( newZaxis );

	rotation.SetToIdentity();
	rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
	rotation.Invert();

	lepton_DILEP_rotated = lepton_DILEP_xyz.Vect();
	
	lepton_DILEP_rotated.Transform(rotation);
	
	costh_PX[i] = lepton_DILEP_rotated.CosTheta();

	phi_PX[i] = lepton_DILEP_rotated.Phi() * 180. / gPI;

	if ( costh_PX[i] < 0. ) phith_PX[i] = phi_PX[i] - 135.;
	if ( costh_PX[i] > 0. ) phith_PX[i] = phi_PX[i] - 45.;

	if ( phith_PX[i] < -180. ) phith_PX[i] = 360. + phith_PX[i];

	//invariant cosalpha

	cosalpha_inv[i] = sqrt( 1. - pow(costh_PX[i], 2.) ) * sin( lepton_DILEP_rotated.Phi() );

	//////////////////////////////////////////////////////
	// diparton helicity frame

	TLorentzVector diparton = parton1_LAB[i]+parton2_LAB[i];
	TVector3 lab_to_diptn = -diparton.BoostVector();

	TLorentzVector parton1_DILEP = parton1_LAB[i];
	parton1_DILEP.Boost(lab_to_dilep);         // parton1 in the dilepton rest frame
	TLorentzVector parton2_DILEP = parton2_LAB[i];
	parton2_DILEP.Boost(lab_to_dilep);         // parton2 in the dilepton rest frame
	TLorentzVector diparton_DILEP = parton1_DILEP + parton2_DILEP;
	TLorentzVector qkn_diptn = dilepton;
	dilepton.Boost(lab_to_diptn);
    
	TVector3 parton1_direction     = parton1_DILEP.Vect().Unit();
	TVector3 parton2_direction     = parton2_DILEP.Vect().Unit();
	TVector3 dilep_direction_parton = -diparton_DILEP.Vect().Unit();
	TVector3 Yaxis_ptn = ( parton1_direction.Cross( parton2_direction ) ).Unit();
	if ( (y[i]-y_gg[i]) < 0 ) Yaxis_ptn = - Yaxis_ptn;

	newZaxis = dilep_direction_parton;
	TVector3 ggHXaxis = newZaxis;
	newYaxis = Yaxis_ptn;
	newXaxis = newYaxis.Cross( newZaxis );

	rotation.SetToIdentity();
	rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
	rotation.Invert();

	lepton_DILEP_rotated = lepton_DILEP_xyz.Vect();
	
	lepton_DILEP_rotated.Transform(rotation);

	costh_ggHX[i] = lepton_DILEP_rotated.CosTheta();

	phi_ggHX[i] = lepton_DILEP_rotated.Phi() * 180. / gPI;

	if ( costh_ggHX[i] < 0. ) phith_ggHX[i] = phi_ggHX[i] - 135.;
	if ( costh_ggHX[i] > 0. ) phith_ggHX[i] = phi_ggHX[i] - 45.;
	
	if ( phith_ggHX[i] < -180. ) phith_ggHX[i] = 360. + phith_ggHX[i];

	//weights
	
	w_HX_tr[i] = 0.75 * ( 1. + costh_HX[i]*costh_HX[i] );
	w_HX_lg[i] = 1.5 * ( 1. - costh_HX[i]*costh_HX[i] );
	w_CS_tr[i] = 0.75 * ( 1. + costh_CS[i]*costh_CS[i] );
	w_CS_lg[i] = 1.5 * ( 1. - costh_CS[i]*costh_CS[i] );
	w_PX_tr[i] = 0.75 * ( 1. + costh_PX[i]*costh_PX[i] );
	w_PX_lg[i] = 1.5 * ( 1. - costh_PX[i]*costh_PX[i] );
	w_ggHX_tr[i] = 0.75 * ( 1. + costh_ggHX[i]*costh_ggHX[i] );
	w_ggHX_lg[i] = 1.5 * ( 1. - costh_ggHX[i]*costh_ggHX[i] );

	axAngle[i] = HXaxis.Angle(ggHXaxis);
	}*/
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // define here some "acceptance" cuts to prevent that unused events are stored in the ntuple (storing is much slower than generating)
    
    bool accepted = true;
    
    accepted = accepted && (xi > xi_min) && (xi < xi_max) ;
    //accepted = accepted && abs(y) < y_max;
    
    // store in the ntuple:
    
    if ( accepted )
      for(int i = 0; i< n_s_val; i++)
	zbar[i]->Fill();
    
    
    if (i_event%n_step == 0) cout << "X" << flush;
    
  } // end of generation loop
  
  cout << endl << endl;
  
  hfile->Write();

  delete pdf_ct;

} // end of main
