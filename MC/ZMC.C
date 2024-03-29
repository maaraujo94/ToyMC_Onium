// code to generate a MC sample for Z boson

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

PDF *pdf_ct = mkPDF("CT14lo", 0);

// how many events to generate. The events will be weighted by the partonic cross section, so many of them will count very little.
const int n_events = 1e8;
//const int n_events = 1e3;

double sqrts = 13000;  // write here collision energy in GeV
double s = sqrts*sqrts;

const double M = 91.1876;

const double xi_min = 0.1;
const double xi_max = 40.;
const double y_max = 4.;

const double sstar_min = pow( xi_min + sqrt(1. + xi_min*xi_min), 2.);
const double sstar_max = s/(M*M);

// defining here all the fit conditions and parameters
const double avgkT2 = 1.;
  
// defs for polarization
const double Mprot = 0.9382720;
const double Mlepton = 0.10566;  // (muon)
const double pbeam = sqrts/2;
const double Ebeam = sqrt(pbeam*pbeam + Mprot*Mprot);
double gPI = TMath::Pi();

TLorentzVector beam1_LAB, beam2_LAB;

// define the (gluon) PDFs:
double func_gluonPDF(double x_gluon, double q2)
{
  // definition of gluon PDF as a function of x (argument) and Q^2 (parameter 0), with function starting at x_min (parameter 1)

  double x_min = q2 / s;

  double pdf = 0.;

  if ( x_gluon > x_min && x_gluon < 1. )
    pdf = pdf_ct->xfxQ2(21, x_gluon, q2) / x_gluon;

  return pdf;
}
double func_downquarkPDF(double x_quark, double q2)
{
  // definition of quark down PDF as a function of x (argument) and Q^2 (parameter 0), with function starting at x_min (parameter 1)

  double x_min = q2 / s;

  double pdf = 0.;

  if ( x_quark > x_min && x_quark < 1. )
    pdf = pdf_ct->xfxQ2(1, x_quark, q2) / x_quark;

  return pdf;
}
double func_upquarkPDF(double x_quark, double q2)
{
  // definition of quark up PDF as a function of x (argument) and Q^2 (parameter 0), with function starting at x_min (parameter 1)

  double x_min = q2 / s;

  double pdf = 0.;

  if ( x_quark > x_min && x_quark < 1. )
    pdf = pdf_ct->xfxQ2(2, x_quark, q2) / x_quark;

  return pdf;
}
double func_downaquarkPDF(double x_quark, double q2)
{
  // definition of anti-quark down PDF as a function of x (argument) and Q^2 (parameter 0), with function starting at x_min (parameter 1)

  double x_min = q2 / s;

  double pdf = 0.;

  if ( x_quark > x_min && x_quark < 1. )
    pdf = pdf_ct->xfxQ2(-1, x_quark, q2) / x_quark;

  return pdf;
}
double func_upaquarkPDF(double x_quark, double q2)
{
  // definition of anti-quark up PDF as a function of x (argument) and Q^2 (parameter 0), with function starting at x_min (parameter 1)

  double x_min = q2 / s;

  double pdf = 0.;

  if ( x_quark > x_min && x_quark < 1. )
    pdf = pdf_ct->xfxQ2(-2, x_quark, q2) / x_quark;

  return pdf;
}

// define the partonic cross section (quark-antiquark)
double xsect_qq(double sstar, double tstar, double ustar)
{
  return (tstar*tstar+ustar*ustar+2*sstar)/(tstar*ustar); 
}

// define the partonic cross section (quark-gluon)
double xsect_qg(double sstar, double tstar, double ustar)
{
  return (sstar*sstar+ustar*ustar+2*tstar)/(sstar*ustar);
}

void ZMC(){

  delete gRandom;
  gRandom = new TRandom3(0);

  //initializing tlorentz
  beam1_LAB.SetPxPyPzE( 0., 0., pbeam, Ebeam);
  beam2_LAB.SetPxPyPzE( 0., 0., -pbeam, Ebeam);
  
  TFile* hfile = new TFile( "MC_res_Z.root", "RECREATE", "ZMC" );

  //defining the parameters to be used here
  double x1, x2, y, y_gg, sstar, cosalpha, pLstar, xi;
  TLorentzVector *muonP = 0; // storing dimuons
  TLorentzVector *muonN = 0;
  
  double costh_CS, phi_CS, phith_CS;
  double costh_HX, phi_HX, phith_HX;
  double costh_PX, phi_PX, phith_PX;
  double costh_ggHX, phi_ggHX,  phith_ggHX;

  double cosalpha_inv, axAngle;
  
  // weights to be applied to the events when any distribution is plotted
  double jac, w_uub, w_ddb, w_ug, w_dg, w_ubg, w_dbg, pick;
  double w_CS_tr, w_HX_tr, w_PX_tr, w_ggHX_tr;
  double w_CS_lg, w_HX_lg, w_PX_lg, w_ggHX_lg;
  
  TTree* zbar = new TTree( "ZMC", Form("ZMC (sqrts = %.1f TeV)", sqrts/1000.) );
    
  zbar->Branch( "x1",       &x1,       "x1/D" );
  zbar->Branch( "x2",       &x2,       "x2/D" );
  zbar->Branch( "sstar",    &sstar,    "sstar/D" );
  zbar->Branch( "cosalpha", &cosalpha, "cosalpha/D" );
  zbar->Branch( "pLstar",   &pLstar,   "pLstar/D" );
  zbar->Branch( "xi",       &xi,       "xi/D" );
  zbar->Branch( "y",        &y,        "y/D"  );
  zbar->Branch( "y_gg",     &y_gg,     "y_gg/D"  );

  zbar->Branch("muonP", &muonP);
  zbar->Branch("muonN", &muonN);
  
  zbar->Branch("costh_CS",   &costh_CS,   "costh_CS/D");
  zbar->Branch("phi_CS",     &phi_CS,     "phi_CS/D"  );
  zbar->Branch("phith_CS",   &phith_CS,   "phith_CS/D");
  zbar->Branch("costh_HX",   &costh_HX,   "costh_HX/D");
  zbar->Branch("phi_HX",     &phi_HX,     "phi_HX/D"  );
  zbar->Branch("phith_HX",   &phith_HX,   "phith_HX/D");
  zbar->Branch("costh_PX",   &costh_PX,   "costh_PX/D");
  zbar->Branch("phi_PX",     &phi_PX,     "phi_PX/D"  );
  zbar->Branch("phith_PX",   &phith_PX,   "phith_PX/D");
  zbar->Branch("costh_ggHX", &costh_ggHX, "costh_ggHX/D");
  zbar->Branch("phi_ggHX",   &phi_ggHX,   "phi_ggHX/D"  );
  zbar->Branch("phith_ggHX", &phith_ggHX, "phith_ggHX/D");
  
  zbar->Branch("cosalpha_inv", &cosalpha_inv, "cosalpha_inv/D");
  zbar->Branch("axAngle", &axAngle, "axAngle/D");
    
  // weights to be applied to the events when any distribution is plotted
  zbar->Branch( "jac",       &jac,       "jac/D"  );
  zbar->Branch( "w_uub",     &w_uub,     "w_uub/D"  );
  zbar->Branch( "w_ddb",     &w_ddb,     "w_ddb/D"  );
  zbar->Branch( "w_ug",      &w_ug,      "w_ug/D"  );
  zbar->Branch( "w_dg",      &w_dg,      "w_dg/D"  );
  zbar->Branch( "w_ubg",     &w_ubg,     "w_ubg/D"  );
  zbar->Branch( "w_dbg",     &w_dbg,     "w_dbg/D"  );
  zbar->Branch( "w_CS_tr",   &w_CS_tr,   "w_CS_tr/D"  );
  zbar->Branch( "w_HX_tr",   &w_HX_tr,   "w_HX_tr/D"  );
  zbar->Branch( "w_PX_tr",   &w_PX_tr,   "w_PX_tr/D"  );
  zbar->Branch( "w_ggHX_tr", &w_ggHX_tr, "w_ggHX_tr/D"  );
  zbar->Branch( "w_CS_lg",   &w_CS_lg,   "w_CS_lg/D"  );
  zbar->Branch( "w_HX_lg",   &w_HX_lg,   "w_HX_lg/D"  );
  zbar->Branch( "w_PX_lg",   &w_PX_lg,   "w_PX_lg/D"  );
  zbar->Branch( "w_ggHX_lg", &w_ggHX_lg, "w_ggHX_lg/D"  );
  
  zbar->Branch( "pick",      &pick,      "pick/D"  );
  
  // counter to show progress of the generation
  const int n_step = n_events/50;
  cout << endl;
  cout << "------------------------------------------------------------" << endl;
  cout << "Progress: ";
  
  /////////////////// CYCLE OF EVENTS ////////////////////////
  double x_min = sstar_min*M*M / s; // min of the parton fractional momentum; this expression derives from x1*x2 = shat / s and |x| < 1
  double x_max = 1.; // is it true that x can always reach 1? CHECK

  for( int i_event = 1; i_event <= n_events; i_event++ ){

    if (i_event%n_step == 0) cout << "X" << flush;

    // generate
    // sstar from uniform
    sstar = gRandom->Uniform(sstar_min, sstar_max);
    double s_hat = M*M * sstar;

    // random generation of either x1 or x2
    // logx -> x from uniform
    pick = gRandom->Uniform(-1,1);
    if(pick > 0) {
      double logx1 = gRandom->Uniform(log(x_min), log(x_max));
      x1 = exp(logx1);
      x2 = s_hat / s / x1;
    }
    else {
      double logx2 = gRandom->Uniform(log(x_min), log(x_max));
      x2 = exp(logx2);
      x1 = s_hat / s / x2;
    }
    
    // cosalpha from uniform
    cosalpha = gRandom->Uniform(-1, 1);

    // quarkonium kinematic variables
    // p, pT, pL, E -> all partonic and /M
    double pstar = (sstar-1)/(2*sqrt(sstar)); 
    pLstar = pstar*cosalpha;
    double Estar = (sstar+1)/(2*sqrt(sstar));
    xi = pstar*sqrt(1-cosalpha*cosalpha);

    // y expressions + pL (pT^hat = pT)
    double y_hat = 0.5*log((Estar + pLstar)/(Estar - pLstar));
    y_gg = 0.5*log(x1/x2);
    y = y_hat + y_gg;
    double pL=sqrt(1+xi*xi)*sinh(y)*M;

    if(xi < xi_min || xi > xi_max || abs(y) > y_max || x1 > 1 || x2 > 1)
      continue;
    
    // u and t variables (partonic /M) for the weights
    double tstar = -sqrt(sstar)*pstar*(1-cosalpha);
    double ustar = -sqrt(sstar)*pstar*(1+cosalpha);
    
    // different weight components
    // from the sigma expression we get an s_hat^2
    double sigma_aux = 1./(s_hat*s_hat);

    // dx1 dx2 dt^hat = Jac * dlogx dsstar dcosalpha^hat
    // Jac = (s^hat - M^2) / (2 s / M^2)
    jac =  (s_hat-M*M)/((2.*s)/(M*M));

    // the PDFs for the partons are symmetrized in x1, x2
    double uubx = 0.5*(func_upquarkPDF(x1, s_hat)*func_upaquarkPDF(x2, s_hat)* xsect_qq(sstar, tstar, ustar)+
		       func_upaquarkPDF(x1, s_hat)*func_upquarkPDF(x2, s_hat)* xsect_qq(sstar, ustar, tstar));
    double ddbx = 0.5*(func_downquarkPDF(x1, s_hat)*func_downaquarkPDF(x2, s_hat)* xsect_qq(sstar, tstar, ustar)+
		       func_downaquarkPDF(x1, s_hat)*func_downquarkPDF(x2, s_hat)* xsect_qq(sstar, ustar, tstar));
    double ugx = 0.5*(func_upquarkPDF(x1, s_hat)*func_gluonPDF(x2, s_hat)* xsect_qg(sstar, tstar, ustar)+
		      func_gluonPDF(x1, s_hat)*func_upquarkPDF(x2, s_hat)* xsect_qg(sstar, ustar, tstar));
    double dgx = 0.5*(func_downquarkPDF(x1, s_hat)*func_gluonPDF(x2, s_hat)* xsect_qg(sstar, tstar, ustar)+
		      func_gluonPDF(x1, s_hat)*func_downquarkPDF(x2, s_hat)* xsect_qg(sstar, ustar, tstar));
    double ubgx = 0.5*(func_upaquarkPDF(x1, s_hat)*func_gluonPDF(x2, s_hat)* xsect_qg(sstar, tstar, ustar)+
		       func_gluonPDF(x1, s_hat)*func_upaquarkPDF(x2, s_hat)* xsect_qg(sstar, ustar, tstar));
    double dbgx = 0.5*(func_downaquarkPDF(x1, s_hat)*func_gluonPDF(x2, s_hat)* xsect_qg(sstar, tstar, ustar)+
		       func_gluonPDF(x1, s_hat)*func_downaquarkPDF(x2, s_hat)* xsect_qg(sstar, ustar, tstar));
     
    // the weight functions are obtained multiplying the below by the Jacobian
    w_uub = sigma_aux * uubx * 8./9.;
    w_ddb = sigma_aux * ddbx * 8./9.    ;
    w_ug  = sigma_aux * ugx  * (-1./3.) ;
    w_dg  = sigma_aux * dgx  * (-1./3.) ;
    w_ubg = sigma_aux * ubgx * (-1./3.) ;
    w_dbg = sigma_aux * dbgx * (-1./3.) ;
       

    // other kinematic variables to be stored
    // for the partons
    
    double Phi1 = 2. * gPI * gRandom->Uniform(1.) - gPI;
    double Phi2 = 2. * gPI * gRandom->Uniform(1.) - gPI;
	  
    double kT1 = abs(gRandom->Gaus(0, avgkT2));
    double kT2 = abs(gRandom->Gaus(0, avgkT2));

    TLorentzVector parton1_LAB, parton2_LAB;

    parton1_LAB.SetPxPyPzE( kT1 * cos(Phi1), kT1 * sin(Phi1), pbeam*x1, sqrt(kT1*kT1 + pbeam*x1*pbeam*x1) );
    parton2_LAB.SetPxPyPzE( kT2 * cos(Phi2), kT2 * sin(Phi2), -pbeam*x2, sqrt(kT2*kT2 + pbeam*x2*pbeam*x2) );
    
    //polarization part
    //generate decay angles
    double costh_gen;
    double sinth_gen;
    double phi_gen;

    costh_gen = -1. + 2. * gRandom->Uniform(1.);
    sinth_gen = sqrt(1-costh_gen*costh_gen);
    phi_gen   = 2. * gPI * gRandom->Uniform(1.);

    //dilepton in lab frame
    double Phi = 2. * gPI * gRandom->Uniform(1.) - gPI;

    // get angles in each frame
    TLorentzVector dilepton;
    dilepton.SetXYZM( xi*M * cos(Phi) , xi*M * sin(Phi), pL, M );

    //lepton in dilepton rest frame
    double p_lepton_DILEP = sqrt( 0.25*M*M - Mlepton*Mlepton );

    TLorentzVector lepton_DILEP;
	
    lepton_DILEP.SetXYZM( p_lepton_DILEP * sinth_gen * cos(phi_gen),
			  p_lepton_DILEP * sinth_gen * sin(phi_gen),
			  p_lepton_DILEP * costh_gen,
			  Mlepton );
	
    // reference directions to calculate angles:
	
    TVector3 lab_to_dilep = -dilepton.BoostVector();
	
    TLorentzVector beam1_DILEP = beam1_LAB;
    beam1_DILEP.Boost(lab_to_dilep);         // beam1 in the dilepton rest frame
    TLorentzVector beam2_DILEP = beam2_LAB;
    beam2_DILEP.Boost(lab_to_dilep);         // beam2 in the dilepton rest frame

    TVector3 beam1_direction     = beam1_DILEP.Vect().Unit();
    TVector3 beam2_direction     = beam2_DILEP.Vect().Unit();
    TVector3 dilep_direction     = dilepton.Vect().Unit();
    TVector3 beam1_beam2_bisect  = ( beam1_direction - beam2_direction ).Unit();
	
    // all polarization frames have the same Y axis = the normal to the plane formed by
    // the directions of the colliding hadrons
	
    TVector3 Yaxis = ( beam1_direction.Cross( beam2_direction ) ).Unit();

    // flip of y axis with rapidity
	
    if ( y < 0 ) Yaxis = - Yaxis;

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

    muonP = lepP;
    muonN = lepN;

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

    costh_CS = lepton_DILEP_rotated.CosTheta();

    phi_CS = lepton_DILEP_rotated.Phi() * 180. / gPI;

    if ( costh_CS < 0. ) phith_CS = phi_CS - 135.;
    if ( costh_CS > 0. ) phith_CS = phi_CS - 45.;
	
    if ( phith_CS < -180. ) phith_CS = 360. + phith_CS;


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

    costh_HX = lepton_DILEP_rotated.CosTheta();

    phi_HX = lepton_DILEP_rotated.Phi() * 180. / gPI;

    if ( costh_HX < 0. ) phith_HX = phi_HX - 135.;
    if ( costh_HX > 0. ) phith_HX = phi_HX - 45.;

    if ( phith_HX < -180. ) phith_HX = 360. + phith_HX;
    
	
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
	
    costh_PX = lepton_DILEP_rotated.CosTheta();

    phi_PX = lepton_DILEP_rotated.Phi() * 180. / gPI;

    if ( costh_PX < 0. ) phith_PX = phi_PX - 135.;
    if ( costh_PX > 0. ) phith_PX = phi_PX - 45.;

    if ( phith_PX < -180. ) phith_PX = 360. + phith_PX;

    //invariant cosalpha

    cosalpha_inv = sqrt( 1. - pow(costh_PX, 2.) ) * sin( lepton_DILEP_rotated.Phi() );

    //////////////////////////////////////////////////////
    // diparton helicity frame

    TLorentzVector diparton = parton1_LAB+parton2_LAB;
    TVector3 lab_to_diptn = -diparton.BoostVector();

    TLorentzVector parton1_DILEP = parton1_LAB;
    parton1_DILEP.Boost(lab_to_dilep);         // parton1 in the dilepton rest frame
    TLorentzVector parton2_DILEP = parton2_LAB;
    parton2_DILEP.Boost(lab_to_dilep);         // parton2 in the dilepton rest frame
    TLorentzVector diparton_DILEP = parton1_DILEP + parton2_DILEP;
    TLorentzVector qkn_diptn = dilepton;
    dilepton.Boost(lab_to_diptn);
    
    TVector3 parton1_direction     = parton1_DILEP.Vect().Unit();
    TVector3 parton2_direction     = parton2_DILEP.Vect().Unit();
    TVector3 dilep_direction_parton = -diparton_DILEP.Vect().Unit();
    TVector3 Yaxis_ptn = ( parton1_direction.Cross( parton2_direction ) ).Unit();
    if ( (y-y_gg) < 0 ) Yaxis_ptn = - Yaxis_ptn;

    newZaxis = dilep_direction_parton;
    TVector3 ggHXaxis = newZaxis;
    newYaxis = Yaxis_ptn;
    newXaxis = newYaxis.Cross( newZaxis );

    rotation.SetToIdentity();
    rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
    rotation.Invert();

    lepton_DILEP_rotated = lepton_DILEP_xyz.Vect();
	
    lepton_DILEP_rotated.Transform(rotation);

    costh_ggHX = lepton_DILEP_rotated.CosTheta();

    phi_ggHX = lepton_DILEP_rotated.Phi() * 180. / gPI;

    if ( costh_ggHX < 0. ) phith_ggHX = phi_ggHX - 135.;
    if ( costh_ggHX > 0. ) phith_ggHX = phi_ggHX - 45.;
	
    if ( phith_ggHX < -180. ) phith_ggHX = 360. + phith_ggHX;

    //weights
	
    w_HX_tr = 0.75 * ( 1. + costh_HX*costh_HX );
    w_HX_lg = 1.5 * ( 1. - costh_HX*costh_HX );
    w_CS_tr = 0.75 * ( 1. + costh_CS*costh_CS );
    w_CS_lg = 1.5 * ( 1. - costh_CS*costh_CS );
    w_PX_tr = 0.75 * ( 1. + costh_PX*costh_PX );
    w_PX_lg = 1.5 * ( 1. - costh_PX*costh_PX );
    w_ggHX_tr = 0.75 * ( 1. + costh_ggHX*costh_ggHX );
    w_ggHX_lg = 1.5 * ( 1. - costh_ggHX*costh_ggHX );

    axAngle = HXaxis.Angle(ggHXaxis);
        
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // Fill tree and end event
    zbar->Fill();
    
    
  } // end of generation loop
  
  cout << endl << endl;
  
  hfile->Write();

  delete pdf_ct;

} // end of main
