// code to generate a quarkonium MC sample for a given state, sqrt(s) and (beta,rho,delta) set

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
//const int n_events = 1e3;

double sqrts = 7000;  // write here collision energy in GeV
double s = sqrts*sqrts;

const double M = 3.097; // for the Jpsi
//const double M = 3.686; // for the Psi(2S)
//const double M = 9.46;  // for the Ups(1S)
//const double M = 10.023;  // for the Ups(2S)
//const double M = 10.355;  // for the Ups(3S)

// generation limits for xi and y
const double xi_min = 1.;
const double xi_max = 50.;
const double y_max = 5.;

const double sstar_min = pow( xi_min + sqrt(1. + xi_min*xi_min), 2.);
const double sstar_max = s/(M*M);

// defining here all the fit conditions and parameters
const double beta = 2.;
const double rho = 2.;
const double delta = 0.;

const double avgkT2 = 1.;
  
// defs for polarization
const double Mprot = 0.9382720;
const double Mlepton = 0.10566;  // (muon)
const double pbeam = sqrts/2;
const double Ebeam = sqrt(pbeam*pbeam + Mprot*Mprot);
const double gPI = TMath::Pi();

TLorentzVector beam1_LAB, beam2_LAB;

// define the PDFs:
double func_gluonPDF(double x_gluon, double q2)
{
  // definition of gluon PDF as a function of x (argument) and Q^2 (parameter 0), with function starting at x_min (parameter 1)

  double x_min = q2 / s;

  double pdf = 0.;

  if ( x_gluon > x_min && x_gluon < 1. )
    pdf = pdf_ct->xfxQ2(21, x_gluon, q2) / x_gluon;

  return pdf;
}
double func_quarkPDF(double x_quark, double q2)
{
  // definition of quark down PDF as a function of x (argument) and Q^2 (parameter 0), with function starting at x_min (parameter 1)

  double x_min = q2 / s;

  double pdf = 0.;

  if ( x_quark > x_min && x_quark < 1. )
    pdf = pdf_ct->xfxQ2(1, x_quark, q2) / x_quark;

  return pdf;
}

// define the partonic cross section
// define as double* x, double* par for any that is used as TF1

// sstar dependent factor
double func_xsect_sstar( double* x, double* par)
{
  double sstar = x[0];

  double xs = 0.;

  if ( sstar > sstar_min && sstar < sstar_max) {
    xs = pow(sstar, -beta);
  }

  return xs;
}

// costheta dependent factor (with sstar as parameter)
double xsect_cosa(double cosa, double sstar)
{
 
  double pstar = (sstar-1)/(2*sqrt(sstar));
  double Estar = (sstar+1)/(2*sqrt(sstar));

  double xc = 0;
  
  if(abs(cosa) < 1) {
    xc = ( pow( Estar+pstar*cosa, -rho ) + pow( Estar-pstar*cosa, -rho ) ) * ( pow( Estar, rho ) / 2.) * pow( 1 + 1e-5 - cosa*cosa, -delta);
  }
  
  return xc;
}

// MAIN
void qqbarMC(){

  delete gRandom;
  gRandom = new TRandom3(0);

  //initializing tlorentz
  beam1_LAB.SetPxPyPzE( 0., 0., pbeam, Ebeam);
  beam2_LAB.SetPxPyPzE( 0., 0., -pbeam, Ebeam);
  
  TF1* xsect_sstar = new TF1( "xsect_sstar", func_xsect_sstar, sstar_min, sstar_max, 0);
  
  TFile* hfile = new TFile( "MC_res.root", "RECREATE", "qqbarMC" );

  //defining the parameters to be used here
  double x1, x2, y, y_gg, sstar, cosalpha, pLstar, xi;
  
  double costh_CS, phi_CS, phith_CS;
  double costh_HX, phi_HX, phith_HX;
  double costh_PX, phi_PX, phith_PX;
  double costh_ggHX, phi_ggHX,  phith_ggHX;

  double cosalpha_inv, axAngle;
  
  // weights to be applied to the events when any distribution is plotted
  double jac, w_gg, w_qg, w_cos, w_sstar, pick;
  double w_CS_tr, w_HX_tr, w_PX_tr, w_ggHX_tr;
  double w_CS_lg, w_HX_lg, w_PX_lg, w_ggHX_lg;
  
  TTree* qqbar = new TTree( "qqbar", Form("qqbar (sqrts = %.1f TeV)", sqrts/1000.) );
  
  qqbar->Branch( "x1",       &x1,       "x1/D" );
  qqbar->Branch( "x2",       &x2,       "x2/D" );
  qqbar->Branch( "sstar",    &sstar,    "sstar/D" );
  qqbar->Branch( "cosalpha", &cosalpha, "cosalpha/D" );
  qqbar->Branch( "pLstar",   &pLstar,   "pLstar/D" );
  qqbar->Branch( "xi",       &xi,       "xi/D" );
  qqbar->Branch( "y",        &y,        "y/D"  );
  qqbar->Branch( "y_gg",     &y_gg,     "y_gg/D"  );
  
  qqbar->Branch("costh_CS",   &costh_CS,   "costh_CS/D");
  qqbar->Branch("phi_CS",     &phi_CS,     "phi_CS/D"  );
  qqbar->Branch("phith_CS",   &phith_CS,   "phith_CS/D");
  qqbar->Branch("costh_HX",   &costh_HX,   "costh_HX/D");
  qqbar->Branch("phi_HX",     &phi_HX,     "phi_HX/D"  );
  qqbar->Branch("phith_HX",   &phith_HX,   "phith_HX/D");
  qqbar->Branch("costh_PX",   &costh_PX,   "costh_PX/D");
  qqbar->Branch("phi_PX",     &phi_PX,     "phi_PX/D"  );
  qqbar->Branch("phith_PX",   &phith_PX,   "phith_PX/D");
  qqbar->Branch("costh_ggHX", &costh_ggHX, "costh_ggHX/D");
  qqbar->Branch("phi_ggHX",   &phi_ggHX,   "phi_ggHX/D"  );
  qqbar->Branch("phith_ggHX", &phith_ggHX, "phith_ggHX/D");

  qqbar->Branch("cosalpha_inv", &cosalpha_inv, "cosalpha_inv/D");
  qqbar->Branch("axAngle", &axAngle, "axAngle/D");
  
  // weights to be applied to the events when any distribution is plotted
  qqbar->Branch( "jac",       &jac,       "jac/D"  );
  qqbar->Branch( "w_gg",      &w_gg,      "w_gg/D"  );
  qqbar->Branch( "w_qg",      &w_qg,      "w_qg/D"  );
  qqbar->Branch( "w_cos",     &w_cos,     "w_cos/D"  );
  qqbar->Branch( "w_sstar",   &w_sstar,   "w_sstar/D"  );
  qqbar->Branch( "w_CS_tr",   &w_CS_tr,   "w_CS_tr/D"  );
  qqbar->Branch( "w_HX_tr",   &w_HX_tr,   "w_HX_tr/D"  );
  qqbar->Branch( "w_PX_tr",   &w_PX_tr,   "w_PX_tr/D"  );
  qqbar->Branch( "w_ggHX_tr", &w_ggHX_tr, "w_ggHX_tr/D"  );
  qqbar->Branch( "w_CS_lg",   &w_CS_lg,   "w_CS_lg/D"  );
  qqbar->Branch( "w_HX_lg",   &w_HX_lg,   "w_HX_lg/D"  );
  qqbar->Branch( "w_PX_lg",   &w_PX_lg,   "w_PX_lg/D"  );
  qqbar->Branch( "w_ggHX_lg", &w_ggHX_lg, "w_ggHX_lg/D"  );
  
  qqbar->Branch( "pick",      &pick,      "pick/D"  );
  
  
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
    // sstar from function
    sstar = xsect_sstar->GetRandom();
    double s_hat = M*M * sstar;

    // random generation of either x1 or x2
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
  
    // cosalpha from uniform;
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
    double pL = sqrt(1+xi*xi)*sinh(y)*M;

    // skip event if outside gen range or unphysical x
    if(xi < xi_min || xi > xi_max || abs(y) > y_max || x1 > 1 || x2 > 1)
      continue;
    
    // different weight components
    // to get the right weight for an event apply all weights that count
    // e.g. a gg event with uniformly-generated cosalpha has w = w_gg * w_cos

    // dx1 dx2 dt^hat = Jac * dlogx dsstar dcosalpha^hat
    //jac =  (s_hat - M*M) / ((2.*s)/(M*M)) ;
    jac = 1./s;
    w_gg = func_gluonPDF( x2, s_hat ) * func_gluonPDF( x1, s_hat );
    w_qg = func_quarkPDF( x2, s_hat ) * func_gluonPDF( x1, s_hat );

    w_cos = xsect_cosa(cosalpha, sstar);
    w_sstar = xsect_sstar->Eval(sstar);
    
    // other kinematic variables to be stored for the partons
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
  
    //get angles in each frame
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
    // define here some "acceptance" cuts to prevent that unused events are stored in the ntuple (storing is much slower than generating)
    
    bool accepted = true;

    accepted = accepted && (xi > xi_min) && (xi < xi_max) ;
    //accepted = accepted && abs(y) < y_max;
    
    // store in the ntuple:
    
    if ( accepted )
      qqbar->Fill();
    
    
  
  } // end of generation loop

  cout << endl << endl;

  hfile->Write();

  delete pdf_ct;

} // end of main
