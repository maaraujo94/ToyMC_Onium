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

const double sqrts = 7000;  // write here collision energy in GeV
const double s = sqrts*sqrts;

//const double M = 3.097; // for the Jpsi
const double M = 9.46;  // for the Ups(1S)

//const double sstar_min = pow( 2. + sqrt(5.), 2.); // min and max of s^hat/M^2; here assuming pT > 2 * M (CHECK!)
const double xi_min = 1.;
const double xi_max = 50.;
const double y_max = 5.;

const double sstar_min = pow( xi_min + sqrt(1. + xi_min*xi_min), 2.);
const double sstar_max = s/(M*M);

const double x_min_gen = M*M * sstar_min / s;
const double x_max_gen = 1.;

// defining here all the fit conditions and parameters
const double beta = 2.;

const double rho = 2.;
const double delta = 0.;

// define the (gluon) PDFs:
// TODO: replace with LHAPDF
double func_gluonPDF(double x_gluon, double q2)
{
  // definition of gluon PDF as a function of x (argument) and Q^2 (parameter 0), with function starting at x_min (parameter 1)

  double x_min = q2 / s;

  double pdf = 0.;

  if ( x_gluon > x_min && x_gluon < x_max_gen )
    pdf = pdf_ct->xfxQ2(21, x_gluon, q2) / x_gluon;

  return pdf;
}
double func_quarkPDF(double x_quark, double q2)
{
  // definition of quark down PDF as a function of x (argument) and Q^2 (parameter 0), with function starting at x_min (parameter 1)

  double x_min = q2 / s;

  double pdf = 0.;

  if ( x_quark > x_min && x_quark < x_max_gen )
    pdf = pdf_ct->xfxQ2(1, x_quark, q2) / x_quark;

  return pdf;
}

// define the the partonic cross section (gluon-gluon)
// define as double* x, double* par for any that is used as TF1

// sstar dependent factor
double func_xsect_sstar( double* x, double* par)
{
   double sstar = x[0];

   double xs = 0.;

   if ( sstar > sstar_min ) {
     xs = pow(sstar, -beta);
   }

   return xs;
}

// costheta dependent factor (with sstar as parameter); better to generate linearly as well
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

void qqbarMC(){

  delete gRandom;
  gRandom = new TRandom3(0);

  TF1* xsect_sstar = new TF1( "xsect_sstar", func_xsect_sstar, sstar_min, sstar_max, 0);
  
  TFile* hfile = new TFile( "MC_res.root", "RECREATE", "qqbarMC" );

  TTree* qqbar = new TTree( "qqbar", "qqbar" );

  double x1;       qqbar->Branch( "x1",       &x1,       "x1/D" );
  double x2;       qqbar->Branch( "x2",       &x2,       "x2/D" );
  double sstar;    qqbar->Branch( "sstar",    &sstar,    "sstar/D" );
  double cosalpha; qqbar->Branch( "cosalpha", &cosalpha, "cosalpha/D" );
  double pLstar;   qqbar->Branch( "pLstar",   &pLstar,   "pLstar/D" );
  double xi;       qqbar->Branch( "xi",       &xi,       "xi/D" );
  double y;        qqbar->Branch( "y",        &y,        "y/D"  );
  double y_gg;     qqbar->Branch( "y_gg",     &y_gg,     "y_gg/D"  );

  // weights to be applied to the events when any distribution is plotted (always w_cos, never w_sstar)
  double w_gg;     qqbar->Branch( "w_gg",     &w_gg,     "w_gg/D"  );
  double w_qg;     qqbar->Branch( "w_qg",     &w_qg,     "w_qg/D"  );
  double w_cos;    qqbar->Branch( "w_cos",    &w_cos,    "w_cos/D"  );
  double w_sstar;  qqbar->Branch( "w_sstar",  &w_sstar,  "w_sstar/D"  );
  
  double pick;     qqbar->Branch( "pick",     &pick,     "pick/D"  );
    
  // counter to show progress of the generation
  const int n_step = n_events/50;
  cout << endl;
  cout << "------------------------------------------------------------" << endl;
  cout << "Progress: ";
  
  /////////////////// CYCLE OF EVENTS ////////////////////////
  for( int i_event = 1; i_event <= n_events; i_event++ ){
    
    // generate
    
    sstar = xsect_sstar->GetRandom();
    
    double s_hat = M*M * sstar;

    double x_min = s_hat / s; // min of the parton fractional momentum; this expression derives from x1*x2 = shat / s and |x| < 1
    double x_max = 1.; // is it true that x can always reach 1? CHECK

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
    
    // also define cosalpha function and extract;
    cosalpha = gRandom->Uniform(-1, 1);
    
    // x1 is necessarily the one generated UNIFORMLY, not the derived one: the Jacobian corrects the metrics of the "integration" over x1 and s_hat:  dx1 dx2 = d(x) d(s_hat) * 1/s/x, Jacobian of the transformation x = x1, s_hat = s * x1*x2
    // because we are generating logx rather than x, we get an additional x1 term in the Jacobian, so we get Jac = 1/(x1*s) * x1 = 1/s
    double fJacobian =  1./s ;

    // different weight components
    // to get the right weight for an event apply all weights that count
    // e.g. a gg event with uniformly-generated cosalpha has w = w_gg * w_cos
    w_gg = fJacobian * func_gluonPDF( x2, s_hat ) * func_gluonPDF( x1, s_hat );
    w_qg = fJacobian * func_quarkPDF( x2, s_hat ) * func_gluonPDF( x1, s_hat );
    w_cos = xsect_cosa(cosalpha, sstar);
    w_sstar = xsect_sstar->Eval(sstar);

    // other kinematic variables to be stored
    double pstar = (sstar-1)/(2*sqrt(sstar));
    
    pLstar = pstar*cosalpha;

    xi = pstar*sqrt(1-cosalpha*cosalpha);

    double Estar = (sstar+1)/(2*sqrt(sstar));

    double exp_y_hat = (Estar + pLstar)/(Estar - pLstar);

    y_gg = 0.5*log(x1/x2);
    
    y = 0.5*log(exp_y_hat * x1/x2);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // define here some "acceptance" cuts to prevent that unused events are stored in the ntuple (storing is much slower than generating)

    bool accepted = true;

    accepted = accepted && (xi > xi_min) && (xi < xi_max) ;
    //accepted = accepted && abs(y) < y_max;
    
    // store in the ntuple:

    if ( accepted ) qqbar->Fill();


    if (i_event%n_step == 0) cout << "X" << flush;
  } // end of generation loop
  
  cout << endl << endl;
  
  hfile->Write();

  delete pdf_ct;

} // end of main
