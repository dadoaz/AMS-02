/////////////////////////////////////////////////////////////////////////
//
// Estimación del número de electrones y positrones esperados en AMS
// 
//
// Daniel Doménech Azorín dadoa@alumni.uv,es
//
// 5.04.2021
// 
////////////////a///////////////////////////////////////////////////////

#include <TMath.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TMinuit.h>
#include "TH1.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"

//Fit function

Double_t Positron(Double_t *x, Double_t *par){
  //-----------------------------------------------------------
  //
  //  Funcion modelo: AMS positron flux, with two terms (two different power laws):
  //    diffusion and source (with exponential cutoff)
  //
  //-----------------------------------------------------------

	Double_t diffusion, source, flux;
	Double_t E       = x[0];
	Double_t phi     = par[0];
	Double_t C_d     = par[1];
 	Double_t gamma_d = par[2];
  	Double_t C_s     = par[3];
  	Double_t gamma_s = par[4];
  	Double_t inv_E_s = par[5];

  	diffusion = C_d*TMath::Power(((E + phi)/7.0),gamma_d);
  	source    = C_s*TMath::Power(((E + phi)/60.0),gamma_s)*TMath::Exp(-(E + phi)*inv_E_s);
  	flux      = (TMath::Power(E, 2)/TMath::Power(E + phi, 2))*(diffusion + source);
  
  	return flux;
}

Double_t Electron(Double_t *x, Double_t *par){
  //-----------------------------------------------------------
  //
  //  Funcion modelo: AMS electron flux, with two terms (two different power laws):
  //    diffusion and source, with a transition amplitude term
  //
  //-----------------------------------------------------------

	Double_t diffusion, source, flux;
 	Double_t E         = x[0];
  	Double_t phi       = par[0];
  	Double_t C_d       = par[1];
  	Double_t gamma_d   = par[2];
  	Double_t C_s       = par[3];
  	Double_t gamma_s   = par[4];
  	Double_t E_tr      = par[5];  
  	Double_t gamma_dif = par[6];

  	diffusion = C_d*TMath::Power(((E + phi)/20.0),gamma_d);
  	source    = C_s*TMath::Power(((E + phi)/300.0),gamma_s);
  	flux      = TMath::Power(E, 2)/TMath::Power(E + phi, 2)*TMath::Power(1+TMath::Power((E+phi)/E_tr,gamma_dif),-1)*(diffusion + source);
  
  	return flux;
}

Double_t Diffusion_p(Double_t *x, Double_t *par){
  
  Double_t diff;
  Double_t E       = x[0];
  Double_t phi     = par[0];
  Double_t C_d     = par[1];
  Double_t gamma_d = par[2];

  diff = C_d*TMath::Power(((E + phi)/7.0),gamma_d);
  
  return (TMath::Power(E, 5)/TMath::Power(E + phi, 2))*diff; 

}

Double_t Diffusion_e(Double_t *x, Double_t *par){
  
  Double_t diff;
  Double_t E          = x[0];
  Double_t phi        = par[0];
  Double_t C_d        = par[1];
  Double_t gamma_d    = par[2];
  Double_t E_tr       = par[3];
  Double_t gamma_dif  = par[4];

  diff = C_d*TMath::Power(((E + phi)/20.0),gamma_d);
  
  return (TMath::Power(E, 5)/TMath::Power(E + phi, 2))*TMath::Power(1+TMath::Power((E+phi)/E_tr,gamma_dif),-1)*diff; 

}

Double_t Source_p(Double_t *x, Double_t *par){
  
  Double_t sou;
  Double_t E       = x[0];
  Double_t phi     = par[0];
  Double_t C_s     = par[1];
  Double_t gamma_s = par[2];
  Double_t E_s     = par[3];

  sou = C_s*TMath::Power(((E + phi)/60.0),gamma_s)*TMath::Exp(-(E + phi)/E_s);
  
  return (TMath::Power(E, 5)/TMath::Power(E + phi, 2))*sou;

}

Double_t Source_e(Double_t *x, Double_t *par){
  
  Double_t sou;
  Double_t E          = x[0];
  Double_t phi        = par[0];
  Double_t C_s        = par[1];
  Double_t gamma_s    = par[2];
  Double_t E_tr       = par[3];
  Double_t gamma_dif  = par[4];

  sou = C_s*TMath::Power(((E + phi)/300.0),gamma_s);
  
  return (TMath::Power(E, 5)/TMath::Power(E + phi, 2))*TMath::Power(1+TMath::Power((E+phi)/E_tr,gamma_dif),-1)*sou;

}

Int_t flujo(){

	gStyle->SetOptStat(111111);
	gStyle->SetOptFit(1111);
	gStyle->SetPalette(kSolar);
	
	TCanvas *c1= new TCanvas("c1","Electron fit",10,10,1000,600);
	//c1->SetGrid();
	c1->SetLogx();

	TF1 *e_fl_sou = new TF1("source electron flux", Source_e, 0.1, 1500.,5);
    e_fl_sou->SetParameters(0.87, 3.96e-06, -3.14, 3.94, -2.14);
    e_fl_sou>SetParNames("#varphi_{e^{-}}", "C_{a}", "#gamma_{a}", "E_{t}", "#Delta #gamma_{t}");
    //fit->SetParameters(1.11, 6.49e-02, -4.10, 6.75e-05, -2.61, 809.);

    TF1 *e_fl_dif = new TF1("diffusion electron flux", Diffusion_e, 0.1, 1500.,5);
    e_fl_dif->SetParameters(0.87, 1.13e-02, -4.31, 3.94, -2.14);
    e_fl_dif>SetParNames("#varphi_{e^{-}}", "C_{b}", "#gamma_{b}", "E_{t}", "#Delta #gamma_{t}");
    //fit->SetParameters(1.11, 6.49e-02, -4.10, 6.75e-05, -2.61, 809.);

    TF1 *p_fl_sou = new TF1("source Positron flux", Source_p, 0.1, 1500.,4);
    p_fl_sou->SetParameters(1.10, 6.80e-05, -2.58, 1.23);
    p_fl_sou>SetParNames("#varphi_{e^{-}}", "C_{a}", "#gamma_{a}", "E_{t}", "#Delta #gamma_{t}");
    //fit->SetParameters(1.11, 6.49e-02, -4.10, 6.75e-05, -2.61, 809.);
    
    TF1 *p_fl_dif = new TF1("diffusion positron flux", Source_p, 0.1, 1500.,3);
    p_fl_dif->SetParameters(1.10, 6.51e-02, -4.07);
    p_fl_dif>SetParNames("#varphi_{e^{-}}", "C_{b}", "#gamma_{b}", "E_{t}", "#Delta #gamma_{t}");
   //fit->SetParameters(1.11, 6.49e-02, -4.10, 6.75e-05, -2.61, 809.);

	
	
	return 0;
}

