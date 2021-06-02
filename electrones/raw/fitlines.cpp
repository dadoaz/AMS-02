/////////////////////////////////////////////////////////////////////////
//
// Ajste de los datos del flujo de electrones medidos por AMS
// separando en dos componentes: término de difusión y término fuente
//
// Daniel Doménech Azorín dadoa@alumni.uv,es
//
// 28.03.2021
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
  Double_t E_s     = par[5];

  diffusion = C_d*TMath::Power(((E + phi)/7.0),gamma_d);
  source    = C_s*TMath::Power(((E + phi)/60.0),gamma_s)*TMath::Exp(-(E + phi)/E_s);
  flux      = (TMath::Power(E, 5)/TMath::Power(E + phi, 2))*(diffusion + source);
  
  return flux;
}

Double_t Electron(Double_t *x, Double_t *par){
  //-----------------------------------------------------------
  //
  //  Funcion modelo: AMS positron flux, with two terms (two different power laws):
  //    diffusion and source (with exponential cutoff)
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
  flux      = TMath::Power(E, 5)/TMath::Power(E + phi, 2)*TMath::Power(1+TMath::Power((E+phi)/E_tr,gamma_dif),-1)*(diffusion + source);
  
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

Double_t flux(Double_t *x, Double_t *par){
  //-----------------------------------------------------------
  //
  //  Funcion modelo: AMS flux in a general form for positrons and electrons, with two terms (two different power laws):
  //    diffusion and source (with exponential cutoff)
  //
  //-----------------------------------------------------------

	Double_t diffusion, source, flux;
	Double_t E         = x[0];
	Double_t phi       = par[0];
	Double_t C_d       = par[1];
 	Double_t gamma_d   = par[2];
  	Double_t C_s       = par[3];
  	Double_t gamma_s   = par[4];
  	Double_t inv_E_s   = par[5];
  	Double_t E_tr      = par[6];  
  	Double_t gamma_dif = par[7];
  	Double_t trans     = par[8]; //This parameter enables (1) or disables (0) the transition term
  	Double_t E3        = par[9]; //This parameter add (3) or not (0) a factor E^3 to the flux function
  	

  	diffusion = C_d*TMath::Power(((E + phi)/7.0),gamma_d);
  	source    = C_s*TMath::Power(((E + phi)/60.0),gamma_s)*TMath::Exp(-(E + phi)*inv_E_s);
  	flux      = (TMath::Power(E, 2)/TMath::Power(E + phi, 2))*(diffusion + source);
  	flux      = flux*(TMath::Power(TMath::Power(1+TMath::Power((E+phi)/E_tr,gamma_dif),-1),trans));
  	flux      = flux*TMath::Power(E,E3);
  
  	return flux;
}

Int_t  fitlines(){  
  
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(1111);
  gStyle->SetPalette(kSolar);

  TCanvas *c1= new TCanvas("c1","Electron fit",10,10,1000,600);
  //c1->SetGrid();
  c1->SetLogx();

  auto* file = new TFile("ssdc_canvas.root");
  TGraphAsymmErrors* graph = (TGraphAsymmErrors*)file->Get("graph1"); //"e+_AMS_PRL2019_ekin_000.txt"
    
  auto mg = new TMultiGraph("mg","mg");

  //TF1 *fit = new TF1("electron fit", Electron, 0.5, 1500.,7);
  //fit->SetParameters(0.85, 1.11e-02, -4.25, 3.9e-06, -3.11, -4.,-2.15);

  TF1 *fit = new TF1("electron fit", flux, 0.5, 1500.,10);
  fit->SetParameters(0.85, 1.11e-02, -4.25, 3.9e-06, -3.11,0., -4.,-2.15, 1, 0);
  //fit->SetParameters(1.11, 6.49e-02, -4.10, 6.75e-05, -2.61, 809.);
  fit->SetParNames("#varphi_{e^{-}}", "C_{a}", "#gamma_{a}", "C_{b}", "#gamma_{b}", "E_{t}", "#Delta #gamma_{t}");

  fit->SetParLimits(0, 0.1, 1.5);
  fit->SetParLimits(1, 0.5e-02, 10e-02);
  fit->SetParLimits(2, -6., -1.);
  fit->SetParLimits(3, 0.1e-06, 10e-05);
  fit->SetParLimits(4, -5., -1.);
  fit->SetParLimits(5, 1.,10.);
  fit->SetParLimits(6, -4.,-1.);

  graph->Fit("electron fit");
  gPad->Modified(); gPad->Update();
  graph->GetXaxis()->SetLimits(0.5, 2000.);
  gPad->Modified(); gPad->Update();
  


  TF1 *diff = new TF1("diff", Diffusion_e, 0.05,2000.,5);
  TF1 *sou = new TF1("sou", Source_e, 0.05,2000.,5);

  diff->SetParameter(0,fit->GetParameter(0));
  diff->SetParameter(1,fit->GetParameter(1));
  diff->SetParameter(2,fit->GetParameter(2));
  diff->SetParameter(3,fit->GetParameter(5));
  diff->SetParameter(4,fit->GetParameter(6));

  sou->SetParameter(0,fit->GetParameter(0));
  sou->SetParameter(1,fit->GetParameter(3));
  sou->SetParameter(2,fit->GetParameter(4));
  sou->SetParameter(3,fit->GetParameter(5));
  sou->SetParameter(4,fit->GetParameter(6));

  graph->Draw("AP plc pmc");
  diff->Draw("same LF plc");
  sou->Draw("same LF plc");

  c1->SaveAs("electron.pdf");
  
  return 0;
}  

  //TGraph *diff_gr = new TGraph(diff);
  //TGraph *sou_gr = new TGraph(sou);

  //mg->Add(diff_test, "pl");
  //mg->Add(sou_gr, "pl");
  //mg->Add(graph, "p");
  
  //mg->Draw("AP plc pmc");
  //gPad->Modified(); gPad->Update();
  //mg->GetXaxis()->SetLimits(0.5, 2000.);
  //mg->SetMinimum(0.);
  //mg->SetMaximum(28.);
  //gPad->Modified(); gPad->Update();
  //mg->Draw("AP plc pmc"); 
  //mg->GetXaxis()->SetLimits(1,3000.);
  //mg->GetYaxis()->SetRangeUser(0.,30.);
 
