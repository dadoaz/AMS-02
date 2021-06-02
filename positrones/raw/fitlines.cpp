/////////////////////////////////////////////////////////////////////////
//
// Ejemplo sencillo de ajuste en ROOT por minimos cuadrados 
// (opcion por defecto) 
// Ajuste a una linea recta
//
// Juan Zu�iga zuniga@ific.uv.es 
//
// 21.01.2013
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

Double_t Diffusion(Double_t *x, Double_t *par){
  
  Double_t diff;
  Double_t E       = x[0];
  Double_t phi     = par[0];
  Double_t C_d     = par[1];
  Double_t gamma_d = par[2];

  diff = C_d*TMath::Power(((E + phi)/7.0),gamma_d);
  
  return (TMath::Power(E, 2)/TMath::Power(E + phi, 2))*diff; 

}

Double_t Source(Double_t *x, Double_t *par){
  
  Double_t sou;
  Double_t E       = x[0];
  Double_t phi     = par[0];
  Double_t C_s     = par[1];
  Double_t gamma_s = par[2];
  Double_t inv_E_s = par[3];

  sou = C_s*TMath::Power(((E + phi)/60.0),gamma_s)*TMath::Exp(-(E + phi)*inv_E_s);
  
  return (TMath::Power(E, 2)/TMath::Power(E + phi, 2))*sou;

}

Int_t  fitlines(){  
  
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(1111);
  gStyle->SetPalette(kSolar);

  TCanvas *c1= new TCanvas("c1","Positron fit",10,10,900,600);
  c1->SetGrid();
  c1->SetLogx();

  auto* file = new TFile("ssdc_canvas.root");
  TGraphAsymmErrors* graph = (TGraphAsymmErrors*)file->Get("graph1"); //"e+_AMS_PRL2019_ekin_000.txt"
  
  TF1 *fit = new TF1("positron fit", Positron, 0., 900., 6);
  fit->SetParameters(1.10, 6.51e-02, -4.07, 6.8e-05, -2.58, 1.23e-3);
  fit->SetParNames("#varphi_{e^{+}}", "C_{d}", "#gamma_{d}", "C_{s}", "#gamma_{s}", "1/E_{s}");

  fit->SetParLimits(0, 0.1, 2.);
  fit->SetParLimits(1, 1e-02, 10e-02);
  fit->SetParLimits(2, -5., -3.);
  fit->SetParLimits(3, 1e-05, 10e-05);
  fit->SetParLimits(4, -4, -1);
  fit->SetParLimits(5, 0.5e-3, 2e-3);

  gPad->Modified(); gPad->Update();
  graph->GetXaxis()->SetLimits(0.5, 2000.);
  gPad->Modified(); gPad->Update();
  graph->Fit("positron fit");

  TF1 *diff = new TF1("diff", Diffusion, 0.05,2000.,3);
  TF1 *sou = new TF1("sou", Source, 0.05,2000.,4);

  diff->SetParameter(0,fit->GetParameter(0));
  diff->SetParameter(1,fit->GetParameter(1));
  diff->SetParameter(2,fit->GetParameter(2));

  sou->SetParameter(0,fit->GetParameter(0));
  sou->SetParameter(1,fit->GetParameter(3));
  sou->SetParameter(2,fit->GetParameter(4));
  sou->SetParameter(3,fit->GetParameter(5));

  graph->Draw("AP plc pmc");
  diff->Draw("same LF plc");
  sou->Draw("same LF plc");

  c1->SaveAs("positron_flux_source+diffusion[debug].pdf");

  return 0;
}  
