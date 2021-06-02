/////////////////////////////////////////////////////////////////////////
//
// Comprobación de que el ajuste del flujo de electrones medidos por AMS
// ofrece el mismo resultados con los parámetros de ajuste obtenidos en AMS
// 
//
// Daniel Doménech Azorín dadoa@alumni.uv,es
//
// 30.03.2021
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

  Double_t diffusion, source, flux, trans;
  Double_t E         = x[0];
  Double_t phi       = par[0];
  Double_t C_d       = par[1];
  Double_t gamma_d   = par[2];
  Double_t C_s       = par[3];
  Double_t gamma_s   = par[4];
  Double_t E_tr      = par[5];  
  Double_t gamma_dif = par[6];

  //Written this way to prove if operation ordering was the reason of failure
  diffusion = (E + phi)/20.0;
  diffusion = TMath::Power(diffusion, gamma_d);
  diffusion = C_d*diffusion;

  source    = (E + phi)/300.0;
  source    = TMath::Power(source, gamma_s);
  source    = C_s*source;
  
  trans     = (E+phi)/E_tr;
  trans     = TMath::Power(trans, gamma_dif);
  trans     = 1/(1+trans);

  flux      = E/(E+phi);
  flux      = flux * flux;
  flux      = E*E*E * flux;
  flux      = flux*trans;
  flux      = flux*(source + diffusion);
  
  return flux;
}

Double_t chiAMS(){

    //Read the data from .root file
    auto* file = new TFile("ssdc_canvas.root");
    TGraphAsymmErrors* graph = (TGraphAsymmErrors*)file->Get("graph1");

    //Let's define the statistic
    Double_t chi  = 0;
    Double_t temp = 0;
    Double_t E    = 0;
    Double_t teo  = 0;
    //Let's introduce the parameters {phi, C_a, gamma_a, C_b, gamma_b, E_tr, gamma_dif}
    Double_t par[7] = {0.87, 1.13e-2, -4.31, 3.96e-6, -3.14, 3.94, -2.14};

    for(Int_t i = 0; i < graph->GetN(); i++){
        E = graph->GetPointX(i);
        teo = Electron(&E, par);

        //cout << E << " " << teo << endl;

        temp = TMath::Power((graph->GetPointY(i) - teo)/(1*graph->GetErrorY(i)),2);
        chi += temp;
        teo = 0;
    }

    return chi;
}