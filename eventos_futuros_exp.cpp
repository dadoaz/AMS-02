/////////////////////////////////////////////////////////////////////////
//
// Estimación de la anisotropía dipolar e+ e- esperada en los experimentos
// calorimétricos CALET, DAMPE y HERD. Modificación del código previo del cálculo de AMS
//
// Daniel Doménech Azorín dadoa@alumni.u.,es
//
// 26.04.2021
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


Double_t Positron(Double_t *x, Double_t *par);
Double_t Electron(Double_t *x, Double_t *par);
Double_t Diffusion_p(Double_t *x, Double_t *par);
Double_t Diffusion_e(Double_t *x, Double_t *par);
Double_t Source_p(Double_t *x, Double_t *par);
Double_t Source_e(Double_t *x, Double_t *par);
void     fill_hist_events(TH1 *hist, TF1 *flux, bool CONST);
void     fill_hist_integral(TH1 *integral, TH1* hist);
void     purity(TH1 *pur, TH1 *source, TH1 *diff);
void     Neff(TH1 *eff, TH1 *sou, TH1 *diff, TH1 *pur);
void     sens(TH1 *sens, TH1 *Neff);
void     print_hist(TH1 *hist, TString name);
Double_t Diffusion_p_E3(Double_t *x, Double_t *par);
Double_t Source_p_E3(Double_t *x, Double_t *par);
Double_t Diffusion_e_E3(Double_t *x, Double_t *par);
Double_t Source_e_E3(Double_t *x, Double_t *par);
void     purity_calorimeters(TH1 *pur, TH1 *p_sou, TH1 *p_dif, TH1 *e_sou, TH1 *e_dif);
void     fill_hist_events_calorimeters(TH1 *hist, TF1 *flux, TString exp );

Int_t eventos_futuros_exp(){

	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1111);
	//gStyle->SetPalette(kCool); //kCool, kCopper, kTemperatureMap, kDeepSea

    //Let's introduce the functions with AMS parameters

	TF1 *se_flux = new TF1("se_flux", Source_e, 0.5, 1500.,5);
    se_flux->SetParameters(0.87, 3.96e-06, -3.14, 3.94, -2.14);
    se_flux->SetParNames("#varphi_{e^{-}}", "C_{a}", "#gamma_{a}", "E_{t}", "#Delta #gamma_{t}");

    TF1 *de_flux = new TF1("de_flux", Diffusion_e, 0.5, 1500.,5);
    de_flux->SetParameters(0.87, 1.13e-02, -4.31, 3.94, -2.14);
    de_flux->SetParNames("#varphi_{e^{-}}", "C_{b}", "#gamma_{b}", "E_{t}", "#Delta #gamma_{t}");
    
    TF1 *sp_flux = new TF1("sp_flux", Source_p, 0.5, 1500.,4);
    sp_flux->SetParameters(1.10, 6.80e-05, -2.58, 1.23e-3);
    sp_flux->SetParNames("#varphi_{e^{-}}", "C_{a}", "#gamma_{a}", "E_{t}", "#Delta #gamma_{t}");
        
    TF1 *dp_flux = new TF1("dp_flux", Diffusion_p, 0.5, 1500.,3);
    dp_flux->SetParameters(1.10, 6.51e-02, -4.07);
    dp_flux->SetParNames("#varphi_{e^{-}}", "C_{b}", "#gamma_{b}", "E_{t}", "#Delta #gamma_{t}");


    //Let's introduce the histograms
    Int_t    N    = 15000;
    Double_t min  = 1;
    Double_t max  = 1500;
    
    // Create histograms for each source
    TH1D *h_p_dif_calet = new TH1D("Diffusion_p_calet", "# expected positrons @ Calet; E [GeV]; Events/0.1GeV", N, min, max);
    TH1D *h_p_sou_calet = new TH1D("Source_p_calet"   , "# expected positrons @ Calet; E [GeV]; Events/0.1GeV", N, min, max);
    TH1D *h_e_dif_calet = new TH1D("Diffusion_e_calet", "# expected electrons @ Calet; E [GeV]; Events/0.1GeV", N, min, max);
    TH1D *h_e_sou_calet = new TH1D("Source_e_calet"   , "# expected electrons @ Calet; E [GeV]; Events/0.1GeV", N, min, max);

    TH1D *h_p_dif_dampe = new TH1D("Diffusion_p_dampe", "# expected positrons @ Dampe; E [GeV]; Events/0.1GeV", N, min, max);
    TH1D *h_p_sou_dampe = new TH1D("Source_p_dampe"   , "# expected positrons @ Dampe; E [GeV]; Events/0.1GeV", N, min, max);
    TH1D *h_e_dif_dampe = new TH1D("Diffusion_e_dampe", "# expected electrons @ Dampe; E [GeV]; Events/0.1GeV", N, min, max);
    TH1D *h_e_sou_dampe = new TH1D("Source_e_dampe"   , "# expected electrons @ Dampe; E [GeV]; Events/0.1GeV", N, min, max);

    TH1D *h_p_dif_herd = new TH1D("Diffusion_p_herd", "# expected positrons @ HERD; E [GeV]; Events/0.1GeV", N, min, max);
    TH1D *h_p_sou_herd = new TH1D("Source_p_herd"   , "# expected positrons @ HERD; E [GeV]; Events/0.1GeV", N, min, max);
    TH1D *h_e_dif_herd = new TH1D("Diffusion_e_herd", "# expected electrons @ HERD; E [GeV]; Events/0.1GeV", N, min, max);
    TH1D *h_e_sou_herd = new TH1D("Source_e_herd"   , "# expected electrons @ HERD; E [GeV]; Events/0.1GeV", N, min, max);

    //Fill histograms with number of expected events

    fill_hist_events_calorimeters(h_p_dif_calet, dp_flux, "Calet");
    fill_hist_events_calorimeters(h_p_sou_calet, sp_flux, "Calet");
    fill_hist_events_calorimeters(h_e_dif_calet, de_flux, "Calet");
    fill_hist_events_calorimeters(h_e_sou_calet, se_flux, "Calet");

    fill_hist_events_calorimeters(h_p_dif_dampe, dp_flux, "Dampe");
    fill_hist_events_calorimeters(h_p_sou_dampe, sp_flux, "Dampe");
    fill_hist_events_calorimeters(h_e_dif_dampe, de_flux, "Dampe");
    fill_hist_events_calorimeters(h_e_sou_dampe, se_flux, "Dampe");

    fill_hist_events_calorimeters(h_p_dif_herd, dp_flux, "Herd");
    fill_hist_events_calorimeters(h_p_sou_herd, sp_flux, "Herd");
    fill_hist_events_calorimeters(h_e_dif_herd, de_flux, "Herd");
    fill_hist_events_calorimeters(h_e_sou_herd, se_flux, "Herd");

    // Merge all background events in one simple histogram

    TH1D *h_bkg_calet = new TH1D("Bkg_ep_calet", "# expected isotropic e^{-} + e^{+}; E[GeV]; Events/0.1Gev", N, min, max);
    TList *l_calet = new TList;
    l_calet->Add(h_e_dif_calet);
    l_calet->Add(h_e_sou_calet);
    l_calet->Add(h_p_dif_calet);

    h_bkg_calet->Merge(l_calet);

    TH1D *h_bkg_dampe = new TH1D("Bkg_ep_dampe", "# expected isotropic e^{-} + e^{+}; E[GeV]; Events/0.1Gev", N, min, max);
    TList *l_dampe = new TList;
    l_dampe->Add(h_e_dif_dampe);
    l_dampe->Add(h_e_sou_dampe);
    l_dampe->Add(h_p_dif_dampe);

    h_bkg_dampe->Merge(l_dampe);

    TH1D *h_bkg_herd = new TH1D("Bkg_ep_herd", "# expected isotropic e^{-} + e^{+}; E[GeV]; Events/0.1Gev", N, min, max);
    TList *l_herd = new TList;
    l_herd->Add(h_e_dif_herd);
    l_herd->Add(h_e_sou_herd);
    l_herd->Add(h_p_dif_herd);

    h_bkg_herd->Merge(l_herd);
    
    // b) Now let's obtain the number of observed CR  above some energy
    
    TH1D *h_bkg_int_calet = new TH1D("Bkg_ep_int_calet", "Survival isotropic e^{+} + e^{-} @ Calet; Minimum E [GeV]; Events", N, min, max);
    TH1D *h_p_sou_int_calet = new TH1D("Source_p_int_calet"   , "Survival source e^{-} @ Calet; Minimum E [GeV]; Events", N, min, max);

    TH1D *h_bkg_int_dampe = new TH1D("Bkg_ep_int_dampe", "Survival isotropic e^{+} + e^{-} @ Dampe; Minimum E [GeV]; Events", N, min, max);
    TH1D *h_p_sou_int_dampe = new TH1D("Source_p_int_dampe"   , "Survival source e^{-} @ Dampe; Minimum E [GeV]; Events", N, min, max);

    TH1D *h_bkg_int_herd = new TH1D("Bkg_ep_int_herd", "Survival isotropic e^{+} + e^{-} @ Herd; Minimum E [GeV]; Events", N, min, max);
    TH1D *h_p_sou_int_herd = new TH1D("Source_p_int_herd"   , "Survival source e^{-} @ Herd; Minimum E [GeV]; Events", N, min, max);

    fill_hist_integral(h_bkg_int_calet, h_bkg_calet);
    fill_hist_integral(h_p_sou_int_calet, h_p_sou_calet);

    fill_hist_integral(h_bkg_int_dampe, h_bkg_dampe);
    fill_hist_integral(h_p_sou_int_dampe, h_p_sou_dampe);

    fill_hist_integral(h_bkg_int_herd, h_bkg_herd);
    fill_hist_integral(h_p_sou_int_herd, h_p_sou_herd);


    //c) Let's calculate the purity of the sample
    
    //Now let's obtain the number of observed CR  above some energy
    TH1D *h_p_pur_calet = new TH1D("Purity_calet", "positron purity ; E [GeV]; p(E)", N, min, max);
    TH1D *h_p_pur_int_calet = new TH1D("Purity_integrated_calet"   , "positron integrated purity; E_{min} [GeV]; p(E)", N, min, max);

    purity(h_p_pur_calet, h_p_sou_calet, h_bkg_calet);
    purity(h_p_pur_int_calet, h_p_sou_int_calet, h_bkg_int_calet);

    TH1D *h_p_pur_dampe = new TH1D("Purity_dampe", "positron purity ; E [GeV]; p(E)", N, min, max);
    TH1D *h_p_pur_int_dampe = new TH1D("Purity_integrated_dampe"   , "positron integrated purity; E_{min} [GeV]; p(E)", N, min, max);

    purity(h_p_pur_dampe, h_p_sou_dampe, h_bkg_dampe);
    purity(h_p_pur_int_dampe, h_p_sou_int_dampe, h_bkg_int_dampe);

    TH1D *h_p_pur_herd = new TH1D("Purity_herd", "positron purity ; E [GeV]; p(E)", N, min, max);
    TH1D *h_p_pur_int_herd = new TH1D("Purity_integrated_herd"   , "positron integrated purity; E_{min} [GeV]; p(E)", N, min, max);

    purity(h_p_pur_herd, h_p_sou_herd, h_bkg_herd);
    purity(h_p_pur_int_herd, h_p_sou_int_herd, h_bkg_int_herd);
    

    //Now let's obtain the number of effective positrons
    TH1D *h_p_Neff_calet = new TH1D("Neff_calet", "# eff positrons ; E [GeV]; Events per 0.1 GeV", N, min, max);
    TH1D *h_p_Neff_int_calet = new TH1D("Neff_int_calet"   , "# integrated eff positrons; E_{min} [GeV]; Events", N, min, max);

    Neff(h_p_Neff_calet, h_p_sou_calet, h_bkg_calet, h_p_pur_calet);
    fill_hist_integral(h_p_Neff_int_calet, h_p_Neff_calet);

    TH1D *h_p_Neff_dampe = new TH1D("Neff_dampe", "# eff positrons ; E [GeV]; Events per 0.1 GeV", N, min, max);
    TH1D *h_p_Neff_int_dampe = new TH1D("Neff_int_dampe"   , "# integrated eff positrons; E_{min} [GeV]; Events", N, min, max);

    Neff(h_p_Neff_dampe, h_p_sou_dampe, h_bkg_dampe, h_p_pur_dampe);
    fill_hist_integral(h_p_Neff_int_dampe, h_p_Neff_dampe);

    TH1D *h_p_Neff_herd = new TH1D("Neff_herd", "# eff positrons ; E [GeV]; Events per 0.1 GeV", N, min, max);
    TH1D *h_p_Neff_int_herd = new TH1D("Neff_int_herd"   , "# integrated eff positrons; E_{min} [GeV]; Events", N, min, max);

    Neff(h_p_Neff_herd, h_p_sou_herd, h_bkg_herd, h_p_pur_herd);
    fill_hist_integral(h_p_Neff_int_herd, h_p_Neff_herd);

    
    // Calculemos ahora la sensibilidad esperada

    TH1D* h_sens_calet = new TH1D("sensitivity_calet", "Expected anisotropy sensibility; E_{min} [GeV]; #delta_{M}^{exp}", N, min, max);

    sens(h_sens_calet, h_p_Neff_int_calet);

    TH1D* h_sens_dampe = new TH1D("sensitivity_dampe", "Expected anisotropy sensibility; E_{min} [GeV]; #delta_{M}^{exp}", N, min, max);

    sens(h_sens_dampe, h_p_Neff_int_dampe);
    
    TH1D* h_sens_herd = new TH1D("sensitivity_herd", "Expected anisotropy sensibility; E_{min} [GeV]; #delta_{M}^{exp}", N, min, max);

    sens(h_sens_herd, h_p_Neff_int_herd);


    //Incluyamos en el .root los flujos escalados a E^3

    TF1 *sp_flux_E3 = new TF1("sp_flux_E3", Source_p_E3, 0.5, 1500.,4);
    sp_flux_E3->SetParameters(1.10, 6.80e-05, -2.58, 1.23e-3);
    sp_flux_E3->SetParNames("#varphi_{e^{-}}", "C_{a}", "#gamma_{a}", "E_{t}", "#Delta #gamma_{t}");
        
    TF1 *dp_flux_E3 = new TF1("dp_flux_E3", Diffusion_p_E3, 0.5, 1500.,3);
    dp_flux_E3->SetParameters(1.10, 6.51e-02, -4.07);
    dp_flux_E3->SetParNames("#varphi_{e^{-}}", "C_{b}", "#gamma_{b}", "E_{t}", "#Delta #gamma_{t}");

    TF1 *se_flux_E3 = new TF1("se_flux_E3", Source_e_E3, 0.5, 1500.,5);
    se_flux_E3->SetParameters(0.87, 3.96e-06, -3.14, 3.94, -2.14);
    se_flux_E3->SetParNames("#varphi_{e^{-}}", "C_{a}", "#gamma_{a}", "E_{t}", "#Delta #gamma_{t}");

    TF1 *de_flux_E3 = new TF1("de_flux_E3", Diffusion_e_E3, 0.5, 1500.,5);
    de_flux_E3->SetParameters(0.87, 1.13e-02, -4.31, 3.94, -2.14);
    de_flux_E3->SetParNames("#varphi_{e^{-}}", "C_{b}", "#gamma_{b}", "E_{t}", "#Delta #gamma_{t}");


    // Hagamos print de algunas magnitudes

    print_hist(h_p_sou_int_calet, "events B source");
    print_hist(h_bkg_int_calet, "events B diff");
    print_hist(h_p_Neff_int_calet, "Neff int");
    print_hist(h_sens_calet, "Sensitivity");


    // Save all histograms in a .root file -> to drawing graphs faster in draw_exp.cpp

    TFile f("Neff_exp.root", "new");
    h_p_Neff_calet->Write();
    h_p_Neff_int_calet->Write();
    h_p_pur_calet->Write();
    h_p_pur_int_calet->Write();
    h_p_sou_calet->Write();
    h_p_sou_int_calet->Write();
    h_bkg_calet->Write();
    h_bkg_int_calet->Write();
    h_sens_calet->Write();
    h_p_dif_calet->Write();
    h_e_sou_calet->Write();
    h_e_dif_calet->Write();

    h_p_Neff_dampe->Write();
    h_p_Neff_int_dampe->Write();
    h_p_pur_dampe->Write();
    h_p_pur_int_dampe->Write();
    h_p_sou_dampe->Write();
    h_p_sou_int_dampe->Write();
    h_bkg_dampe->Write();
    h_bkg_int_dampe->Write();
    h_sens_dampe->Write();
    h_p_dif_dampe->Write();
    h_e_sou_dampe->Write();
    h_e_dif_dampe->Write();

    h_p_Neff_herd->Write();
    h_p_Neff_int_herd->Write();
    h_p_pur_herd->Write();
    h_p_pur_int_herd->Write();
    h_p_sou_herd->Write();
    h_p_sou_int_herd->Write();
    h_bkg_herd->Write();
    h_bkg_int_herd->Write();
    h_sens_herd->Write();
    h_p_dif_herd->Write();
    h_e_sou_herd->Write();
    h_e_dif_herd->Write();

    de_flux_E3->Write();
    se_flux_E3->Write();
    dp_flux_E3->Write();
    sp_flux_E3->Write();


    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

    return (TMath::Power(E, 2)/TMath::Power(E + phi, 2))*diff; 

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
    
    return (TMath::Power(E, 2)/TMath::Power(E + phi, 2))*TMath::Power(1+TMath::Power((E+phi)/E_tr,gamma_dif),-1)*diff; 
  
}

Double_t Source_p(Double_t *x, Double_t *par){
  
    Double_t sou;
    Double_t E       = x[0];
    Double_t phi     = par[0];
    Double_t C_s     = par[1];
    Double_t gamma_s = par[2];
    Double_t inv_E_s = par[3];
  
    sou = C_s*TMath::Power(((E + phi)/60.0),gamma_s)*TMath::Exp(-(E + phi)*inv_E_s);
    
    return (TMath::Power(E, 2)/TMath::Power(E + phi, 2))*sou;

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

    return (TMath::Power(E, 2)/TMath::Power(E + phi, 2))*TMath::Power(1+TMath::Power((E+phi)/E_tr,gamma_dif),-1)*sou;

}

void fill_hist_events(TH1 *hist, TF1 *flux, bool CONST){
    
    //Let's read the acceptance and time exposition functions in AMS

    //Acceptance
    auto* file_acc = new TFile("fAcc.root");
    TF1* acc = (TF1*)file_acc->Get("fAcc"); 
	
    //Time exposition
    auto* file_exp = new TFile("fExpo.root");
    TF1* exp = (TF1*)file_exp->Get("fExpo");

    Double_t xx;
    Double_t fx;
    
    //Now let's fill the histograms with the expected number of events N = \phi * T * A * \Delta E
    //Define the case of evaluation
    // CONST == true -> acceptance and time exposition constant  a) constant. Caution! We have converted acceptance to cm^2 to m^2 by *1e-4
    // CONST == false  -> acceptance and time eposition with energy dependence  b) scalated to a total exposition of 17 years (mission lifetime)

    Int_t N = hist->GetNbinsX();

    for(Int_t i = 0; i < N; i++){
        
        xx = hist->GetBinCenter(i);
        
        if(CONST){
          fx = flux->Eval(xx)*0.04*3.86e8*hist->GetBinWidth(i);
        }
        else{
          fx = flux->Eval(xx)*acc->Eval(xx)*exp->Eval(xx)*hist->GetBinWidth(i)*1e-4*17./9;
        }
        
        hist->SetBinContent(i, fx);
    }
}

void fill_hist_integral(TH1 *integral, TH1* hist){

    Int_t N = integral->GetNbinsX();

    if( N == hist->GetNbinsX() ){

      for(Int_t i = 0; i < N; i++){
		    integral->SetBinContent(i, hist->Integral(i,N));
      }
    }
}

void purity(TH1 *pur, TH1 *source, TH1 *diff){

    Int_t N = pur->GetNbinsX();

    if(( N == source->GetNbinsX()) && ( N == diff->GetNbinsX() )){
      
      for(Int_t i = 0; i < N; i++){
    	
        if(source->GetBinContent(i) + diff->GetBinContent(i) != 0){
          Double_t xx = source->GetBinContent(i)/(source->GetBinContent(i) + diff->GetBinContent(i));
          pur->SetBinContent(i, xx);
        }
      }
    } 
}

void Neff(TH1 *eff, TH1 *sou, TH1 *diff, TH1 *pur){

    Int_t N = eff->GetNbinsX();

    if ( (N == pur->GetNbinsX()) && (N== sou->GetNbinsX()) && ( N == diff->GetNbinsX())) {
      
        for(Int_t i = 0; i < N; i++){

            Double_t xsou   = sou->GetBinContent(i);
            Double_t xdiff  = diff->GetBinContent(i);
            Double_t xpur   = pur->GetBinContent(i);

            eff->SetBinContent(i, (xsou + xdiff) * xpur * xpur);
        }
    }
}

void sens(TH1 *sens, TH1 *Neff){

    Double_t N = sens->GetNbinsX();

    if( N == Neff->GetNbinsX()){
        
        for(Int_t i = 0; i < N; i++){
            
            Double_t xx = TMath::Sqrt( TMath::ChisquareQuantile(1./2, 3) * 3./ Neff->GetBinContent(i) );
            sens->SetBinContent(i, xx);
        }
    }
}

void print_hist(TH1 *hist, TString name){
  
  cout << "====================" << endl;
  cout << name << endl;
  cout << "\t" << hist->GetBinContent(1) << "\t" << hist->GetBinCenter(1) << endl;
  cout << "\t" << hist->GetBinContent(100) << "\t" << hist->GetBinCenter(100) << endl;
  cout << "\t" << hist->GetBinContent(200) << "\t" << hist->GetBinCenter(200) << endl;
  cout << "\t" << hist->GetBinContent(500) << "\t" << hist->GetBinCenter(500) << endl;
  cout << "\t" << hist->GetBinContent(1000) << "\t" << hist->GetBinCenter(1000) << endl;
  cout << "\t" << hist->GetBinContent(2000) << "\t" << hist->GetBinCenter(2000) << endl;
  cout << "\t" << hist->GetBinContent(5000) << "\t" << hist->GetBinCenter(5000) << endl;
  cout << "\t" << hist->GetBinContent(10000) << "\t" << hist->GetBinCenter(10000) << endl;
}

Double_t Diffusion_p_E3(Double_t *x, Double_t *par){
  
    Double_t diff;
    Double_t E       = x[0];
    Double_t phi     = par[0];
    Double_t C_d     = par[1];
    Double_t gamma_d = par[2];

    diff = C_d*TMath::Power(((E + phi)/7.0),gamma_d);

    return (TMath::Power(E, 5)/TMath::Power(E + phi, 2))*diff; 

}

Double_t Source_p_E3(Double_t *x, Double_t *par){
  
    Double_t sou;
    Double_t E       = x[0];
    Double_t phi     = par[0];
    Double_t C_s     = par[1];
    Double_t gamma_s = par[2];
    Double_t inv_E_s = par[3];
  
    sou = C_s*TMath::Power(((E + phi)/60.0),gamma_s)*TMath::Exp(-(E + phi)*inv_E_s);
    
    return (TMath::Power(E, 5)/TMath::Power(E + phi, 2))*sou;

}

Double_t Diffusion_e_E3(Double_t *x, Double_t *par){
    
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

Double_t Source_e_E3(Double_t *x, Double_t *par){
  
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

void purity_calorimeters(TH1 *pur, TH1 *p_sou, TH1 *p_dif, TH1 *e_sou, TH1 *e_dif){

    Int_t N = pur->GetNbinsX();

    if(( N == p_sou->GetNbinsX()) && ( N == p_dif->GetNbinsX() ) && ( N == e_sou->GetNbinsX()) && ( N == e_dif->GetNbinsX() )){
      
      for(Int_t i = 0; i < N; i++){
    	
        if(p_sou->GetBinContent(i) + p_dif->GetBinContent(i) + e_sou->GetBinContent(i) + e_dif->GetBinContent(i) != 0){
          Double_t xx = p_sou->GetBinContent(i)/(p_sou->GetBinContent(i) + p_dif->GetBinContent(i) + e_sou->GetBinContent(i) + e_dif->GetBinContent(i));
          pur->SetBinContent(i, xx);
        }
      }
    } 
}

void fill_hist_events_calorimeters(TH1 *hist, TF1 *flux, TString exp ){

    Double_t xx;
    Double_t fx;
    Double_t eff_year = 365*24*3600*0.9; /*Effective years at mission in seconds (0.9 livetime)*/
    
    //Switch cannot work with strings, let's convert it into integeers
    Int_t t = 0;
    if( exp == "Calet"){t = 1;}
    if( exp == "Dampe"){t = 2;}
    if( exp == "Herd") {t = 3;}

    //Now let's fill the histograms with the expected number of events N = \phi * T * A * \Delta E
    
    Int_t N = hist->GetNbinsX();
        
    for(Int_t i = 0; i < N; i++){
    
        xx = hist->GetBinCenter(i);
        
        switch (t){

            case 1:
                fx = flux->Eval(xx)*0.1*6*eff_year*hist->GetBinWidth(i);
                break;
            case 2:
                fx = flux->Eval(xx)*0.3*6*eff_year*hist->GetBinWidth(i);
                break;

            case 3:
                fx = flux->Eval(xx)*3*10*eff_year*hist->GetBinWidth(i);
                break;
            default:
                cout << "Tonto";
                break;
        }
        
        hist->SetBinContent(i, fx);
    }


    
}

