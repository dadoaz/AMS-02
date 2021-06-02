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


Double_t Positron(Double_t *x, Double_t *par);
Double_t Electron(Double_t *x, Double_t *par);
Double_t Diffusion_p(Double_t *x, Double_t *par);
Double_t Diffusion_e(Double_t *x, Double_t *par);
Double_t Source_p(Double_t *x, Double_t *par);
Double_t Source_e(Double_t *x, Double_t *par);
void     fill_hist_events(TH1 *hist, TF1 *flux, bool CONST);
void     fill_hist_integral(TH1 *integral, TH1* hist);
void     purity(TH1 *pur, TH1 *source, TH1 *diff);
void     Ntot(TH1 *tot, TH1 *sou, TH1 *diff);
void     Neff(TH1 *eff, TH1 *sou, TH1 *diff, TH1 *pur);
void     sens(TH1 *sens, TH1 *Neff);
void     print_hist(TH1 *hist, TString name);
Double_t Diffusion_p_E3(Double_t *x, Double_t *par);
Double_t Source_p_E3(Double_t *x, Double_t *par);


Int_t eventos(){

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
    TH1D *h_p_dif = new TH1D("Diffusion_p", "# expected positrons; E [GeV]; Events/0.1GeV", N, min, max);
    TH1D *h_p_sou = new TH1D("Source_p"   , "# expected positrons; E [GeV]; Events/0.1GeV", N, min, max);

    TH1D *h_e_dif = new TH1D("Diffusion_e", "# expected electrons; E [GeV]; Events/0.1GeV", N, min, max);
    TH1D *h_e_sou = new TH1D("Source_e"   , "# expected electrons; E [GeV]; Events/0.1GeV", N, min, max);

    //Fill histograms with number of expected events
    bool CONST = false;
    fill_hist_events(h_p_dif, dp_flux, CONST);
    fill_hist_events(h_p_sou, sp_flux, CONST);

    fill_hist_events(h_e_dif, de_flux, CONST);
    fill_hist_events(h_e_sou, se_flux, CONST);
    
    TCanvas *c1= new TCanvas("c1","Electron fit",10,10,900,600);
	c1->SetGrid();
	c1->SetLogx();
    c1->SetLogy();

    //Create a legend
    TLegend* leg = new TLegend(0.6, 0.8, 0.9, 0.9);
    leg->AddEntry(h_p_dif, "Diffusion");
    leg->AddEntry(h_p_sou, "Source");

    h_p_dif->SetMinimum(1);
	h_p_dif->SetFillStyle(3144);
    h_p_dif->SetFillColor(kGray+3);
    
    h_p_sou->SetFillStyle(3001);
    h_p_sou->SetFillColor(kOrange-3);

    h_p_dif->Draw("");
    h_p_sou->Draw("same");
    leg->Draw("same");
    

	  //Now let's obtain the number of observed CR  above some energy
    TH1D *h_p_dif_int = new TH1D("Diffusion_p_int", "Survival positrons; Minimum E [GeV]; Events", N, min, max);
    TH1D *h_p_sou_int = new TH1D("Source_p_int"   , "Survival positrons; Minimum E [GeV]; Events", N, min, max);

    fill_hist_integral(h_p_dif_int, h_p_dif);
    fill_hist_integral(h_p_sou_int, h_p_sou);

    TCanvas *c2= new TCanvas("c2","Electron fit",10,10,900,600);
    c2->SetGrid();
    c2->SetLogx();
    c2->SetLogy();
    
    h_p_dif_int->SetMinimum(1);
	h_p_dif_int->SetFillStyle(3144);
    h_p_dif_int->SetFillColor(kGray+3);

    h_p_sou_int->SetFillStyle(3001);
    h_p_sou_int->SetFillColor(kOrange-3);

    h_p_dif_int->Draw();
    h_p_sou_int->Draw("same");
    leg->Draw("same");



    //c) Let's calculate the purity of the sample
    
    //Now let's obtain the number of observed CR  above some energy
    TH1D *h_p_pur = new TH1D("Purity", "positron purity ; E [GeV]; p(E)", N, min, max);
    TH1D *h_p_pur_int = new TH1D("Purity integrated"   , "positron integrated purity; E_{min} [GeV]; p(E)", N, min, max);

    purity(h_p_pur, h_p_sou, h_p_dif);
    purity(h_p_pur_int, h_p_sou_int, h_p_dif_int);

    TCanvas *c3= new TCanvas("c3","Electron fit",10,10,900,600);
    c3->SetGrid();
    c3->SetLogx();
    //c3->SetLogy();

    h_p_pur->SetFillColor(kCyan-7);
    h_p_pur->SetFillStyle(3001);
    h_p_pur->Draw();//"pfc plc"
    
    
    TCanvas *c7= new TCanvas("c7","Electron fit",10,10,900,600);
    c7->SetGrid();
    c7->SetLogx();
    //c7->SetLogy();

    h_p_pur_int->SetFillColor(kMagenta-7);
    h_p_pur_int->SetFillStyle(3001);
    h_p_pur_int->Draw();

    TString d = ""; 
    //c1->SaveAs("positron_expected_events_Benergydep" + d + ".pdf");    //Benergydep
    //c2->SaveAs("positron_expected_survival_Benergydep" + d + ".pdf");  
    //c3->SaveAs("positron_purity_B" + d + ".pdf");    //Aconst
    c7->SaveAs("positron_purity_integrated_B" + d + ".pdf");


    //Now let's obtain the number of total positrons
    TH1D *h_p_Ntot = new TH1D("Ntot", "# tot positrons ; E [GeV]; Events per 0.1 GeV", N, min, max);
    TH1D *h_p_Ntot_int = new TH1D("Ntot int"   , "# integrated tot positrons; E_{min} [GeV]; Events", N, min, max);

    Ntot(h_p_Ntot, h_p_sou, h_p_dif);
    fill_hist_integral(h_p_Ntot_int, h_p_Ntot);


    //Now let's obtain the number of effective positrons
    TH1D *h_p_Neff = new TH1D("Neff", "# eff positrons ; E [GeV]; Events per 0.1 GeV", N, min, max);
    TH1D *h_p_Neff_int = new TH1D("Neff int"   , "# integrated eff positrons; E_{min} [GeV]; Events", N, min, max);

    Neff(h_p_Neff, h_p_sou, h_p_dif, h_p_pur);
    fill_hist_integral(h_p_Neff_int, h_p_Neff);



    TCanvas *c4= new TCanvas("c4","Electron fit",10,10,900,600);
    c4->SetGrid();
    c4->SetLogx();
    c4->SetLogy();

    h_p_Neff->SetMinimum(1);
    h_p_Neff->SetFillColor(kCyan-7);
    h_p_Neff->SetFillStyle(3001);
    h_p_Neff->Draw();
    
    TCanvas *c5= new TCanvas("c5","Electron fit",10,10,900,600);
    c5->SetGrid();
    c5->SetLogx();
    c5->SetLogy();
    
    h_p_Neff_int->SetMinimum(1);
    h_p_Neff_int->SetFillColor(kMagenta-7);
    h_p_Neff_int->SetFillStyle(3001);
    h_p_Neff_int->Draw();

    
    c4->SaveAs("positron_Neff_B" + d + ".pdf");    //Aconst
    c5->SaveAs("positron_Neff_int_B" + d + ".pdf");
	  

    // Calculemos ahora la sensibilidad esperada

    TH1D* h_sens = new TH1D("sensitivity", "Expected anisotropy sensibility; E_{min} [GeV]; #delta_{M}^{exp}", N, min, max);

    sens(h_sens, h_p_Neff_int);

    TCanvas *c6= new TCanvas("c6","Electron fit",10,10,900,600);
    c6->SetGrid();
    c6->SetLogx();
    c6->SetLogy();
    
    h_sens->SetMaximum(1);
    h_sens->SetFillColor(kMagenta-7);
    h_sens->SetFillStyle(3001);
    h_sens->Draw();


    //Incluyamos en el .root los flujos escalados a E^3

    TF1 *sp_flux_E3 = new TF1("sp_flux_E3", Source_p_E3, 0.5, 1500.,4);
    sp_flux_E3->SetParameters(1.10, 6.80e-05, -2.58, 1.23e-3);
    sp_flux_E3->SetParNames("#varphi_{e^{-}}", "C_{a}", "#gamma_{a}", "E_{t}", "#Delta #gamma_{t}");
        
    TF1 *dp_flux_E3 = new TF1("dp_flux_E3", Diffusion_p_E3, 0.5, 1500.,3);
    dp_flux_E3->SetParameters(1.10, 6.51e-02, -4.07);
    dp_flux_E3->SetParNames("#varphi_{e^{-}}", "C_{b}", "#gamma_{b}", "E_{t}", "#Delta #gamma_{t}");


    // Hagamos print de algunas magnitudes

    print_hist(h_p_sou_int, "events B source");
    print_hist(h_p_dif_int, "events B diff");
    print_hist(h_p_Ntot_int, "Ntot int");
    print_hist(h_p_Neff_int, "Neff int");
    print_hist(h_sens, "Sensitivity");

    TFile f("Neff.root", "new");
    h_p_Ntot->Write();
    h_p_Ntot_int->Write();
    h_p_Neff->Write();
    h_p_Neff_int->Write();
    h_p_pur->Write();
    h_p_pur_int->Write();
    h_p_sou->Write();
    h_p_sou_int->Write();
    h_p_dif->Write();
    h_p_dif_int->Write();
    h_sens->Write();
    dp_flux->Write();
    sp_flux->Write();
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

void Ntot(TH1 *tot, TH1 *sou, TH1 *diff){

    Int_t N = tot->GetNbinsX();

    if ((N== sou->GetNbinsX()) && (N == diff->GetNbinsX())){
      
        for(Int_t i = 0; i < N; i++){

            Double_t xsou   = sou->GetBinContent(i);
            Double_t xdiff  = diff->GetBinContent(i);

            tot->SetBinContent(i, (xsou + xdiff));
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

  Int_t v[10] = {1, 5, 10, 16, 20, 50, 100, 200, 500, 1000};

  for(Int_t i = 0; i < 10; i++){
      
      Int_t x = hist->FindFixBin(v[i]);
      cout << "\t" << hist->GetBinContent(x) << "\t" << hist->GetBinCenter(x) << endl;
  }
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