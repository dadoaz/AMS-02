/////////////////////////////////////////////////////////////////////////
//
// Estimación de la anisotropía de positrones en diferentes modelos de púlsares
// 
//
// Daniel Doménech Azorín dadoa@alumni.uv,es
//
// 4.15.2021
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
#include <vector>

Double_t Positron(Double_t *x, Double_t *par);
Double_t Electron(Double_t *x, Double_t *par);
Double_t Diffusion_p(Double_t *x, Double_t *par);
Double_t Diffusion_e(Double_t *x, Double_t *par);
Double_t Source_p(Double_t *x, Double_t *par);
Double_t Source_e(Double_t *x, Double_t *par);
Double_t p_frac(Double_t x);
Double_t pos_frac(Double_t *x, Double_t *par);
void     print_hist(TH1 *hist, TString name);
void     Hoop_correction(TH1 *result, TH1 *hist, TF1 *p_frac, TF1 *assym);
void     a_Hoop_correction (TH1 *result, TH1 *hist, TH1 *pur);
void     a_Man_correction (TH1 *result, TF1 *func, TH1 *pur);


struct SumTF1 { 

   SumTF1(const std::vector<TF1 *> & flist) : fFuncList(flist) {}
   
   double operator() (const double * x, const double *p) {
      double result = 0;
      for (unsigned int i = 0; i < fFuncList.size(); ++i) 
         result += fFuncList[i]->EvalPar(x,p); 
      return result; 
   } 
   
   std::vector<TF1*> fFuncList; 
      
};





Int_t plot_sens(){

    gStyle->SetOptStat(0);
	gStyle->SetOptFit(1111);
    gStyle->SetPalette(kSolar);

    ////Read the TH1 of processed histograms from pulsar models
    auto* file1 = new TFile("Assym_sens.root");
    
    TH1D* a_Hoop_sum_corr = (TH1D*)file1->Get("assym_Hoop_sum");
    TH1D* a_Hoop_Gem_corr = (TH1D*)file1->Get("assym_Hoop_Gem");
    TH1D* a_Hoop_Mon_corr = (TH1D*)file1->Get("assym_Hoop_Mon");
    TH1D* a_Man_sum_corr  = (TH1D*)file1->Get("assym_Man_sum");
    TH1D* a_Man_Gem_corr  = (TH1D*)file1->Get("assym_Man_Gem");
    TH1D* a_Man_Mon_corr  = (TH1D*)file1->Get("assym_Man_Mon");
    TH1D* h_sens_ams      = (TH1D*)file1->Get("sensitivity");

    //Read TF1 of original pulsar model functions data. Both separate between Geminga & Monogem 
    auto* file2 = new TFile("HooperLinden.root");
    auto* file3 = new TFile("Manconi.root");
    
    TF1* Hoop_Gem = (TF1*)file2->Get("fGeminga");
    TF1* Hoop_Mon = (TF1*)file2->Get("fMonogem");
    TF1* Man_Gem = (TF1*)file3->Get("fGeminga");
    TF1* Man_Mon = (TF1*)file3->Get("fMonogem");

    //Create the sum function of both contributions of Geminga & Monogem
    //We are not using it 

    std::vector<TF1 *> v;
    v.push_back((TF1*)file2->Get("fGeminga"));
    v.push_back((TF1*)file2->Get("fMonogem"));

    std::vector<TF1 *> w;
    w.push_back((TF1*)file3->Get("fGeminga"));
    w.push_back((TF1*)file3->Get("fMonogem"));


    TF1* Hoop_sum = new TF1("Hoop_sum", SumTF1(v), 0.5, 1500, 0);
    TF1* Man_sum =  new TF1("Man_sum",  SumTF1(w), 0.5, 1500, 0);

    //Plot the results

    TCanvas* c1 = new TCanvas("c1", "Assym sens models", 50, 50, 900, 600);

    c1->SetLogx();
    c1->SetLogy();
    c1->SetGrid();

    
    h_sens->SetMinimum(1e-4);
    h_sens->GetXaxis()->SetRangeUser(1,1499.8);
    h_sens->Draw("][ same");
    h_sens->SetTitle("Expected e^{+} anisotropy sensitivity");

    
    h_sens->SetLineColor(kAzure+1);
    h_sens->SetLineWidth(2);
    h_sens->SetFillColor(kAzure+1);
    a_Hoop_sum_corr->SetLineColor(kMagenta);
    a_Hoop_sum_corr->SetLineWidth(3);
    a_Man_sum_corr->SetLineColor(kViolet+4);
    a_Man_sum_corr->SetLineWidth(3);
    
    a_Man_sum_corr->Draw("][ same");
    a_Hoop_sum_corr->Draw("][ same");

    TLegend* leg = new TLegend(0.1, 0.75, 0.4, 0.9);
    
    
    leg->AddEntry(a_Hoop_sum_corr, "Hooper & Linden");
    leg->AddEntry(a_Man_sum_corr, "Manconi");
    leg->AddEntry(h_sens, "AMS (2028)");

    leg->Draw("same");

    c1->SaveAs("models.pdf");

    print_hist(h_sens, "AMS");
    print_hist(a_Hoop_sum_corr, "Hoop");
    print_hist(a_Man_sum_corr, "Man");

    TCanvas* c2 = new TCanvas("c2", "Geminga", 50, 50, 900, 600);

    c2->SetLogx();
    c2->SetLogy();
    c2->SetGrid();

    a_Man_sum_corr->SetLineColor(kBlue);
    a_Man_Mon_corr->SetLineColor(kAzure+1);
    Man_Mon->SetLineColor(kAzure+1);
    Man_Mon->SetLineStyle(3);
    a_Man_Gem_corr->SetLineColor(kAzure+10);
    Man_Gem->SetLineColor(kAzure+10);
    Man_Gem->SetLineStyle(3);
    
    a_Man_sum_corr->SetLineWidth(4);
    a_Man_Mon_corr->SetLineWidth(3);
    a_Man_Gem_corr->SetLineWidth(3);

    a_Man_Gem_corr->SetMaximum(1e-1);
    a_Man_Gem_corr->SetMinimum(1e-4);
    a_Man_Gem_corr->GetXaxis()->SetRangeUser(1,1499.8);
    a_Man_Gem_corr->SetTitle("Manconi");

    
    a_Man_Gem_corr->Draw("][ same");
    a_Man_Mon_corr->Draw("][ same");
    Man_Gem->Draw("C same");
    Man_Mon->Draw("C same");
    a_Man_sum_corr->Draw("][ same");

    TLegend* leg2 = new TLegend(0.1, 0.7, 0.4, 0.9);

    leg2->AddEntry(a_Man_sum_corr, "Total (source)");
    leg2->AddEntry(a_Man_Mon_corr, "Monogem (source)");
    leg2->AddEntry(Man_Mon, "Monogem (e^{+})");
    leg2->AddEntry(a_Man_Gem_corr, "Geminga (source)");
    leg2->AddEntry(Man_Gem, "Geminga (e^{+})");
    

    leg2->Draw("same");

    c2->SaveAs("Hooper.pdf");

    TCanvas* c3 = new TCanvas("c3", "Geminga", 50, 50, 900, 600);

    c3->SetLogx();
    c3->SetLogy();
    c3->SetGrid();

    a_Hoop_sum_corr->SetLineColor(kRed);
    a_Hoop_Mon_corr->SetLineColor(kOrange+1);
    Hoop_Mon->SetLineColor(kOrange+1);
    Hoop_Mon->SetLineStyle(3);
    a_Hoop_Gem_corr->SetLineColor(kOrange);
    Hoop_Gem->SetLineColor(kOrange);
    Hoop_Gem->SetLineStyle(3);
    
    a_Hoop_sum_corr->SetLineWidth(4);
    a_Hoop_Mon_corr->SetLineWidth(3);
    a_Hoop_Gem_corr->SetLineWidth(3);

    a_Hoop_Gem_corr->SetMaximum(1e-1);
    a_Hoop_Gem_corr->SetMinimum(1e-4);
    a_Hoop_Gem_corr->GetXaxis()->SetRangeUser(1,1499.8);
    a_Hoop_Gem_corr->SetTitle("Hooper & Linden");

    
    a_Hoop_Gem_corr->Draw("][ same");
    a_Hoop_Mon_corr->Draw("][ same");
    Hoop_Gem->Draw("C same");
    Hoop_Mon->Draw("C same");
    a_Hoop_sum_corr->Draw("][ same");

    TLegend* leg3 = new TLegend(0.1, 0.7, 0.4, 0.9);

    leg3->AddEntry(a_Hoop_sum_corr, "Total (source)");
    leg3->AddEntry(a_Hoop_Mon_corr, "Monogem (source)");
    leg3->AddEntry(Hoop_Mon, "Monogem (e^{+})");
    leg3->AddEntry(a_Hoop_Gem_corr, "Geminga (source)");
    leg3->AddEntry(Hoop_Gem, "Geminga (e^{+})");
    

    leg3->Draw("same");

    c3->SaveAs("Manconi.pdf");
    
    
    return 0;
    
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

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


Double_t p_frac(Double_t x){

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

    Double_t p = sp_flux->Eval(x) + dp_flux->Eval(x);
    Double_t e = se_flux->Eval(x) + de_flux->Eval(x);
    
    return p / (p + e);
}

Double_t pos_frac(Double_t *x, Double_t *par){
    
    Double_t E = x[0];

    return p_frac(E);
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

void  Hoop_correction(TH1 *result, TH1 *hist, TF1 *p_frac, TF1 *assym){

    Double_t x, y;
    Int_t N = hist->GetNbinsX();

    TH1D *test = (TH1D*)hist->Clone("test");

    for(Int_t i = 0; i < N; i++){

        x = hist->GetBinCenter(i);

        y = hist->GetBinContent(i) * assym->Eval(x) / (2 * p_frac->Eval(x));

        test->SetBinContent(i, y);
    }

    if( N == result->GetNbinsX() ){
        
        for(Int_t i = 0; i < N; i++){
            
            y = test->Integral(i, N);
            y = y/hist->Integral(i, N);

            result-> SetBinContent(i,y);

        } 
    }
}

void a_Hoop_correction (TH1 *result, TH1 *hist, TH1 *pur){
    
    Double_t x, y;
    Int_t N = hist->GetNbinsX();

    if( N == result->GetNbinsX() ){
        
        for(Int_t i = 0; i < N; i++){
          
            y = hist->GetBinContent(i) / pur->GetBinContent(i);

            result-> SetBinContent(i,y);

        } 
    }

}
void a_Man_correction (TH1 *result, TF1 *func, TH1 *pur){

    Double_t x, y;
    Int_t N = result->GetNbinsX();

    if( N == pur->GetNbinsX() ){

        for(Int_t i = 0; i < N; i++){

                x = result->GetBinCenter(i);

                y = func->Eval(x) / pur->GetBinContent(i);

                result-> SetBinContent(i,y);

        } 
    }
}
