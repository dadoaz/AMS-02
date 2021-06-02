/////////////////////////////////////////////////////////////////////////
//
// Plotting all historams from caluclus of anisotropy sensitivity
// for models, AMS and calorimetric experiments
//
// Daniel Doménech Azorín dadoa@alumni.uv.es
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

void     print_hist(TH1 *hist, TString name);

void plot_sens_exp() {

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

    // Read the histogramos from AMS analysis of positron dipolar anisotropy
    auto* file2 = new TFile("Neff.root");
    TH1D* h_bkg_ams = (TH1D*)file2->Get("Diffusion_p");
    TH1D* h_p_sou_ams = (TH1D*)file2->Get("Source_p");
    TH1D* h_bkg_int_ams = (TH1D*)file2->Get("Diffusion_p_int");
    TH1D* h_p_sou_int_ams = (TH1D*)file2->Get("Source_p_int");
    TH1D* h_pur_ams = (TH1D*)file2->Get("Purity");
    TH1D* h_pur_int_ams = (TH1D*)file2->Get("Purity integrated");
    TH1D* h_Neff_ams = (TH1D*)file2->Get("Neff");
    TH1D* h_Neff_int_ams = (TH1D*)file2->Get("Neff int");
    TH1D* h_sens_ams = (TH1D*)file2->Get("sensitivity");


    //Read the histograms from calorimeter experiments of dipolar anisotropy
    auto* file3 = new TFile("Neff_exp.root");

    TH1D* h_bkg_calet = (TH1D*)file3->Get("Bkg_ep_calet");
    TH1D* h_p_sou_calet = (TH1D*)file3->Get("Source_p_calet");
    TH1D* h_bkg_int_calet = (TH1D*)file3->Get("Bkg_ep_int_calet");
    TH1D* h_p_sou_int_calet = (TH1D*)file3->Get("Source_p_int_calet");
    TH1D* h_pur_calet = (TH1D*)file3->Get("Purity_calet");
    TH1D* h_pur_int_calet = (TH1D*)file3->Get("Purity_integrated_calet");
    TH1D* h_Neff_calet = (TH1D*)file3->Get("Neff_calet");
    TH1D* h_Neff_int_calet = (TH1D*)file3->Get("Neff_int_calet");
    TH1D* h_sens_calet = (TH1D*)file3->Get("sensitivity_calet");
    TH1D* h_e_sou_calet = (TH1D*)file3->Get("Source_e_calet");
    TH1D* h_e_dif_calet = (TH1D*)file3->Get("Diffusion_e_calet");
    TH1D* h_p_dif_calet = (TH1D*)file3->Get("Diffusion_p_calet");

    TH1D* h_bkg_dampe = (TH1D*)file3->Get("Bkg_ep_dampe");
    TH1D* h_p_sou_dampe = (TH1D*)file3->Get("Source_p_dampe");
    TH1D* h_bkg_int_dampe = (TH1D*)file3->Get("Bkg_ep_int_dampe");
    TH1D* h_p_sou_int_dampe = (TH1D*)file3->Get("Source_p_int_dampe");
    TH1D* h_pur_dampe = (TH1D*)file3->Get("Purity_dampe");
    TH1D* h_pur_int_dampe = (TH1D*)file3->Get("Purity_integrated_dampe");
    TH1D* h_Neff_dampe = (TH1D*)file3->Get("Neff_dampe");
    TH1D* h_Neff_int_dampe = (TH1D*)file3->Get("Neff_int_dampe");
    TH1D* h_sens_dampe = (TH1D*)file3->Get("sensitivity_dampe");
    TH1D* h_e_sou_dampe = (TH1D*)file3->Get("Source_e_dampe");
    TH1D* h_e_dif_dampe = (TH1D*)file3->Get("Diffusion_e_dampe");
    TH1D* h_p_dif_dampe = (TH1D*)file3->Get("Diffusion_p_dampe");

    TH1D* h_bkg_herd = (TH1D*)file3->Get("Bkg_ep_herd");
    TH1D* h_p_sou_herd = (TH1D*)file3->Get("Source_p_herd");
    TH1D* h_bkg_int_herd = (TH1D*)file3->Get("Bkg_ep_int_herd");
    TH1D* h_p_sou_int_herd = (TH1D*)file3->Get("Source_p_int_herd");
    TH1D* h_pur_herd = (TH1D*)file3->Get("Purity_herd");
    TH1D* h_pur_int_herd = (TH1D*)file3->Get("Purity_integrated_herd");
    TH1D* h_Neff_herd = (TH1D*)file3->Get("Neff_herd");
    TH1D* h_Neff_int_herd = (TH1D*)file3->Get("Neff_int_herd");
    TH1D* h_sens_herd = (TH1D*)file3->Get("sensitivity_herd");
    TH1D* h_e_sou_herd = (TH1D*)file3->Get("Source_e_herd");
    TH1D* h_e_dif_herd = (TH1D*)file3->Get("Diffusion_e_herd");
    TH1D* h_p_dif_herd = (TH1D*)file3->Get("Diffusion_p_herd");

    

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
                                        /* DRAW */
    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    TCanvas *c1= new TCanvas("c1","c1",10,10,900,600);
	c1->SetGrid();
	c1->SetLogx();
    c1->SetLogy();

    //Create a legend
    TLegend* leg = new TLegend(0.6, 0.75, 0.9, 0.9);
    leg->AddEntry(h_p_sou_calet, "e^{+} Source");
    leg->AddEntry(h_p_dif_calet, "e^{+} Diffusion");
    leg->AddEntry(h_e_dif_calet, "e^{-} power law A");
    leg->AddEntry(h_e_sou_calet, "e^{-} power law B");
    leg->AddEntry(h_bkg_calet,   "Total isotropic background");

    h_e_dif_calet->SetMinimum(1);
	//h_e_dif_calet->SetFillStyle(3144);
    h_e_dif_calet->SetLineColor(kCyan);
    h_e_dif_calet->SetLineWidth(2);
    
    //h_e_sou_calet->SetFillStyle(3144);
    h_e_sou_calet->SetLineColor(kGreen);
    h_e_sou_calet->SetLineWidth(2);

    //h_p_dif_calet->SetFillStyle(3144);
    h_p_dif_calet->SetLineColor(kOrange+10);
    h_p_dif_calet->SetLineWidth(2);
    h_p_dif_calet->SetLineWidth(2);

    h_p_sou_calet->SetFillStyle(3001);
    h_p_sou_calet->SetFillColor(kOrange-3);
    

    h_bkg_calet->SetMinimum(1);
    h_bkg_calet->SetLineWidth(4);
    h_bkg_calet->SetLineColor(kBlack);
    
    h_bkg_calet->SetTitle("Expected e^{+} e^{-} events @ CALET");
    h_bkg_calet->Draw("][");
    h_e_dif_calet->Draw("][ same");
    h_e_sou_calet->Draw("][ same");
    h_p_dif_calet->Draw("][ same");
    h_p_sou_calet->Draw("][ same");
    leg->Draw("same");

    TCanvas *c10= new TCanvas("c10","c10",10,10,900,600);
	c10->SetGrid();
	c10->SetLogx();
    c10->SetLogy();

    //Create a legend
    TLegend* leg10 = new TLegend(0.6, 0.75, 0.9, 0.9);
    leg10->AddEntry(h_p_sou_dampe, "e^{+} Source");
    leg10->AddEntry(h_p_dif_dampe, "e^{+} Diffusion");
    leg10->AddEntry(h_e_dif_dampe, "e^{-} power law A");
    leg10->AddEntry(h_e_sou_dampe, "e^{-} power law B");
    leg10->AddEntry(h_bkg_dampe,   "Total isotropic background");

    h_e_dif_dampe->SetMinimum(1);
	//h_e_dif_dampe->SetFillStyle(3144);
    h_e_dif_dampe->SetLineColor(kCyan);
    h_e_dif_dampe->SetLineWidth(2);
    
    //h_e_sou_dampe->SetFillStyle(3144);
    h_e_sou_dampe->SetLineColor(kGreen);
    h_e_sou_dampe->SetLineWidth(2);

    //h_p_dif_dampe->SetFillStyle(3144);
    h_p_dif_dampe->SetLineColor(kOrange+10);
    h_p_dif_dampe->SetLineWidth(2);
    h_p_dif_dampe->SetLineWidth(2);

    h_p_sou_dampe->SetFillStyle(3001);
    h_p_sou_dampe->SetFillColor(kOrange-3);
    

    h_bkg_dampe->SetMinimum(1);
    h_bkg_dampe->SetLineWidth(4);
    h_bkg_dampe->SetLineColor(kBlack);
    
    h_bkg_dampe->SetTitle("Expected e^{+} e^{-} events @ DAMPE");
    h_bkg_dampe->Draw("][");
    h_e_dif_dampe->Draw("][ same");
    h_e_sou_dampe->Draw("][ same");
    h_p_dif_dampe->Draw("][ same");
    h_p_sou_dampe->Draw("][ same");
    leg10->Draw("same");    

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    TCanvas *c2= new TCanvas("c2","c2",10,10,900,600);
    c2->SetGrid();
    c2->SetLogx();
    c2->SetLogy();
    
    h_bkg_int_calet->SetMinimum(1);
	h_bkg_int_calet->SetFillStyle(3144);
    h_bkg_int_calet->SetFillColor(kGray+3);

    h_p_sou_int_calet->SetFillStyle(3001);
    h_p_sou_int_calet->SetFillColor(kOrange-3);

    h_bkg_int_calet->Draw();
    h_p_sou_int_calet->Draw("same");
    //leg->Draw("same");

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    TCanvas *c3= new TCanvas("c3","c3",10,10,900,600);
    c3->SetGrid();
    c3->SetLogx();
    //c3->SetLogy();

    TLegend* leg30 = new TLegend(0.1, 0.8, 0.4, 0.9);
    leg30->AddEntry(h_pur_calet, "CALET - DAMPE - HERD");
    leg30->AddEntry(h_pur_ams, "AMS");
    

    h_pur_calet->SetLineWidth(2);
    h_pur_calet->SetFillColor(kCyan-7);
    h_pur_calet->SetFillStyle(3001);
    
    h_pur_ams->SetFillColorAlpha(kRed, 0.4);
    h_pur_ams->SetLineColor(kRed);
    h_pur_ams->SetLineWidth(2);
    h_pur_ams->SetTitle("Differental purity");
    h_pur_ams->Draw("][");
    h_pur_calet->SetMaximum(1);
    h_pur_calet->Draw("same ][");//"pfc plc"
    h_pur_dampe->Draw("same ][");
    h_pur_herd->Draw("same  ][");
    leg30->Draw("same");
    
    
    TCanvas *c4= new TCanvas("c4","c4",10,10,900,600);
    c4->SetGrid();
    c4->SetLogx();
    //c4->SetLogy();

    h_pur_int_calet->SetFillColor(kMagenta-7);
    h_pur_int_calet->SetFillStyle(3001);
    h_pur_int_calet->Draw();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    TCanvas *c5= new TCanvas("c5","c5",10,10,900,600);
    c5->SetGrid();
    c5->SetLogx();
    c5->SetLogy();

    h_Neff_calet->SetMinimum(1);
    h_Neff_calet->SetFillColor(kCyan-7);
    h_Neff_calet->SetFillStyle(3001);
    h_Neff_calet->Draw();
    
    TCanvas *c6= new TCanvas("c6","c6",10,10,900,600);
    c6->SetGrid();
    c6->SetLogx();
    c6->SetLogy();
    
    h_Neff_int_calet->SetMinimum(1);
    h_Neff_int_calet->SetFillColor(kMagenta-7);
    h_Neff_int_calet->SetFillStyle(3001);
    h_Neff_int_calet->Draw();
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////

    TCanvas *c7= new TCanvas("c7","c7",10,10,900,700);
    c7->SetGrid();
    c7->SetLogx();
    c7->SetLogy();
   
    TLegend* leg7 = new TLegend(0.1, 0.78, 0.4, 0.9);
    leg7->AddEntry(h_sens_calet, "CALET - 6 years");
    leg7->AddEntry(h_sens_dampe, "DAMPE - 6 years");
    leg7->AddEntry(h_sens_herd, "HERD - 10 years");
    
    h_sens_calet->GetXaxis()->SetRangeUser(10, 1499);
    h_sens_calet->SetMinimum(1e-3);
    h_sens_calet->SetMaximum(1);
    
    h_sens_calet->SetMaximum(1);
    h_sens_calet->SetMinimum(1e-4);
    h_sens_calet->SetFillColor(kMagenta-7);
    h_sens_calet->SetFillStyle(3001);
    h_sens_calet->Draw();
    h_sens_dampe->SetFillColor(kCyan-7);
    h_sens_dampe->SetFillStyle(3001);
    h_sens_dampe->Draw("same");
    h_sens_herd->SetFillColor(kYellow-7);
    h_sens_herd->SetFillStyle(3001);
    h_sens_herd->Draw("same");
    leg7->Draw("same");

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    TCanvas* c22 = new TCanvas("c22", "Assym sens models", 50, 50, 900, 600);

    c22->SetLogx();
    c22->SetLogy();
    c22->SetGrid();

    h_sens_calet->SetMaximum(5e-1);
    h_sens_calet->SetMinimum(1e-3);
    h_sens_calet->GetXaxis()->SetRangeUser(10,1499.8);
    
    float alpha = 0.55;
    
    h_sens_calet->SetTitle("Expected dipole anisotropy sensitivity");
    h_sens_calet->SetLineColor(kMagenta);
    h_sens_calet->SetLineWidth(2);
    //h_sens_calet->SetFillColor(kMagenta);
    h_sens_calet->SetFillColorAlpha(kMagenta, alpha);
    h_sens_calet->SetFillStyle(1001);

    h_sens_dampe->SetLineColor(kCyan);
    h_sens_dampe->SetLineWidth(2);
    h_sens_dampe->SetFillColorAlpha(kCyan, alpha);
    //h_sens_dampe->SetFillColor(kCyan);
    h_sens_dampe->SetFillStyle(1001);

    h_sens_herd->SetLineColor(kYellow);
    h_sens_herd->SetLineWidth(2);
    h_sens_herd->SetFillColorAlpha(kYellow, alpha);
    //h_sens_herd->SetFillColor(kYellow);
    h_sens_herd->SetFillStyle(1001);
    
    
    h_sens_ams->SetLineColor(kGreen);
    h_sens_ams->SetLineWidth(2);
    h_sens_ams->SetFillColorAlpha(kGreen, alpha);
    //h_sens_ams->SetFillColor(kGreen);
    h_sens_ams->SetFillStyle(1001);

    a_Hoop_sum_corr->SetLineColor(kViolet+4);
    a_Hoop_sum_corr->SetLineWidth(3);
    a_Man_sum_corr->SetLineColor(kRed);
    a_Man_sum_corr->SetLineWidth(3);

    //h_sens_calet->SetFillStyle(0);
    //h_sens_dampe->SetFillStyle(0);
    //h_sens_herd->SetFillStyle(0);
    //h_sens_ams->SetFillStyle(0);
    
    h_sens_herd->GetXaxis()->SetRangeUser(10, 1499);
    h_sens_herd->SetMinimum(1e-3);
    h_sens_herd->SetMaximum(1);
    
    h_sens_calet->Draw("][ same");
    h_sens_dampe->Draw("][ same");
    h_sens_ams->Draw("][ same");
    h_sens_herd->Draw("][ same");
    
    
    a_Man_sum_corr->Draw("][ same");
    a_Hoop_sum_corr->Draw("][ same");

    TLegend* leg22 = new TLegend(0.1, 0.65, 0.4, 0.9);
    
    
    leg22->AddEntry(h_sens_calet, "CALET - 6 years");
    leg22->AddEntry(h_sens_dampe, "DAMPE - 6 years");
    leg22->AddEntry(h_sens_ams, "AMS - 17 years");
    leg22->AddEntry(h_sens_herd, "HERD - 10 years");
    leg22->AddEntry(a_Hoop_sum_corr, "Hooper & Linden");
    leg22->AddEntry(a_Man_sum_corr, "Manconi");

    leg22->Draw("same");

    c22->SaveAs("calo/models.pdf");



    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    TH1D *h_tot_calet = new TH1D("Ntot_calet", "", 15000, 1, 1500);
    h_tot_calet->Add(h_p_sou_int_calet, h_bkg_int_calet, 1., 1.);

    TH1D *h_tot_dampe = new TH1D("Ntot_dampe", "", 15000, 1, 1500);
    h_tot_dampe->Add(h_p_sou_int_dampe, h_bkg_int_dampe, 1., 1.);
    
    TH1D *h_tot_herd = new TH1D("Ntot_herd", "", 15000, 1, 1500);
    h_tot_herd->Add(h_p_sou_int_herd, h_bkg_int_herd, 1., 1.);
    
    TH1D *h_tot_ams = new TH1D("Ntot_ams", "", 15000, 1, 1500);
    h_tot_ams->Add(h_p_sou_int_ams, h_bkg_int_ams, 1., 1.);

    print_hist(h_sens_calet, "Calet N sou");
    print_hist(h_sens_dampe, "Dampe N sou");
    print_hist(h_sens_herd, "Herd N sou");



    TString d = ""; 
    //c1->SaveAs("positron_expected_events_Benergydep" + d + ".pdf");    //Benergydep
    //c2->SaveAs("positron_expected_survival_Benergydep" + d + ".pdf");  
    //c3->SaveAs("positron_purity_B" + d + ".pdf");    //Aconst
    //c4->SaveAs("positron_purity_integrated_B" + d + ".pdf");
    //c5->SaveAs("positron_Neff_B_Calet" + d + ".pdf");    //Aconst
    //c6->SaveAs("positron_Neff_int_B_Calet" + d + ".pdf");
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
