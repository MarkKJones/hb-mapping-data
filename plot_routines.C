#include <iostream>
#include <fstream>
#include <TString.h>
#include <TLegend.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TMath.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TObjArray.h>
#include <iterator>
#include <vector>
#include <algorithm>
using namespace std;

void Get_two_col_data_from_file(TString filename,vector<double> &,vector<double> & );
void Plot_hb_bi(vector<double> & ,vector<double> & );
void Set_probe(vector<double> & ,vector<double> & );
void Plot_zscan(TString current, TString nhole,vector<double> & ,vector<double> & );
void Plot_zscan_pos(TString current,vector<double> & ,vector<double> & ,vector<double> & ,vector<double> & );
Double_t field_corr(Double_t field );
void run_plot_zscan(TString c, TString n);
// Global parameters
vector<double> nominal_field_corfile_data;
vector<double> corr_field_corfile_data;
vector<double> HB_central_current_data;
vector<double> HB_central_field_data;
vector<double> HB_central_zscan_data_dis;
vector<double> HB_central_zscan_data;
vector<double> HB_central_zscan_tosca_dis;
vector<double> HB_central_zscan_tosca;
TF1 *neg_field_corr;
TF1 *pos_field_corr;
TF1 *cur_field;
Bool_t probe_not_set = kTRUE;
//
void run_plot_zscan(TString current, TString nhole) {
  if (probe_not_set) {
      Get_two_col_data_from_file("Lakeshore-probe.dat",nominal_field_corfile_data,corr_field_corfile_data);
      Set_probe(nominal_field_corfile_data,corr_field_corfile_data);
  }
  TString dfile;
  dfile="Central_zscan_files/central_zscan_"+current+"a_hole"+nhole+".dat";
  Get_two_col_data_from_file(dfile,HB_central_zscan_data,HB_central_zscan_tosca);
  if (HB_central_zscan_data.size()==57 && HB_central_zscan_tosca.size()==57) {
  Plot_zscan(current,nhole,HB_central_zscan_data,HB_central_zscan_tosca);
  } else {
    cout << " vector size below 57" << HB_central_zscan_data.size() << " " << HB_central_zscan_tosca.size() << endl;
  }
  HB_central_zscan_data.clear();
  HB_central_zscan_tosca.clear();
}
//
void plot_zscan_hole() {
   Int_t colh;
   vector<double> tosca_int;
   vector<double>  tosca_max;
   vector<double>  tosca_effl;
   vector<double>  data_int;
   vector<double>  data_max;
   vector<double>  data_effl;
   vector<double>  rat_int;
   vector<double>  rat_max;
   vector<double>  rat_effl;
   vector<double>  cur;
  Float_t col1,col2,col3,col4,col5;
   TString dnegfile="Central_zscan_files/central_zscan_neg_results.dat";
   TTree *T1 = new TTree("ntuple","data");
   Long64_t nlines = T1->ReadFile(dnegfile,"col1:colh/I:col2/F:col3:col4:col5",' ');
    T1->SetBranchAddress("col1",&col1);
    T1->SetBranchAddress("colh",&colh);
   T1->SetBranchAddress("col2",&col2);
   T1->SetBranchAddress("col3",&col3);
   T1->SetBranchAddress("col4",&col4);
   T1->SetBranchAddress("col5",&col5);
   TGraph *gr_rat_int[5];
   TGraph *gr_data_int[5];
      TGraph *gr_tosca_int[5];
   TGraph *gr_rat_max[5];
   TGraph *gr_data_max[5];
      TGraph *gr_tosca_max[5];
   TGraph *gr_rat_effl[5];
   TGraph *gr_data_effl[5];
      TGraph *gr_tosca_effl[5];
   for (Int_t nh = 1; nh <6;nh++) {
   for (Int_t i = 0; i <nlines;i++) {
     T1->GetEntry(i);
     if (colh ==nh) {
     cout << col1 << " "  << colh << " "  << col2/col1  << " "  << col3 << " "  << col4 << " "  << col5 << " " << endl; 
     cur.push_back(col1);
     data_max.push_back(TMath::Abs(col2/col1));
     data_int.push_back(TMath::Abs(col3/col1));
     data_effl.push_back(TMath::Abs(col3/col2)*100.);
     tosca_max.push_back(TMath::Abs(col4/col1));
     tosca_int.push_back(TMath::Abs(col5/col1));
     tosca_effl.push_back(TMath::Abs(col5/col4)*100.);
     rat_int.push_back(TMath::Abs(col3/col5));
     rat_max.push_back(TMath::Abs(col2/col4));
     rat_effl.push_back(TMath::Abs(col3/col2/(col5/col4)));
     }
     cout << cur.size() << endl;
     gr_rat_int[nh-1] = new TGraph(cur.size(),&(cur[0]),&(rat_int[0]));  
     gr_data_int[nh-1] = new TGraph(cur.size(),&(cur[0]),&(data_int[0])); 
     gr_tosca_int[nh-1] = new TGraph(cur.size(),&(cur[0]),&(tosca_int[0])); 
     gr_rat_effl[nh-1] = new TGraph(cur.size(),&(cur[0]),&(rat_effl[0]));  
     gr_data_effl[nh-1] = new TGraph(cur.size(),&(cur[0]),&(data_effl[0])); 
     gr_tosca_effl[nh-1] = new TGraph(cur.size(),&(cur[0]),&(tosca_effl[0])); 
     gr_rat_max[nh-1] = new TGraph(cur.size(),&(cur[0]),&(rat_max[0]));  
     gr_data_max[nh-1] = new TGraph(cur.size(),&(cur[0]),&(data_max[0])); 
     gr_tosca_max[nh-1] = new TGraph(cur.size(),&(cur[0]),&(tosca_max[0])); 
     cur.clear();
     data_max.clear();
     data_int.clear();
     data_effl.clear();
     tosca_max.clear();
     tosca_int.clear();
     tosca_effl.clear();
     rat_int.clear();
     rat_max.clear();
     rat_effl.clear();
   }
   }
      T1->Delete();
 //
 gROOT->Reset();
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.5,"Y");
 gStyle->SetTitleOffset(.8,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.1);
    TMultiGraph *mg_int = new TMultiGraph();
 TLegend *leg_int = new TLegend(.3,.2,.66,.4);
 TCanvas *can = new TCanvas("can","Int Bdl ",800,800);
   can->Divide(1,2);
   for (Int_t nh = 1; nh <6;nh++) {
       leg_int->AddEntry(gr_data_int[nh-1],Form(" hole %d",nh),"p");
       mg_int->Add(gr_data_int[nh-1]);  
       gr_data_int[nh-1]->SetMarkerStyle(22);
       gr_data_int[nh-1]->SetMarkerSize(1.5);
       gr_data_int[nh-1]->SetMarkerColor(nh);
   }
 TString htitle="HB Integral Bdl/Current";
 mg_int->Draw("AP");
 mg_int->SetTitle(htitle);
 mg_int->GetXaxis()->SetTitle("Current");
 mg_int->GetYaxis()->SetTitle("Int Bdl/I (Tm/I)");
 leg_int->Draw();
   //
}
//
void plot_tosca_comparison() {
 gROOT->Reset();
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.5,"Y");
 gStyle->SetTitleOffset(.8,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.1);
   vector <double> HB_AB_cur;
   vector <double> HB_AB_tosca_int;
   vector<double>  HB_AB_tosca_max;
   vector<double>  HB_AB_tosca_effl;
   vector <double> tosca_rat_cur;
   vector <double> tosca_rat_int;
   vector<double>  tosca_rat_max;
   vector<double>  tosca_rat_effl;
   TTree *T2 = new TTree("ntuple","data");
  Float_t col1,col2,col3;
Long64_t nlines = T2->ReadFile("Central_zscan_files/tosca_HB-AB-cases.dat","col1:col2:col3",' ');
    T2->SetBranchAddress("col1",&col1);
   T2->SetBranchAddress("col2",&col2);
   T2->SetBranchAddress("col3",&col3);
   for (Int_t i = 0; i <nlines;i++) {
     T2->GetEntry(i);
     HB_AB_cur.push_back(col1);
     HB_AB_tosca_max.push_back(TMath::Abs(col2/col1));
     HB_AB_tosca_int.push_back(TMath::Abs(col3/col1));
     HB_AB_tosca_effl.push_back(TMath::Abs(col3/col2)*100.);
     cout << col1 << " " << col2 << " " << col2/col1 << endl;
   }
    T2->Delete();
//
   vector <double> HB_V10_cur;
   vector <double> HB_V10_tosca_int;
   vector<double>  HB_V10_tosca_max;
   vector<double>  HB_V10_tosca_effl;
   TTree *T3 = new TTree("ntuple","data");
   nlines = T3->ReadFile("Central_zscan_files/tosca_HB-V10-cases.dat","col1:col2:col3",' ');
    T3->SetBranchAddress("col1",&col1);
   T3->SetBranchAddress("col2",&col2);
   T3->SetBranchAddress("col3",&col3);
   for (Int_t i = 0; i <nlines;i++) {
     T3->GetEntry(i);
     HB_V10_cur.push_back(col1);
     HB_V10_tosca_max.push_back(TMath::Abs(col2/col1));
     HB_V10_tosca_int.push_back(TMath::Abs(col3/col1));
     HB_V10_tosca_effl.push_back(TMath::Abs(col3/col2)*100.);
     cout << col1 << " " << col2 << " " << col2/col1 << endl;
   }
    T3->Delete();
//
    for (UInt_t i = 0; i <HB_AB_cur.size();i++) {
    for (UInt_t j = 0; j <HB_V10_cur.size();j++) {
      if (HB_AB_cur[i] == HB_V10_cur[j]) {
     tosca_rat_cur.push_back(HB_V10_cur[j]);
     tosca_rat_max.push_back(HB_V10_tosca_max[j]/HB_AB_tosca_max[i]);
     tosca_rat_int.push_back(HB_V10_tosca_int[j]/HB_AB_tosca_int[i]);
     tosca_rat_effl.push_back(HB_V10_tosca_effl[j]/HB_AB_tosca_effl[i]);
     cout << tosca_rat_cur[tosca_rat_max.size()-1] << " " << tosca_rat_max[tosca_rat_max.size()-1] << " " <<  tosca_rat_int[tosca_rat_max.size()-1]<<  " " << tosca_rat_effl[tosca_rat_max.size()-1] << endl;
      }
    } 
      }      
//
 TCanvas *can = new TCanvas("can","Int Bdl ",800,800);
   can->Divide(1,3);
    TGraph *gr_HB_AB_tosca_int = new TGraph(HB_AB_cur.size(),&(HB_AB_cur[0]),&(HB_AB_tosca_int[0]));
    TGraph *gr_HB_V10_tosca_int = new TGraph(HB_V10_cur.size(),&(HB_V10_cur[0]),&(HB_V10_tosca_int[0]));
    TF1 *fit_V10_tosca_int;
    gr_HB_V10_tosca_int->Fit("pol2","","",0,4000);
    TGraph *gr_HB_AB_tosca_max = new TGraph(HB_AB_cur.size(),&(HB_AB_cur[0]),&(HB_AB_tosca_max[0]));
    TGraph *gr_HB_V10_tosca_max = new TGraph(HB_V10_cur.size(),&(HB_V10_cur[0]),&(HB_V10_tosca_max[0]));
    TGraph *gr_HB_AB_tosca_effl = new TGraph(HB_AB_cur.size(),&(HB_AB_cur[0]),&(HB_AB_tosca_effl[0]));
    TGraph *gr_HB_V10_tosca_effl = new TGraph(HB_V10_cur.size(),&(HB_V10_cur[0]),&(HB_V10_tosca_effl[0]));
    TMultiGraph *mg_int = new TMultiGraph();
    TMultiGraph *mg_max = new TMultiGraph();
    TMultiGraph *mg_effl = new TMultiGraph();
 TLegend *leg_int = new TLegend(.3,.2,.66,.4);
 TLegend *leg_max = new TLegend(.3,.2,.66,.4);
 leg_int->AddEntry(gr_HB_AB_tosca_int,"Tosca HB-AB","p");
 leg_int->AddEntry(gr_HB_V10_tosca_int,"Tosca HB-V10","p");
 mg_int->Add(gr_HB_AB_tosca_int);  
 mg_int->Add(gr_HB_V10_tosca_int);  
 TString htitle="HB Integral Bdl/Current";
 gr_HB_AB_tosca_int->SetMarkerStyle(22);
 gr_HB_AB_tosca_int->SetMarkerSize(1.5);
 gr_HB_AB_tosca_int->SetMarkerColor(1);
 gr_HB_V10_tosca_int->SetMarkerStyle(21);
 gr_HB_V10_tosca_int->SetMarkerSize(1.5);
 gr_HB_V10_tosca_int->SetMarkerColor(2);
 can->cd(1);
 mg_int->Draw("AP");
 mg_int->SetTitle(htitle);
 mg_int->GetXaxis()->SetTitle("Current");
 mg_int->GetYaxis()->SetTitle("Int Bdl/I (Tm/I)");
 leg_int->Draw();
 can->cd(2);
 mg_max->Add(gr_HB_AB_tosca_max);  
 mg_max->Add(gr_HB_V10_tosca_max);  
 leg_max->AddEntry(gr_HB_AB_tosca_max,"Tosca HB-AB","p");
 leg_max->AddEntry(gr_HB_V10_tosca_max,"Tosca HB-V10","p");
 htitle="HB Max B";
 gr_HB_AB_tosca_max->SetMarkerStyle(22);
 gr_HB_AB_tosca_max->SetMarkerSize(1.5);
 gr_HB_AB_tosca_max->SetMarkerColor(1);
 gr_HB_V10_tosca_max->SetMarkerStyle(21);
 gr_HB_V10_tosca_max->SetMarkerSize(1.5);
 gr_HB_V10_tosca_max->SetMarkerColor(2);
    gr_HB_V10_tosca_max->Fit("pol3","","",500,4000);
 mg_max->Draw("AP");
 mg_max->SetTitle(htitle);
 mg_max->GetXaxis()->SetTitle("Current");
 mg_max->GetYaxis()->SetTitle("Max B (T)");
 leg_max->Draw();
 can->cd(3);
 mg_effl->Add(gr_HB_AB_tosca_effl);  
 mg_effl->Add(gr_HB_V10_tosca_effl);  
 htitle="HB Eff length";
 gr_HB_AB_tosca_effl->SetMarkerStyle(22);
 gr_HB_AB_tosca_effl->SetMarkerSize(1.5);
 gr_HB_AB_tosca_effl->SetMarkerColor(1);
 gr_HB_V10_tosca_effl->SetMarkerStyle(21);
 gr_HB_V10_tosca_effl->SetMarkerSize(1.5);
 gr_HB_V10_tosca_effl->SetMarkerColor(2);
    gr_HB_V10_tosca_effl->Fit("pol3","","",500,4000);
 mg_effl->Draw("AP");
 mg_effl->SetTitle(htitle);
 mg_effl->GetXaxis()->SetTitle("Current");
 mg_effl->GetYaxis()->SetTitle("Eff Length (cm)");
 //
   vector <double> HB_cur;
   vector <double> HB_mom;
   vector <double> HB_ratio;
   vector <double> HB_kpp_cur;
   Double_t fac=299.8/(3./180.*3.14159);
    fit_V10_tosca_int = gr_HB_V10_tosca_int->GetFunction("pol2");
   for ( Int_t i=0;i<100;i++) {
     HB_cur.push_back(i*4000./100.);
     HB_mom.push_back(fac*HB_cur[i]*(fit_V10_tosca_int->Eval(HB_cur[i],0.,0.))/1000.);
     HB_kpp_cur.push_back(330.05*HB_mom[i]-3.1784*HB_mom[i]*HB_mom[i]+0.4018*HB_mom[i]*HB_mom[i]*HB_mom[i]);
   }
    TGraph *gr_HB_mom = new TGraph(HB_cur.size(),&(HB_mom[0]),&(HB_cur[0]));
    TGraph *gr_HB_cur = new TGraph(HB_cur.size(),&(HB_cur[0]),&(HB_mom[0]));
    TGraph *gr_KPP_mom = new TGraph(HB_cur.size(),&(HB_mom[0]),&(HB_kpp_cur[0]));
    TMultiGraph *mg_cur = new TMultiGraph();
 TLegend *leg_cur = new TLegend(.3,.2,.66,.4);
 TF1* HB_cur_mom_func;
 TF1* HB_mom_cur_func;
 TF1* HB_kpp_cur_mom_func;
 mg_cur->Add(gr_HB_cur);  
 leg_cur->AddEntry(gr_HB_mom,"Tosca HB-V10","p");   
 mg_cur->Add(gr_KPP_mom);  
 leg_cur->AddEntry(gr_KPP_mom,"KPP","p"); 
 TF1 *fit_pol = new TF1("fit_pol","[0]*x+[1]*x^2+[2]*x^3"); 
  fit_pol->SetParameters(329.,-2.5,0.31); 
  gr_HB_mom->Fit(fit_pol,"","",0.,11.);
 gr_HB_cur->Fit(fit_pol,"","",0.,4000.);
 gr_KPP_mom->Fit(fit_pol,"","",0.,4000.);
 gr_HB_cur->SetLineColor(2);
 HB_cur_mom_func = gr_HB_mom->GetFunction("fit_pol");
 HB_mom_cur_func = gr_HB_cur->GetFunction("fit_pol");
 HB_kpp_cur_mom_func = gr_KPP_mom->GetFunction("fit_pol");
 cout << " Mom at current = 972.25 = " << HB_mom_cur_func->Eval(972.25,0,0) << endl;
 cout << " Mom at current = 500 = " << HB_mom_cur_func->Eval(500.,0,0) << endl;
 cout << " KPP Current at mom = 3.0 = " << HB_kpp_cur_mom_func->Eval(3.,0,0) << endl;
 cout << " V10 Current at mom = 3.0 = " << HB_cur_mom_func->Eval(3.0,0,0) << endl;
 TCanvas *can1 = new TCanvas("can1","Momentum ",800,800);
   can1->Divide(1,2);
   can1->cd(1);
   gr_HB_cur->Draw();
   can1->cd(2);
   gr_HB_mom->Draw();
//
}
//
void plot_zscan_results() {
  TString dposfile="Central_zscan_files/central_zscan_pos_results.dat";
  TString dnegfile="Central_zscan_files/central_zscan_neg_results.dat";
  vector<double> Pos_tosca_int;
  vector<double> Pos_tosca_max;
  vector<double> Pos_tosca_effl;
  vector<double> Pos_data_int;
  vector<double> Pos_data_max;
  vector<double> Pos_data_effl;
  vector<double> Pos_rat_int;
  vector<double> Pos_rat_max;
  vector<double> Pos_rat_effl;
  vector<double> Pos_cur;
  Float_t col1,col2,col3,col4,col5;
   TTree *T = new TTree("ntuple","data");
   Long64_t nlines = T->ReadFile(dposfile,"col1:col2:col3:col4:col5",' ');
    T->SetBranchAddress("col1",&col1);
   T->SetBranchAddress("col2",&col2);
   T->SetBranchAddress("col3",&col3);
   T->SetBranchAddress("col4",&col4);
   T->SetBranchAddress("col5",&col5);
   for (Int_t i = 0; i <nlines;i++) {
     T->GetEntry(i);
     cout << col1 << " "  << col2/col1 << " "  << col3 << " "  << col4 << " "  << col5 << " " << endl; 
     Pos_cur.push_back(col1);
     Pos_data_max.push_back(col2/col1);
     Pos_data_int.push_back(col3/col1);
     Pos_data_effl.push_back(col3/col2*100);
     Pos_tosca_max.push_back(col4/col1);
     Pos_rat_max.push_back(col2/col4);
     Pos_tosca_int.push_back(col5/col1);
     Pos_rat_int.push_back(col3/col5);
     Pos_rat_effl.push_back(col3/col2/(col5/col4));
     Pos_tosca_effl.push_back(col5/col4*100);
   }
   T->Delete();
   Int_t colh;
   vector <double>  Neg_tosca_int;
   vector<double>  Neg_tosca_max;
   vector<double>  Neg_tosca_effl;
   vector<double>  Neg_data_int;
   vector<double>  Neg_data_max;
   vector<double>  Neg_data_effl;
   vector<double>  Neg_rat_int;
   vector<double>  Neg_rat_max;
   vector<double>  Neg_rat_effl;
   vector<double>  Neg_cur;
   TTree *T1 = new TTree("ntuple","data");
   nlines = T1->ReadFile(dnegfile,"col1:colh/I:col2/F:col3:col4:col5",' ');
    T1->SetBranchAddress("col1",&col1);
    T1->SetBranchAddress("colh",&colh);
   T1->SetBranchAddress("col2",&col2);
   T1->SetBranchAddress("col3",&col3);
   T1->SetBranchAddress("col4",&col4);
   T1->SetBranchAddress("col5",&col5);
   for (Int_t i = 0; i <nlines;i++) {
     T1->GetEntry(i);
     if (colh ==5) {
     cout << col1 << " "  << colh << " "  << col2/col1  << " "  << col3 << " "  << col4 << " "  << col5 << " " << endl; 
     Neg_cur.push_back(col1);
     Neg_data_max.push_back(TMath::Abs(col2/col1));
     Neg_data_int.push_back(TMath::Abs(col3/col1));
     Neg_data_effl.push_back(TMath::Abs(col3/col2)*100.);
     Neg_tosca_max.push_back(TMath::Abs(col4/col1));
     Neg_tosca_int.push_back(TMath::Abs(col5/col1));
     Neg_tosca_effl.push_back(TMath::Abs(col5/col4)*100.);
     Neg_rat_int.push_back(TMath::Abs(col3/col5));
     Neg_rat_max.push_back(TMath::Abs(col2/col4));
     Neg_rat_effl.push_back(TMath::Abs(col3/col2/(col5/col4)));
     }
   }
   T1->Delete();
//
   vector <double> HB_AB_cur;
   vector <double> HB_AB_tosca_int;
   vector<double>  HB_AB_tosca_max;
   vector<double>  HB_AB_tosca_effl;
   TTree *T2 = new TTree("ntuple","data");
   nlines = T2->ReadFile("Central_zscan_files/tosca_HB-AB-cases.dat","col1:col2:col3",' ');
    T2->SetBranchAddress("col1",&col1);
   T2->SetBranchAddress("col2",&col2);
   T2->SetBranchAddress("col3",&col3);
   for (Int_t i = 0; i <nlines;i++) {
     T2->GetEntry(i);
     HB_AB_cur.push_back(col1);
     HB_AB_tosca_max.push_back(TMath::Abs(col2/col1));
     HB_AB_tosca_int.push_back(TMath::Abs(col3/col1));
     HB_AB_tosca_effl.push_back(TMath::Abs(col3/col2)*100.);
     cout << col1 << " " << col2 << " " << col2/col1 << endl;
   }
    T2->Delete();
  // plot
  gROOT->Reset();
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.5,"Y");
 gStyle->SetTitleOffset(.8,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.1);
 TCanvas *can = new TCanvas("can","Int Bdl ",800,800);
   can->Divide(1,2);
    can->cd(1);
    TGraph *gr_pos_data_int = new TGraph(Pos_cur.size(),&(Pos_cur[0]),&(Pos_data_int[0]));
    TGraph *gr_pos_tosca_int = new TGraph(Pos_cur.size(),&(Pos_cur[0]),&(Pos_tosca_int[0]));
    TGraph *gr_pos_rat_int = new TGraph(Pos_cur.size(),&(Pos_cur[0]),&(Pos_rat_int[0]));
    TGraph *gr_neg_rat_int = new TGraph(Neg_cur.size(),&(Neg_cur[0]),&(Neg_rat_int[0]));
    TGraph *gr_neg_data_int = new TGraph(Neg_cur.size(),&(Neg_cur[0]),&(Neg_data_int[0]));
    TGraph *gr_neg_tosca_int = new TGraph(Neg_cur.size(),&(Neg_cur[0]),&(Neg_tosca_int[0]));
    TGraph *gr_HB_AB_tosca_int = new TGraph(HB_AB_cur.size(),&(HB_AB_cur[0]),&(HB_AB_tosca_int[0]));
    TMultiGraph *mg_int = new TMultiGraph();
 TLegend *leg_int = new TLegend(.3,.2,.66,.4);
    TMultiGraph *mg_int_rat = new TMultiGraph();
 TLegend *leg_int_rat = new TLegend(.3,.2,.66,.4);
 leg_int->AddEntry(gr_pos_data_int,"Data Pos ","p");
 leg_int->AddEntry(gr_pos_tosca_int,"Tosca Pos ","p");
 leg_int->AddEntry(gr_neg_data_int,"Data Neg ","p");
 leg_int->AddEntry(gr_neg_tosca_int,"Tosca Neg ","p");
 leg_int->AddEntry(gr_HB_AB_tosca_int,"Tosca HB-AB","p");
 mg_int->Add(gr_pos_data_int);  
 mg_int->Add(gr_neg_data_int);  
 mg_int->Add(gr_neg_tosca_int);  
 mg_int->Add(gr_pos_tosca_int);  
 mg_int->Add(gr_HB_AB_tosca_int);  
 TString htitle="HB Integral Bdl/Current";
 gr_neg_data_int->SetMarkerStyle(21);
 gr_neg_data_int->SetMarkerSize(1.5);
 gr_neg_data_int->SetMarkerColor(3);
 gr_neg_tosca_int->SetMarkerStyle(22);
 gr_neg_tosca_int->SetMarkerSize(1.5);
 gr_neg_tosca_int->SetMarkerColor(4);
 gr_pos_data_int->SetMarkerStyle(21);
 gr_pos_tosca_int->SetMarkerStyle(22);
 gr_HB_AB_tosca_int->SetMarkerStyle(22);
 gr_HB_AB_tosca_int->SetMarkerSize(1.5);
 gr_HB_AB_tosca_int->SetMarkerColor(6);
 mg_int->Draw("AP");
 mg_int->SetTitle(htitle);
 mg_int->GetXaxis()->SetTitle("Current");
 mg_int->GetYaxis()->SetTitle("Int Bdl/I (Tm/I)");
 leg_int->Draw();
    can->cd(2);
  mg_int_rat->Add(gr_pos_rat_int);  
  mg_int_rat->Add(gr_neg_rat_int);  
 htitle="HB Integral Bdl/Current Ratio data/tosca";
 leg_int_rat->AddEntry(gr_pos_rat_int,"Pos ","p");
 leg_int_rat->AddEntry(gr_neg_rat_int,"Neg ","p");
    gr_pos_rat_int->SetMarkerStyle(22);
    gr_pos_rat_int->SetMarkerSize(1.5);
    gr_neg_rat_int->SetMarkerStyle(21);
    gr_neg_rat_int->SetMarkerSize(1.5);
 mg_int_rat->Draw("AP");
 mg_int_rat->SetTitle(htitle);
 mg_int_rat->GetXaxis()->SetTitle("Current");
 mg_int_rat->GetYaxis()->SetTitle("Ratio Data/Tosca");
 leg_int_rat->Draw();
 //
 TCanvas *can_max = new TCanvas("can_max","Max B ",800,800);
   can_max->Divide(1,2);
    can_max->cd(1);
    TGraph *gr_pos_data_max = new TGraph(Pos_cur.size(),&(Pos_cur[0]),&(Pos_data_max[0]));
    TGraph *gr_pos_tosca_max = new TGraph(Pos_cur.size(),&(Pos_cur[0]),&(Pos_tosca_max[0]));
    TGraph *gr_pos_rat_max = new TGraph(Pos_cur.size(),&(Pos_cur[0]),&(Pos_rat_max[0]));
    TGraph *gr_neg_rat_max = new TGraph(Neg_cur.size(),&(Neg_cur[0]),&(Neg_rat_max[0]));
    TGraph *gr_neg_data_max = new TGraph(Neg_cur.size(),&(Neg_cur[0]),&(Neg_data_max[0]));
    TGraph *gr_neg_tosca_max = new TGraph(Neg_cur.size(),&(Neg_cur[0]),&(Neg_tosca_max[0]));
    TGraph *gr_HB_AB_tosca_max = new TGraph(HB_AB_cur.size(),&(HB_AB_cur[0]),&(HB_AB_tosca_max[0]));
    TMultiGraph *mg_max = new TMultiGraph();
 TLegend *leg_max = new TLegend(.3,.2,.66,.4);
    TMultiGraph *mg_max_rat = new TMultiGraph();
 TLegend *leg_max_rat = new TLegend(.3,.2,.66,.4);
 leg_max->AddEntry(gr_pos_data_max,"Data Pos ","p");
 leg_max->AddEntry(gr_pos_tosca_max,"Tosca Pos ","p");
 leg_max->AddEntry(gr_neg_data_max,"Data Neg ","p");
 leg_max->AddEntry(gr_neg_tosca_max,"Tosca Neg ","p");
 leg_max->AddEntry(gr_HB_AB_tosca_max,"Tosca HB-AB","p");
 mg_max->Add(gr_pos_data_max);  
 mg_max->Add(gr_neg_data_max);  
 mg_max->Add(gr_neg_tosca_max);  
 mg_max->Add(gr_pos_tosca_max);  
 mg_max->Add(gr_HB_AB_tosca_max);  
 htitle="HB  Max B/I";
 gr_neg_data_max->SetMarkerStyle(21);
 gr_neg_data_max->SetMarkerSize(1.5);
 gr_neg_data_max->SetMarkerColor(3);
 gr_neg_tosca_max->SetMarkerStyle(22);
 gr_neg_tosca_max->SetMarkerSize(1.5);
 gr_neg_tosca_max->SetMarkerColor(4);
 gr_pos_data_max->SetMarkerStyle(21);
 gr_pos_tosca_max->SetMarkerStyle(22);
 gr_pos_data_max->SetMarkerSize(1.5);
 gr_pos_tosca_max->SetMarkerSize(1.5);
 gr_pos_tosca_max->SetMarkerColor(2);
 gr_HB_AB_tosca_max->SetMarkerStyle(22);
 gr_HB_AB_tosca_max->SetMarkerSize(1.5);
 gr_HB_AB_tosca_max->SetMarkerColor(6);
 mg_max->Draw("AP");
 mg_max->SetTitle(htitle);
 mg_max->GetXaxis()->SetTitle("Current");
 mg_max->GetYaxis()->SetTitle("Max B/I (T/A)");
 leg_max->Draw();
    can_max->cd(2);
  mg_max_rat->Add(gr_pos_rat_max);  
  mg_max_rat->Add(gr_neg_rat_max);  
 htitle=" HB Max  B  Ratio data/tosca";
 leg_max_rat->AddEntry(gr_pos_rat_max,"Pos ","p");
 leg_max_rat->AddEntry(gr_neg_rat_max,"Neg ","p");
    gr_pos_rat_max->SetMarkerStyle(22);
    gr_pos_rat_max->SetMarkerSize(1.5);
    gr_neg_rat_max->SetMarkerStyle(21);
    gr_neg_rat_max->SetMarkerSize(1.5);
 mg_max_rat->Draw("AP");
 mg_max_rat->SetTitle(htitle);
 mg_max_rat->GetXaxis()->SetTitle("Current");
 mg_max_rat->GetYaxis()->SetTitle("Ratio Data/Tosca");
 leg_max_rat->Draw();
 //
 //
 TCanvas *can_effl = new TCanvas("can_effl","Eff length ",800,800);
   can_effl->Divide(1,2);
    can_effl->cd(1);
    TGraph *gr_pos_data_effl = new TGraph(Pos_cur.size(),&(Pos_cur[0]),&(Pos_data_effl[0]));
    TGraph *gr_pos_tosca_effl = new TGraph(Pos_cur.size(),&(Pos_cur[0]),&(Pos_tosca_effl[0]));
    TGraph *gr_pos_rat_effl = new TGraph(Pos_cur.size(),&(Pos_cur[0]),&(Pos_rat_effl[0]));
    TGraph *gr_neg_rat_effl = new TGraph(Neg_cur.size(),&(Neg_cur[0]),&(Neg_rat_effl[0]));
    TGraph *gr_neg_data_effl = new TGraph(Neg_cur.size(),&(Neg_cur[0]),&(Neg_data_effl[0]));
    TGraph *gr_neg_tosca_effl = new TGraph(Neg_cur.size(),&(Neg_cur[0]),&(Neg_tosca_effl[0]));
    TGraph *gr_HB_AB_tosca_effl = new TGraph(HB_AB_cur.size(),&(HB_AB_cur[0]),&(HB_AB_tosca_effl[0]));
    TMultiGraph *mg_effl = new TMultiGraph();
 TLegend *leg_effl = new TLegend(.3,.2,.66,.4);
    TMultiGraph *mg_effl_rat = new TMultiGraph();
 TLegend *leg_effl_rat = new TLegend(.3,.2,.66,.4);
 leg_effl->AddEntry(gr_pos_data_effl,"Data Pos ","p");
 leg_effl->AddEntry(gr_pos_tosca_effl,"Tosca Pos ","p");
 leg_effl->AddEntry(gr_neg_data_effl,"Data Neg ","p");
 leg_effl->AddEntry(gr_neg_tosca_effl,"Tosca Neg ","p");
 leg_effl->AddEntry(gr_HB_AB_tosca_effl,"Tosca HB-AB","p");
 mg_effl->Add(gr_pos_data_effl);  
 mg_effl->Add(gr_neg_data_effl);  
 mg_effl->Add(gr_neg_tosca_effl);  
 mg_effl->Add(gr_pos_tosca_effl);  
 mg_effl->Add(gr_HB_AB_tosca_effl);  
 htitle="HB  Eff Length (cm)";
 gr_neg_data_effl->SetMarkerStyle(21);
 gr_neg_data_effl->SetMarkerSize(1.5);
 gr_neg_data_effl->SetMarkerColor(3);
 gr_neg_tosca_effl->SetMarkerStyle(22);
 gr_neg_tosca_effl->SetMarkerSize(1.5);
 gr_neg_tosca_effl->SetMarkerColor(4);
 gr_pos_data_effl->SetMarkerStyle(21);
 gr_pos_tosca_effl->SetMarkerStyle(22);
 gr_pos_data_effl->SetMarkerSize(1.5);
 gr_pos_tosca_effl->SetMarkerSize(1.5);
 gr_pos_tosca_effl->SetMarkerColor(2);
 gr_HB_AB_tosca_effl->SetMarkerStyle(22);
 gr_HB_AB_tosca_effl->SetMarkerSize(1.5);
 gr_HB_AB_tosca_effl->SetMarkerColor(6);
 mg_effl->Draw("AP");
 mg_effl->SetTitle(htitle);
 mg_effl->GetXaxis()->SetTitle("Current");
 mg_effl->GetYaxis()->SetTitle("Eff Length (cm)");
 leg_effl->Draw();
    can_effl->cd(2);
  mg_effl_rat->Add(gr_pos_rat_effl);  
  mg_effl_rat->Add(gr_neg_rat_effl);  
 htitle="HB Eff length Ratio data/tosca";
 leg_effl_rat->AddEntry(gr_pos_rat_effl,"Pos ","p");
 leg_effl_rat->AddEntry(gr_neg_rat_effl,"Neg ","p");
    gr_pos_rat_effl->SetMarkerStyle(22);
    gr_pos_rat_effl->SetMarkerSize(1.5);
    gr_neg_rat_effl->SetMarkerStyle(21);
    gr_neg_rat_effl->SetMarkerSize(1.5);
 mg_effl_rat->Draw("AP");
 mg_effl_rat->SetTitle(htitle);
 mg_effl_rat->GetXaxis()->SetTitle("Current");
 mg_effl_rat->GetYaxis()->SetTitle("Ratio Data/Tosca");
 leg_effl_rat->Draw();
 //
}
//
void run_plot_zscan_pos(TString current) {
  if (probe_not_set) {
      Get_two_col_data_from_file("Lakeshore-probe.dat",nominal_field_corfile_data,corr_field_corfile_data);
      Set_probe(nominal_field_corfile_data,corr_field_corfile_data);
  }
  TString dfile;
  dfile="Central_zscan_files/central_zscan_"+current+"A_pos.dat";
  Float_t col1,col2,col3;
   TTree *T = new TTree("ntuple","data");
   Long64_t nlines = T->ReadFile(dfile,"col1:col2:col3",',');
   T->SetBranchAddress("col1",&col1);
   T->SetBranchAddress("col2",&col2);
   T->SetBranchAddress("col3",&col3);
   for (Int_t i = 0; i <nlines;i++) {
     T->GetEntry(i);
     if (col2!=0) {
       HB_central_zscan_data_dis.push_back(col1);
       HB_central_zscan_data.push_back(col2);
     }
     if (col3!=0) {
       HB_central_zscan_tosca_dis.push_back(col1);
       HB_central_zscan_tosca.push_back(col3);
     }
   }
   T->Delete();
  Plot_zscan_pos(current,HB_central_zscan_data_dis,HB_central_zscan_data,HB_central_zscan_tosca_dis,HB_central_zscan_tosca);
  HB_central_zscan_data.clear();
  HB_central_zscan_data_dis.clear();
  HB_central_zscan_tosca.clear();
  HB_central_zscan_tosca_dis.clear();
}
//
void Plot_zscan_pos(TString cur,vector<double> &zpos_data ,vector<double> &uncorr_field_data ,vector<double> &zpos_tosca ,vector<double> &field_tosca ) {
 gROOT->Reset();
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.5,"Y");
 gStyle->SetTitleOffset(.8,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.05);
 vector<double> corr_field_data;
 vector<double> corr_field_tosca_ratio;
 Double_t data_int=0;
 Double_t tosca_int=0;
 Double_t data_max=0;
 Double_t tosca_max=0;
 for (UInt_t i=0; i<uncorr_field_data.size(); i++) {
   corr_field_data.push_back(uncorr_field_data[i]/10. + field_corr(uncorr_field_data[i]/10.));
   if ( TMath::Abs(corr_field_data[i]) > TMath::Abs(data_max)) data_max =corr_field_data[i];
   //    cout << zpos[i] << " " << corr_field_data[i] << " " << uncorr_field_data[i]/10. << " " <<  corr_field_data[i]/field_tosca[i] << " " << uncorr_field_data[i]/field_tosca[i]<< endl;
 }
 for (UInt_t i=0; i<corr_field_data.size()-1; i++) {
   data_int+=(corr_field_data[i+1]+corr_field_data[i])/2.*TMath::Abs(zpos_data[i+1]-zpos_data[i])/100.;
 }
 for (UInt_t i=0; i<zpos_tosca.size(); i++) {
   field_tosca[i] = field_tosca[i]/10.;
   if ( TMath::Abs(field_tosca[i]) > TMath::Abs(tosca_max)) tosca_max =field_tosca[i];
 }
 for (UInt_t i=0; i<zpos_tosca.size()-1; i++) {
   tosca_int+=(field_tosca[i+1]+field_tosca[i])/2*TMath::Abs(zpos_tosca[i+1]-zpos_tosca[i])/100.;
 }
 for (UInt_t i=0; i<zpos_tosca.size(); i++) {
  for (UInt_t j=0; j<zpos_data.size(); j++) {
    if (zpos_tosca[i]==zpos_data[j]) corr_field_tosca_ratio.push_back((1-corr_field_data[j]/field_tosca[i])*100.);
  }
 }
cout << " Max data = " << data_max << " Max Tosca = " << tosca_max << endl; 
cout << " Int data = " << data_int << " Max Tosca = " << tosca_int << endl; 
cout << " Eff L data = " << 100*data_int/data_max << " Eff L Tosca = " << 100*tosca_int/tosca_max << endl; 
 cout << cur.Atof() << " " << data_max << " " << data_int << " " << 100*data_int/data_max << " " << tosca_max << " " << tosca_int << " " << 100*tosca_int/tosca_max << " " << data_max/tosca_max << " " << data_int/tosca_int << " " << data_int/tosca_int/(data_max/tosca_max) << endl;
 Double_t effL_data=100*data_int/data_max;
 Double_t effL_tosca=100*tosca_int/tosca_max;
  TString ofile;
  ofile="Plotfile/central_zscan_"+cur+"A_pos.pdf";
 TCanvas *cscan = new TCanvas("cscan"," Central zcan ",800,800);
 cscan->Divide(1,2);
 cscan->cd(1);
 TGraph *gr_data = new TGraph(corr_field_data.size(),&(zpos_data[0]),&(corr_field_data[0]));
 TGraph *gr_tosca = new TGraph(field_tosca.size(),&(zpos_tosca[0]),&(field_tosca[0]));
 TGraph *gr_ratio = new TGraph(corr_field_tosca_ratio.size(),&(zpos_data[0]),&(corr_field_tosca_ratio[0]));
 TMultiGraph *mg_scan = new TMultiGraph();
 TLegend *leg_scan = new TLegend(.3,.2,.66,.4);
 leg_scan->AddEntry(gr_data,Form("Data max = %5.4f T, Int = %5.4f Tm",data_max,data_int),"p");
 leg_scan->AddEntry(gr_tosca,Form("Tosca max = %5.4f T, Int = %5.4f Tm",tosca_max,tosca_int),"p");
 gr_data->SetMarkerStyle(21);
 gr_tosca->SetMarkerStyle(22);
 gr_data->SetMarkerSize(.5);
 gr_tosca->SetMarkerSize(.5);
 gr_tosca->SetMarkerColor(2);
 mg_scan->Add(gr_data);
 mg_scan->Add(gr_tosca);
 mg_scan->Draw("AP");
 TString htitle=" Current = "+cur+" A ";
 mg_scan->SetTitle(htitle);
 mg_scan->GetXaxis()->SetTitle("Z pos (cm)");
 mg_scan->GetYaxis()->SetTitle("Field (T)");
 leg_scan->SetTextSize(.035);
 leg_scan->Draw();
 cscan->cd(2);
 gr_ratio->Draw("AP");
 gr_ratio->SetMarkerStyle(22);
 gr_ratio->SetTitle(htitle);
 gr_ratio->SetMarkerSize(.5);
 gr_ratio->GetXaxis()->SetTitle("Z pos (cm)");
 gr_ratio->GetYaxis()->SetTitle("Measured Field/ Tosca model (%)");
} 
//
//
void Plot_zscan(TString cur, TString nh,vector<double> &uncorr_field_data ,vector<double> &field_tosca ) {
 gROOT->Reset();
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.5,"Y");
 gStyle->SetTitleOffset(.8,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.05);
 vector<double> corr_field_data;
 vector<double> corr_field_tosca_ratio;
 vector<double> zpos; // data has 57 positions with pos = 27 at z=-3.413
 Double_t data_int=0;
 Double_t tosca_int=0;
 Double_t data_max=0;
 Double_t tosca_max=0;
 Int_t step;
 for (UInt_t i=0; i<uncorr_field_data.size(); i++) {
   zpos.push_back(2.54*((i+1)-27.)-3.413);
   corr_field_data.push_back(uncorr_field_data[i]/10. + field_corr(uncorr_field_data[i]/10.));
   field_tosca[i] = field_tosca[i]/10.;
   corr_field_tosca_ratio.push_back((1-corr_field_data[i]/field_tosca[i])*100.);
   data_int+=corr_field_data[i]*2.54/100.;
   tosca_int+=field_tosca[i]*2.54/100.;
   step = i+1-25;
   if ( TMath::Abs(corr_field_data[i]) > TMath::Abs(data_max)) data_max =corr_field_data[i];
   if ( TMath::Abs(field_tosca[i]) > TMath::Abs(tosca_max)) tosca_max =field_tosca[i];
   //    cout << zpos[i] << " " << corr_field_data[i] << " " << uncorr_field_data[i]/10. << " " <<  corr_field_data[i]/field_tosca[i] << " " << uncorr_field_data[i]/field_tosca[i]<< endl;
 }
cout << " Max data = " << data_max << " Max Tosca = " << tosca_max << endl; 
cout << " Int data = " << data_int << " Max Tosca = " << tosca_int << endl; 
cout << " Eff L data = " << 100*data_int/data_max << " Eff L Tosca = " << 100*tosca_int/tosca_max << endl; 
 cout << cur.Atof() << " " << nh.Atof() << " " << data_max << " " << data_int << " " << 100*data_int/data_max << " " << tosca_max << " " << tosca_int << " " << 100*tosca_int/tosca_max << " " << data_max/tosca_max << " " << data_int/tosca_int << " " << data_int/tosca_int/(data_max/tosca_max) << endl;
 Double_t effL_data=100*data_int/data_max;
 Double_t effL_tosca=100*tosca_int/tosca_max;
  TString ofile;
  ofile="Plotfile/central_zscan_"+cur+"a_hole"+nh+".pdf";
 TCanvas *cscan = new TCanvas("cscan"," Central zcan ",800,800);
 cscan->Divide(1,2);
 cscan->cd(1);
 TGraph *gr_data = new TGraph(corr_field_data.size(),&(zpos[0]),&(corr_field_data[0]));
 TGraph *gr_tosca = new TGraph(field_tosca.size(),&(zpos[0]),&(field_tosca[0]));
 TGraph *gr_ratio = new TGraph(corr_field_tosca_ratio.size(),&(zpos[0]),&(corr_field_tosca_ratio[0]));
 TMultiGraph *mg_scan = new TMultiGraph();
 TLegend *leg_scan = new TLegend(.3,.7,.66,.9);
 leg_scan->AddEntry(gr_data,Form("Data max = %5.4f T, Int = %5.4f Tm",data_max,data_int),"p");
 leg_scan->AddEntry(gr_tosca,Form("Tosca max = %5.4f T, Int = %5.4f Tm",tosca_max,tosca_int),"p");
 gr_data->SetMarkerStyle(21);
 gr_tosca->SetMarkerStyle(22);
 gr_data->SetMarkerSize(.5);
 gr_tosca->SetMarkerSize(.5);
 gr_tosca->SetMarkerColor(2);
 mg_scan->Add(gr_data);
 mg_scan->Add(gr_tosca);
 mg_scan->Draw("AP");
 TString htitle=" Current = "+cur+" A , hole  = "+nh;
 mg_scan->SetTitle(htitle);
 mg_scan->GetXaxis()->SetTitle("Z pos (cm)");
 mg_scan->GetYaxis()->SetTitle("Field (T)");
 leg_scan->SetTextSize(.035);
 leg_scan->Draw();
 cscan->cd(2);
 gr_ratio->Draw("AP");
 gr_ratio->SetMarkerStyle(22);
 gr_ratio->SetTitle(htitle);
 gr_ratio->SetMarkerSize(.5);
 gr_ratio->GetXaxis()->SetTitle("Z pos (cm)");
 gr_ratio->GetYaxis()->SetTitle("Measured Field/ Tosca model (%)");
} 
//
void run_code()
{
  
  Get_two_col_data_from_file("Lakeshore-probe.dat",nominal_field_corfile_data,corr_field_corfile_data);
  Set_probe(nominal_field_corfile_data,corr_field_corfile_data);
  Get_two_col_data_from_file("HB-central-B-I.dat",HB_central_current_data,HB_central_field_data);
  Plot_hb_bi(HB_central_current_data,HB_central_field_data);
}
//
void Plot_hb_pos_neg_diff() {
vector<double> HB_neg_down_current;
vector<double> HB_neg_down_field;
vector<double> HB_neg_down_bi;
vector<double> HB_neg_down_resid;
vector<double> HB_neg_up_current;
vector<double> HB_neg_up_field;
vector<double> HB_neg_up_bi;
vector<double> HB_neg_up_resid;
vector<double> HB_negpos_up_diff;
vector<double> HB_negpos_up_diffper;
vector<double> HB_negpos_down_diff;
vector<double> HB_negpos_down_diffper;
vector<double> HB_pos_up_current;
vector<double> HB_pos_up_field;
vector<double> HB_pos_up_bi;
vector<double> HB_pos_up_resid;
vector<double> HB_pos_down_current;
vector<double> HB_pos_down_field;
vector<double> HB_pos_down_bi;
vector<double> HB_pos_down_resid;
vector<double> HB_pos_up_down_diff;
vector<double> HB_pos_up_down_diffper;
vector<double> HB_neg_up_down_diff;
vector<double> HB_neg_up_down_diffper;
TF1 *pos_field_up_fit;
TF1 *pos_field_down_fit;
TF1 *neg_field_up_fit;
TF1 *neg_field_down_fit;
 gROOT->Reset();
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.12);
 //
  Get_two_col_data_from_file("Lakeshore-probe.dat",nominal_field_corfile_data,corr_field_corfile_data);
   Set_probe(nominal_field_corfile_data,corr_field_corfile_data);
   Get_two_col_data_from_file("neg-down-current.dat",HB_neg_down_field,HB_neg_down_current);
   Get_two_col_data_from_file("neg-up-current.dat",HB_neg_up_field,HB_neg_up_current);
   Get_two_col_data_from_file("pos-up-current.dat",HB_pos_up_field,HB_pos_up_current);
   Get_two_col_data_from_file("pos-down-current.dat",HB_pos_down_field,HB_pos_down_current);
   //
 for (UInt_t i=0; i<HB_neg_up_current.size(); i++) {
   HB_neg_up_field[i]=HB_neg_up_field[i]/10. + field_corr(HB_neg_up_field[i]/10.);
   HB_neg_up_bi.push_back(HB_neg_up_field[i]/HB_neg_up_current[i]);
   //    cout << HB_neg_up_field[i] << " " << HB_neg_up_current[i] << endl;
}
 for (UInt_t i=0; i<HB_neg_down_current.size(); i++) {
   HB_neg_down_field[i]=HB_neg_down_field[i]/10. + field_corr(HB_neg_down_field[i]/10.);
   HB_neg_down_bi.push_back(HB_neg_down_field[i]/HB_neg_down_current[i]);
 }
 for (UInt_t i=0; i<HB_pos_up_current.size(); i++) {
   HB_pos_up_field[i]=HB_pos_up_field[i]/10. + field_corr(HB_pos_up_field[i]/10.);
   HB_pos_up_bi.push_back(HB_pos_up_field[i]/HB_pos_up_current[i]);
 }
 for (UInt_t i=0; i<HB_pos_down_current.size(); i++) {
   HB_pos_down_field[i]=HB_pos_down_field[i]/10. + field_corr(HB_pos_down_field[i]/10.);
   
   if (HB_pos_down_current[i] !=0)  HB_pos_down_bi.push_back(HB_pos_down_field[i]/HB_pos_down_current[i]);
 }
 for (UInt_t i=0; i<HB_pos_up_current.size(); i++) {
   HB_negpos_up_diff.push_back(10000*(HB_pos_up_field[i]-TMath::Abs(HB_neg_up_field[i])));
   HB_negpos_up_diffper.push_back(100.*TMath::Abs((HB_pos_up_field[i]-TMath::Abs(HB_neg_up_field[i]))/HB_pos_up_field[i]));
 }
 for (UInt_t i=0; i<HB_pos_down_current.size(); i++) {
   HB_negpos_down_diff.push_back(10000*(HB_pos_down_field[i]-TMath::Abs(HB_neg_down_field[i])));
   HB_negpos_down_diffper.push_back(100.*TMath::Abs((HB_pos_down_field[i]-TMath::Abs(HB_neg_down_field[i]))/HB_pos_down_field[i]));
 }
 for (UInt_t i=0; i<HB_pos_up_current.size(); i++) {
   HB_pos_up_down_diff.push_back(10000*(HB_pos_up_field[i]-HB_pos_down_field[i]));
   HB_pos_up_down_diffper.push_back(100*HB_pos_up_down_diff[i]/(10000*HB_pos_up_field[i]));
 }
 for (UInt_t i=0; i<HB_neg_up_current.size(); i++) {
   HB_neg_up_down_diff.push_back(10000*(HB_neg_up_field[i]-HB_neg_down_field[i]));
   HB_neg_up_down_diffper.push_back(TMath::Abs(100*HB_neg_up_down_diff[i]/(10000*HB_neg_up_field[i])));
 }
 //
TCanvas *cupdiff = new TCanvas("cupdiff","HB Ramp up Pos-Neg Field",800,800);
 cupdiff->Divide(1,2);
 cupdiff->cd(1);
 TGraph *gr_updiff = new TGraph(HB_neg_up_current.size(),&(HB_pos_up_current[0]),&(HB_negpos_up_diff[0]));
 gr_updiff->Draw("AP");
 gr_updiff->SetMarkerStyle(21);
 gr_updiff->SetTitle("HB Ramping Up Difference between Pos and Neg  Field");
 gr_updiff->GetXaxis()->SetTitle("Current (A)");
 gr_updiff->GetYaxis()->SetTitle("Field Difference (G)");
 cupdiff->cd(2);
 TGraph *gr_updiffper = new TGraph(HB_neg_up_current.size(),&(HB_pos_up_current[0]),&(HB_negpos_up_diffper[0]));
 gr_updiffper->Draw("AP");
 gr_updiffper->SetMarkerStyle(21);
 gr_updiffper->SetTitle("HB Ramping Up Difference between Pos and Neg  Field");
 gr_updiffper->GetXaxis()->SetTitle("Current (A)");
 gr_updiffper->GetYaxis()->SetTitle("Field Difference (%)");
 gr_updiffper->SetMaximum(.25);
 gr_updiffper->SetMinimum(.0);
 //
 //
TCanvas *cdowndiff = new TCanvas("cdowndiff","HB Ramp down Pos-Neg Field",800,800);
 cdowndiff->Divide(1,2);
 cdowndiff->cd(1);
 TGraph *gr_downdiff = new TGraph(HB_neg_down_current.size(),&(HB_pos_down_current[0]),&(HB_negpos_down_diff[0]));
 gr_downdiff->Draw("AP");
 gr_downdiff->SetMarkerStyle(21);
 gr_downdiff->SetTitle("HB Ramping Down Difference between Pos and Neg  Field");
 gr_downdiff->GetXaxis()->SetTitle("Current (A)");
 gr_downdiff->GetYaxis()->SetTitle("Field Difference (G)");
 cdowndiff->cd(2);
 TGraph *gr_downdiffper = new TGraph(HB_neg_down_current.size(),&(HB_pos_down_current[0]),&(HB_negpos_down_diffper[0]));
 gr_downdiffper->Draw("AP");
 gr_downdiffper->SetMarkerStyle(21);
 gr_downdiffper->SetTitle("HB Ramping Down Difference between Pos and Neg  Field");
 gr_downdiffper->GetXaxis()->SetTitle("Current (A)");
 gr_downdiffper->GetYaxis()->SetTitle("Field Difference (%)");
 gr_downdiffper->SetMaximum(.25);
 gr_downdiffper->SetMinimum(.0);
 //
TCanvas *cneg = new TCanvas("cneg","HB Neg B/I versus I",800,800);
 cneg->Divide(1,2);
 cneg->cd(1);
 TGraph *grnep_up = new TGraph(HB_neg_up_current.size(),&(HB_neg_up_current[0]),&(HB_neg_up_bi[0]));
 grnep_up->Draw("AP");
 grnep_up->SetMarkerStyle(21);
 grnep_up->SetTitle("HB Neg Field/Current versus Current");
 grnep_up->GetXaxis()->SetTitle("Current (A)");
 grnep_up->GetYaxis()->SetTitle("Field/Current (T/A)");
 grnep_up->SetMarkerSize(.5);
 grnep_up->SetMaximum(.00070);
 grnep_up->SetMinimum(.00060);
 grnep_up->Fit("pol3","","",-4000.,-100.);
 neg_field_up_fit = grnep_up->GetFunction("pol3");
 vector<double> HB_neg_up_frac;
 Double_t par1=neg_field_up_fit->GetParameter(0);
 for (UInt_t i=0;i<HB_neg_up_current.size();i++) {
   HB_neg_up_resid.push_back((neg_field_up_fit->Eval(HB_neg_up_current[i],0.,0.)-HB_neg_up_bi[i])/HB_neg_up_bi[i]);
   HB_neg_up_frac.push_back(100.*(1.-HB_neg_up_bi[i]/par1));
 }
 cneg->cd(2);
  TGraph *gr_neg_frac = new TGraph(HB_neg_up_current.size(),&(HB_neg_up_current[0]),&(HB_neg_up_frac[0]));
 gr_neg_frac->Draw("AP");
 gr_neg_frac->SetMarkerStyle(21);
 gr_neg_frac->SetTitle("HB Neg ormalized Field/Current versus Current");
 gr_neg_frac->GetXaxis()->SetTitle("Current (A)");
 gr_neg_frac->GetYaxis()->SetTitle("1-Field/Current/p0 (%)");
 gr_neg_frac->SetMarkerSize(.5);
 gr_neg_frac->SetMaximum(5.);
 gr_neg_frac->SetMinimum(0.);
 /* TGraph *gr_neg_up_resid = new TGraph(HB_neg_up_current.size(),&(HB_neg_up_current[0]),&(HB_neg_up_resid[0]));
 gr_neg_up_resid->Draw("AP");
 gr_neg_up_resid->SetTitle("HB Neg B/I FIT Residual versus Current");
 gr_neg_up_resid->GetXaxis()->SetTitle("Current (A)");
 gr_neg_up_resid->GetYaxis()->SetTitle("[B/I(FIT) - B/I]/ B/I");
 gr_neg_up_resid->SetMarkerStyle(21);
 gr_neg_up_resid->SetMarkerSize(.5);
 gr_neg_up_resid->SetMaximum(.002);
 gr_neg_up_resid->SetMinimum(-.002);
 */
  //
TCanvas *cpos = new TCanvas("cpos","HB Pos B/I versus I",800,800);
 cpos->Divide(1,2);
 cpos->cd(1);
 TGraph *grposup = new TGraph(HB_pos_up_current.size(),&(HB_pos_up_current[0]),&(HB_pos_up_bi[0]));
 grposup->Draw("AP");
 grposup->SetMarkerStyle(21);
 grposup->SetTitle("HB Pos Field/Current versus Current");
 grposup->GetXaxis()->SetTitle("Current (A)");
 grposup->GetYaxis()->SetTitle("Field/Current (T/A)");
 grposup->SetMarkerSize(.5);
 grposup->SetMaximum(.00070);
 grposup->SetMinimum(.00060);
 grposup->Fit("pol3","","",100.,4000.);
 pos_field_up_fit = grposup->GetFunction("pol3");
 vector<double> HB_pos_up_frac;
 Double_t parpos=pos_field_up_fit->GetParameter(0);
 for (UInt_t i=0;i<HB_pos_up_current.size();i++) {
   //cout << pos_field_up_fit->Eval(pos_up_cur[i],0.,0.) << " " << pos_up_field[i] << endl;
   HB_pos_up_resid.push_back((pos_field_up_fit->Eval(HB_pos_up_current[i],0.,0.)-HB_pos_up_bi[i])/HB_pos_up_bi[i]);
   HB_pos_up_frac.push_back(100.*(1.-HB_pos_up_bi[i]/parpos));
 }
 cpos->cd(2);
  TGraph *gr_pos_frac = new TGraph(HB_pos_up_current.size(),&(HB_pos_up_current[0]),&(HB_pos_up_frac[0]));
 gr_pos_frac->Draw("AP");
 gr_pos_frac->SetMarkerStyle(21);
 gr_pos_frac->SetTitle("HB Pos ormalized Field/Current versus Current");
 gr_pos_frac->GetXaxis()->SetTitle("Current (A)");
 gr_pos_frac->GetYaxis()->SetTitle("1-Field/Current/p0 (%)");
 gr_pos_frac->SetMarkerSize(.5);
 gr_pos_frac->SetMaximum(5.);
 gr_pos_frac->SetMinimum(0.);
 /*
 TGraph *gr_pos_up_resid = new TGraph(HB_pos_up_current.size(),&(HB_pos_up_current[0]),&(HB_pos_up_resid[0]));
 gr_pos_up_resid->Draw("AP");
 gr_pos_up_resid->SetTitle("HB Pos B/I FIT Residual versus Current");
 gr_pos_up_resid->GetXaxis()->SetTitle("Current (A)");
 gr_pos_up_resid->GetYaxis()->SetTitle("[B/I(FIT) - B/I]/ B/I");
 gr_pos_up_resid->SetMarkerStyle(21);
 gr_pos_up_resid->SetMarkerSize(.5);
 gr_pos_up_resid->SetMaximum(.002);
 gr_pos_up_resid->SetMinimum(-.002);
 */
//
 TGraph *grposup2 = new TGraph(HB_pos_up_current.size(),&(HB_pos_up_current[0]),&(HB_pos_up_bi[0]));
 TGraph *grposdown = new TGraph(HB_pos_down_current.size(),&(HB_pos_down_current[0]),&(HB_pos_down_bi[0]));
 TGraph *grposupdowndiff = new TGraph(HB_pos_down_field.size(),&(HB_pos_down_field[0]),&(HB_pos_up_down_diff[0]));
 TGraph *grposupdowndiffper = new TGraph(HB_pos_down_field.size(),&(HB_pos_down_field[0]),&(HB_pos_up_down_diffper[0]));
  TMultiGraph *mg1 = new TMultiGraph();
TCanvas *cposupdown = new TCanvas("cposupdown","HB Pos Up/Down",800,800);
 cposupdown->Divide(1,3);
 cposupdown->cd(1);
 grposup2->SetMarkerStyle(21);
 grposup2->SetMarkerSize(1.5);
 grposdown->SetMarkerStyle(22);
 grposdown->SetMarkerSize(1.5);
 grposdown->SetMarkerColor(2);
 mg1->Add(grposup2);
 mg1->Add(grposdown);
 mg1->Draw("AP");
 mg1->SetMaximum(.0007);
 mg1->SetMinimum(0.0006);
 mg1->SetTitle("HB Pos Field/Current versus Current");
 mg1->GetXaxis()->SetTitle("Current (A)");
 mg1->GetYaxis()->SetTitle("Field/Current (T/A)");
 cposupdown->cd(2);
 grposupdowndiff->SetMarkerStyle(22);
 grposupdowndiff->SetMarkerSize(1.5);
 grposupdowndiff->SetMarkerColor(1);
 grposupdowndiff->Draw("AP");
 grposupdowndiff->SetMaximum(0);
 grposupdowndiff->SetMinimum(-20.);
 grposupdowndiff->SetTitle("HB Pos (Up -Down) versus Field");
 grposupdowndiff->GetXaxis()->SetTitle("Field (T)");
 grposupdowndiff->GetYaxis()->SetTitle("(Up-Down) (G)");
 cposupdown->cd(3);
 grposupdowndiffper->SetMarkerStyle(22);
 grposupdowndiffper->SetMarkerSize(1.5);
 grposupdowndiffper->SetMarkerColor(1);
 grposupdowndiffper->Draw("AP");
 grposupdowndiffper->SetMaximum(0.);
 grposupdowndiffper->SetMinimum(-.5);
 grposupdowndiffper->SetTitle("HB Pos (Up -Down)/Down versus Field");
 grposupdowndiffper->GetXaxis()->SetTitle("Field (T)");
 grposupdowndiffper->GetYaxis()->SetTitle("(Up-Down)/(Down) (%)");
 //
//
 TGraph *grnegup2 = new TGraph(HB_neg_up_current.size(),&(HB_neg_up_current[0]),&(HB_neg_up_bi[0]));
 TGraph *grnegdown = new TGraph(HB_neg_down_current.size(),&(HB_neg_down_current[0]),&(HB_neg_down_bi[0]));
 TGraph *grnegupdowndiff = new TGraph(HB_neg_down_field.size(),&(HB_neg_down_field[0]),&(HB_neg_up_down_diff[0]));
 TGraph *grnegupdowndiffper = new TGraph(HB_neg_down_field.size(),&(HB_neg_down_field[0]),&(HB_neg_up_down_diffper[0]));
  TMultiGraph *mg_neg = new TMultiGraph();
TCanvas *cnegupdown = new TCanvas("cnegupdown","HB Neg Up/Down",800,800);
 cnegupdown->Divide(1,3);
 cnegupdown->cd(1);
 grnegup2->SetMarkerStyle(21);
 grnegup2->SetMarkerSize(1.5);
 grnegdown->SetMarkerStyle(22);
 grnegdown->SetMarkerSize(1.5);
 grnegdown->SetMarkerColor(2);
 mg_neg->Add(grnegup2);
 mg_neg->Add(grnegdown);
 mg_neg->Draw("AP");
 mg_neg->SetMaximum(.0007);
 mg_neg->SetMinimum(0.0006);
 mg_neg->SetTitle("HB Neg Field/Current versus Current");
 mg_neg->GetXaxis()->SetTitle("Current (A)");
 mg_neg->GetYaxis()->SetTitle("Field/Current (T/A)");
 cnegupdown->cd(2);
 grnegupdowndiff->SetMarkerStyle(22);
 grnegupdowndiff->SetMarkerSize(1.5);
 grnegupdowndiff->SetMarkerColor(1);
 grnegupdowndiff->Draw("AP");
 grnegupdowndiff->SetMaximum(20.);
 grnegupdowndiff->SetMinimum(0.);
 grnegupdowndiff->SetTitle("HB Neg (Up -Down) versus Field");
 grnegupdowndiff->GetXaxis()->SetTitle("Field (T)");
 grnegupdowndiff->GetYaxis()->SetTitle("(Up-Down) (G)");
 cnegupdown->cd(3);
 grnegupdowndiffper->SetMarkerStyle(22);
 grnegupdowndiffper->SetMarkerSize(1.5);
 grnegupdowndiffper->SetMarkerColor(1);
 grnegupdowndiffper->Draw("AP");
  grnegupdowndiffper->SetMaximum(0.5);
 grnegupdowndiffper->SetMinimum(0);
 grnegupdowndiffper->SetTitle("HB Neg (Up -Down)/Down versus Field");
 grnegupdowndiffper->GetXaxis()->SetTitle("Field (T)");
 grnegupdowndiffper->GetYaxis()->SetTitle("(Up-Down)/(Down) (%)");
}
//
Double_t field_corr(Double_t field )
{
  Double_t fieldcorr;
  fieldcorr=0.;
   if (field<0) fieldcorr=neg_field_corr->Eval(field,0.,0.);
   if (field >0) fieldcorr=-pos_field_corr->Eval(field,0.,0.);
  return fieldcorr;
}
//
void Get_two_col_data_from_file(TString filename,vector<double> &x,vector<double> &y )
{
 ifstream corfile;
 corfile.open(filename);
 TString curline;
 Int_t nline=0;
 while (corfile.good())
   {
     curline.ReadLine(corfile,kTRUE);
    TString sc=curline.Data();
    Int_t chi,ncomma=sc.CountChar(',');
    if (ncomma ==1) {
        chi=sc.Index(',');
        TString temp(sc(0,chi));
	x.push_back(temp.Atof());
        temp=sc(chi+1,sc.Sizeof());
	y.push_back(temp.Atof());
	//  cout << " data " << x[nline] << " "<< y[nline] <<endl;
         nline++;
    }
   }
}
void Plot_hb_bi(vector<double> &x ,vector<double> &y ) {
 gROOT->Reset();
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.5,"Y");
 gStyle->SetTitleOffset(.8,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.2);
 //
 Int_t vecsize=x.size();
 const Int_t npar=1000;
 Double_t xvec[npar];
 Double_t yvec[npar];
 Double_t bvec[npar];
 Double_t resid[npar];
 for (UInt_t i=0; i<x.size(); i++) {
   xvec[i]=x[i];
   yvec[i]=y[i]/10.;
   bvec[i]=yvec[i]+field_corr(yvec[i]);
   cout << " dataos " << xvec[i] << " "<< yvec[i] << " " << field_corr(yvec[i]) << " " << field_corr(yvec[i])/yvec[i] <<endl;
   yvec[i]=bvec[i]/x[i];
 }
 //
  const int nftot=5;
  Double_t tosca_cur[nftot]={1200.,2000.,3000.,3500,4000.};//A
  Double_t tosca_field[nftot]={8.01,13.342,19.879,22.933,25.771}; //kG
  Double_t tosca_curratio[nftot];
  for (int nf=0;nf<nftot;nf++) {
    tosca_curratio[nf]=tosca_field[nf]/tosca_cur[nf]/10.;// T/A
     for (UInt_t i=0; i<x.size(); i++) {
       if (abs(tosca_cur[nf]-xvec[i]) < 10.) {
	 cout << xvec[i] << " " << yvec[i] << " " << tosca_curratio[nf] << " " << yvec[i]/tosca_curratio[nf] << endl;
       }
        }
  }
 //
TCanvas *can4 = new TCanvas("can4","HB B/I versus I",800,800);
 can4->Divide(1,1);
 can4->cd(1);
 TGraph *grcur = new TGraph(vecsize,bvec,xvec);
 grcur->SetMarkerStyle(21);
 grcur->SetTitle("HB I versus B");
 grcur->GetYaxis()->SetTitle("Current [A]");
 grcur->GetXaxis()->SetTitle("B [T]");
 grcur->SetMarkerSize(1.0);
 grcur->SetMaximum(4000);
 grcur->SetMinimum(0);
 grcur->Draw("AP");
   grcur->Fit("pol3","QO","",0.,2.5);
  cur_field = grcur->GetFunction("pol3");
  Double_t ftest=0.051*6400/299.8/.78;
  cout << " current for b field = " << ftest << " " << cur_field->Eval(ftest,0.,0.) << endl;
   //
TCanvas *can = new TCanvas("can","HB B/I versus I",800,800);
 can->Divide(1,1);
 can->cd(1);
 TGraph *grposup = new TGraph(vecsize,xvec,yvec);

 grposup->SetName("grpopsup");
 grposup->Draw("AP");
 grposup->Draw("AP");
 grposup->SetMarkerStyle(21);
 grposup->SetTitle("HB B/I versus I");
 grposup->GetXaxis()->SetTitle("Current [A]");
 grposup->GetYaxis()->SetTitle("B/I [T/A]");
 grposup->SetMarkerSize(1.0);
 grposup->SetMaximum(.00070);
 grposup->SetMinimum(.00060);
  grposup->Fit("pol3","QO","",100.,4000.);
  TF1 *fit_data;
  fit_data = grposup->GetFunction("pol3");
 for (int i=0;i<vecsize;i++) {
   //cout << neg_field_corr>Eval(x[i],0.,0.) << " " << y[i] << endl;
      resid[i]=(fit_data->Eval(xvec[i],0.,0.)-yvec[i])/yvec[i];
}
 /*
TCanvas *can2 = new TCanvas("can2","HB B/I versus I",800,800);
 can2->cd(2);
 TGraph *gr2 = new TGraph(vecsize,xvec,resid);
 gr2->Draw("AP");
 gr2->SetTitle("HB B/I FIT Residual versus Current");
 gr2->GetXaxis()->SetTitle("Current (A)");
 gr2->GetYaxis()->SetTitle("[B/I(FIT) - B/I]/ B/I");
 gr2->SetMarkerStyle(21);
 gr2->SetMarkerSize(.5);
 gr2->SetMaximum(.002);
 gr2->SetMinimum(-.002);
 */
 TGraph *grtosca = new TGraph(nftot,tosca_cur,tosca_curratio);
 grtosca->SetName("grtosca");
 TLegend *leg = new TLegend(.25,.25,.5,.5);
TCanvas *can2 = new TCanvas("can2","Q1 (T/m)/Current versus Current",800,800);
 can2->Divide(1,1);
 can2->cd(1);
 grposup->Draw("AP");
 grtosca->Draw("P same");
  grtosca->SetMarkerStyle(22);
  grtosca->SetMarkerColor(kRed);
 grtosca->SetMarkerSize(1.5);
 leg->AddEntry("grposup","Data","p");
 leg->AddEntry("grtosca","Tosca","p");
 leg->Draw();
}
//
void Set_probe(vector<double> &x ,vector<double> &y ) {
 gROOT->Reset();
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.12);
 //
 probe_not_set = kFALSE;
 Int_t vecsize=x.size();
 const Int_t npar=1000;
 Double_t xvec[npar];
 Double_t yvec[npar];
 Double_t resid[npar];
 for (UInt_t i=0; i<x.size(); i++) {
   xvec[i]=x[i]/10.;
   yvec[i]=y[i]/10000.;
 }
 //
 Double_t fit_val=0.;
TCanvas *cprobe = new TCanvas("cprobe","Lakeshore probe B versus B_corr",800,800);
 cprobe->Divide(1,2);
 cprobe->cd(1);
 TGraph *grneg = new TGraph(vecsize,xvec,yvec);
 TGraph *grpos = new TGraph(vecsize,xvec,yvec);
 grneg->Draw("AP");
 grneg->SetMarkerStyle(21);
 grneg->SetTitle("Lakeshore probe B correction versus B");
 grneg->GetXaxis()->SetTitle("Field (T)");
 grneg->GetYaxis()->SetTitle("Field correction (T)");
 grneg->SetMarkerSize(.5);
 grneg->SetMaximum(.03);
 grneg->SetMinimum(-.03);
 grneg->Fit("pol8","Q","",-3.0,0.0);
 neg_field_corr = grneg->GetFunction("pol8");
 grpos->Fit("pol6","Q","",0.0,3.0);
 pos_field_corr = grpos->GetFunction("pol6");
 for (int i=0;i<vecsize;i++) {
   if (x[i] <0) fit_val=neg_field_corr->Eval(xvec[i],0.,0.);
   if (x[i] >0) fit_val=pos_field_corr->Eval(xvec[i],0.,0.);
   //   cout << fit_val << " " << yvec[i] << endl;
      resid[i]=(fit_val-yvec[i])/yvec[i];
 }
 cprobe->cd(2);
 TGraph *gr_resid = new TGraph(vecsize,xvec,resid);
 gr_resid->Draw("AP");
 gr_resid->SetTitle("Field Corr FIT Residual versus Current");
 gr_resid->GetXaxis()->SetTitle("Field (T)");
 gr_resid->GetYaxis()->SetTitle("(FIT-Meas)/Meas");
 gr_resid->SetMarkerStyle(21);
 gr_resid->SetMarkerSize(.5);
 gr_resid->SetMaximum(.1);
 gr_resid->SetMinimum(-.1);
 //
}
