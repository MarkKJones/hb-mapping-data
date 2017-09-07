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
void Get_three_col_data_from_file(TString filename,vector<double> &,vector<double> & ,vector<double> & );
void Plot_hb_bi(vector<double> & ,vector<double> & );
void Set_probe(vector<double> & ,vector<double> & );
Double_t field_corr(Double_t field );

// Global parameters
vector<double> nominal_field_corfile_data;
vector<double> corr_field_corfile_data;
vector<double> HB_central_current_data;
vector<double> HB_central_field_data;
TF1 *neg_field_corr;
TF1 *pos_field_corr;
TF1 *cur_field;
//
void run_code()
{
  
  Get_two_col_data_from_file("Lakeshore-probe.dat",nominal_field_corfile_data,corr_field_corfile_data);
  Get_two_col_data_from_file("HB-central-B-I.dat",HB_central_current_data,HB_central_field_data);
  Set_probe(nominal_field_corfile_data,corr_field_corfile_data);
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
 for (int i=0;i<HB_neg_up_current.size();i++) {
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
 for (int i=0;i<HB_pos_up_current.size();i++) {
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
//
void Get_three_col_data_from_file(TString filename,vector<double> &x,vector<double> &y,vector<double> &dy )
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
    if (ncomma ==2) {
        chi=sc.Index(',');
        TString temp(sc(0,chi));
	x.push_back(temp.Atof());
        temp=sc(chi+1,sc.Sizeof());
        chi=temp.Index(',');
        temp=sc(0,chi);
	y.push_back(temp.Atof());
        temp=sc(chi+1,sc.Sizeof());
	dy.push_back(temp.Atof());
	//cout << " data " << x[nline] << " "<< y[nline]<< " " <<  dy[nline]<<endl;
         nline++;
    }
   }
}
//
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
 Double_t fit_val;
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
