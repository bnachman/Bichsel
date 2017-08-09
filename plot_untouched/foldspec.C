//
// Root Macro to plot Energy Loss spectra from foldtest.for for a set of thin silicon thicknesses 
// and zoom in to the threshold leading edge 
//

//
// Utility to read in spectrum data from one file 
//
int ReadSpec(string* Fname, float Eval[], float Hval[], float Hint[]) {   

  FILE* infile = fopen(Fname->c_str(), "rb");
  if (!infile) {
    printf ("Cannot open %s'\n", Fname->c_str() );
    return 1;
  } else {
    printf ("Input file = %s'\n", Fname->c_str() );
  }
  int ind = 0;
  int n=6;
  int L;
  float E, adc, H, asp, ass;  
  while(n==6) {
    n = fscanf(infile, "%i %f %f %f %f %f", &L, &E, &adc, &H, &asp, &ass);
    if(n==6) {
      Eval[ind] = E;
      Hval[ind] = H;
      Hint[ind] = asp;
      ind++;
    }
  } 
  printf("Read %i lines from %s\n",ind,Fname->c_str());
  fclose(infile);  
  return ind;
}

//
//  Main plotting routine 
//
void foldspec() { 
  //
  // Plot inititializations 
  //
  gStyle->SetMarkerStyle(1);
  gStyle->SetMarkerColor(kCyan);
  gStyle->SetTitleSize(0.06,"xyz");
  gStyle->SetLabelSize(0.04,"xyz");
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadLeftMargin(0.1);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.05);
  gStyle->SetPadGridY(0);
  gStyle->SetPadGridX(0);

  const int nFiles = 6; 
  const std::string specFile[nFiles] = {
   "results/foldspec_muon_1000mev_10micron.dat",
   "results/foldspec_muon_1000mev_15micron.dat",
   "results/foldspec_muon_1000mev_20micron.dat",
   "results/foldspec_muon_1000mev_25micron.dat",
   "results/foldspec_muon_1000mev_30micron.dat",
   "results/foldspec_muon_1000mev_35micron.dat"
  };

  TCanvas *c2 = new TCanvas("c2","Bichsel E-loss Spectra",100,10,700,900);
  c2->Divide(1,2);
  
  TGraph* gr[nFiles]; 
  //int symbol[nFiles] = { 20, 21, 22, 23 };  
  int color[nFiles]  = {  1,  2,  3,  4,  6,  7 };  
  int thmu[nFiles]   = { 10, 15, 20, 25, 30, 35 };  
  //
  //  Loop through the files
  //
  float EeV[200];        //  Input E-loss bin value (eV) 
  float Hdat[200];       //  Input E-loss spectrum probability (/eV) 
  float Hintegral[200];  //  Input E-loss spectrum probability integral(/eV) 
  float Electron[200];   //  E-loss in number of electrons  
  float Hscale[200];     //  E-loss spectrum probability  

  const float Emax  = 3500.;
  const float eHole = 3.70; //  <Eloss> in eV per e/hole pair at 0 degree C

  c2->cd(1);
  
  TLegend* leg =  new TLegend(0.65,0.5,0.85,0.87);

  for(int i=0; i<nFiles; i++) {
    string Ftemp = specFile[i];
    int nLine = ReadSpec(&Ftemp,EeV,Hdat,Hintegral);
    for(int j=0; j<nLine; j++) {
      Electron[j] = EeV[j]/eHole;        // No. of electrons 
      Hscale[j]   = eHole*1.0E3*Hdat[j]; // Prob(E) per 1000 electrons  
    } 
    int nPlot = 0 ;
    for(int j=0; j<nLine; j++) {
      if(Electron[j]<Emax) nPlot++; 
    }   
    gr[i] = new TGraph(nPlot,Electron,Hscale); 
    gr[i]->SetMarkerStyle(1);
    gr[i]->SetMarkerColor(color[i]);
    gr[i]->SetLineColor(color[i]);
    if(i==0) { 
      gr[i]->SetTitle("Energy Loss spectrum");
      gr[i]->GetXaxis()->SetTitle("E Loss (electrons)");
      gr[i]->GetYaxis()->SetTitle("Prob(Ne) / ke");
      gr[i]->GetYaxis()->SetTitleOffset(0.7);
      gr[i]->SetMinimum(0);
      gr[i]->Draw("AL"); 
    } else { 
      gr[i]->Draw("LSAME");
    }
    char thickness[20];
    sprintf(thickness,"%2i micron\n",thmu[i]); 
    leg->AddEntry(gr[i],thickness,"L");
  }
  leg->Draw();
  TText* t1=new TText(1100.,1.6,"1 GeV muon");
  t1->Draw();
  
  //
  // Integrated probability
  //
  c2->cd(2);
  TGraph* gr2[nFiles];
  gStyle->SetPadTopMargin(0.01);


  const float Ecut = 3000.; // Range to plot for Eloss threshold  
  TLegend* leg2 =  new TLegend(0.65,0.5,0.85,0.87);
  const float ymax = 20.;
  float Thresh[200];
  float Ineff[200];
  
  for(int i=0; i<nFiles; i++) {
    string Ftemp = specFile[i]; 
    int nLine = ReadSpec(&Ftemp,EeV,Hdat,Hintegral);
    for(int j=0; j<nLine; j++) {
      Thresh[j] = EeV[j]/eHole;        // No. of electrons 
      Ineff[j]  = 100.*Hintegral[j];   // Inefficiency (%) 
    } 
    int nPlot = 0; 
    for(int j=0; j<nLine; j++) {
      if(Thresh[j]<Ecut) nPlot++; 
    }   
    gr2[i] = new TGraph(nPlot,Thresh,Ineff); 
    gr2[i]->SetMarkerStyle(1);
    gr2[i]->SetMarkerColor(color[i]);
    gr2[i]->SetLineColor(color[i]);
    if(i==0) { 
      gr2[i]->SetTitle("");
      gr2[i]->GetXaxis()->SetTitle("Threshold (electrons)");
      gr2[i]->GetYaxis()->SetTitle("Inefficiency (%)");
      gr2[i]->GetYaxis()->SetTitleOffset(0.7);
      gr2[i]->SetMinimum(0);
      gr2[i]->SetMaximum(ymax);
      gr2[i]->Draw("AL"); 
    } else { 
      gr2[i]->Draw("LSAME");
    }
  }
  leg->Draw();
} 
