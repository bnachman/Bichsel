//
// Root Macro to plot spectra for energy Loss per micron ("dE/dX") spectra from foldtest.for 
// for various silicon thicknesses as an illustration why "dE/dX" is a dubious concept (which implies
// scalable behavior at any thickness), as Hans Bichsel always complained. 

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
void dEdX() { 
  //
  // Plot inititializations 
  //
  gStyle->SetMarkerStyle(1);
  gStyle->SetMarkerColor(kCyan);
  gStyle->SetTitleSize(0.06,"xyz");
  gStyle->SetLabelSize(0.04,"xyz");
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.06);
  gStyle->SetPadGridY(0);
  gStyle->SetPadGridX(0);

  const int nFiles = 4; 
  const std::string specFile[nFiles] = {
   "results/foldspec_muon_1000mev_1micron.dat",
   "results/foldspec_muon_1000mev_5micron.dat",
   "results/foldspec_muon_1000mev_20micron.dat",
   "results/foldspec_muon_1000mev_300micron.dat"
  };

  TCanvas *c2 = new TCanvas("c2","Bichsel dE/dx Spectra",100,10,900,600);
  
  TGraph* gr[nFiles]; 
  int color[nFiles]  = {  1,  7,  4,  2 };  
  int thmu[nFiles]   = { 1, 5, 20, 300 };  
  //
  //  Loop through the files
  //
  float EeV[400];        //  Input E-loss bin value (eV) 
  float Hdat[400];       //  Input E-loss spectrum probability (/eV) 
  float Hintegral[400];  //  Input E-loss spectrum probability integral(/eV) 
  float Electron[400];   //  E-loss in number of electrons  
  float Hscale[400];     //  E-loss spectrum probability  

  const float Emax  = 120.;
  const float eHole = 3.70; //  <Eloss> in eV per e/hole pair at 0 degree C

  c2->cd(1);
  
  TLegend* leg =  new TLegend(0.65,0.65,0.85,0.87);

  for(int i=0; i<nFiles; i++) {
    string Ftemp = specFile[i];
    int nLine = ReadSpec(&Ftemp,EeV,Hdat,Hintegral);
    for(int j=0; j<nLine; j++) {
      Electron[j] = EeV[j]/eHole/thmu[i]; // No. of electrons per micron 
      Hscale[j]   = eHole*Hdat[j]*thmu[i];  // Prob(E)   
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
      gr[i]->SetTitle("dE/dx spectrum");
      gr[i]->GetXaxis()->SetTitle("E Loss (electrons)/micron");
      gr[i]->GetYaxis()->SetTitle("Prob(E)");
      gr[i]->GetYaxis()->SetTitleOffset(0.9);
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
  TText* t1=new TText(30.,0.04,"1 GeV muon");
  t1->Draw();
  
} 
