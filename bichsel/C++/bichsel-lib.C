#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <set>
#include <map>
#include <sys/time.h>
#include <math.h>
#include <stdlib.h>
#include <ctime>
//const int PROB=133; // 40MHz/75KHz times 2 pixels
const int MXWAIT=63;
const int NBINS=100;
const int NCOL=112;

//****************************************************************************************
// Global Variables: were mydeclared by Bichsel in a fortran common block
// Bichsel used implicit typing: names starting with i,j,k,l,m,n are integer; others real. 
// NOte that all array indices run fomr 1 to max and are dimensioned as [max+1]. 0 index not used. 
//****************************************************************************************
//variables from common /barray/ 
float f[1405+1],h[1405+1],E[1705+1],DI[1705+1],dE[1705+1],xn;
//variables from common /ut/   
int N1,NU;
float F0,H0,U,um,EX,CZ0,CN,CM0,CM1,CM2,CM3,CM4;
//variables from common /IND/  
int N2,N2P,MIE,MIF,MIH,LEF,LEH,nume,lemx;
float D1;
//variables from common /ABS/  
float atnu,dEdx,exth,thi,xi,rkap;
const float ZA=14., AW=28.086, rho=2.329; // silicon atom
const float pi=3.14159265359, Ry=13.6058; // pi and Rydberg
//variables from common /NML/  
float saxk,Etop,bemx,FSG,zi,su[8+1],nels[2+1];
//variables from common /MEAN/ 
float CMA,CMB,CMD,D2,D3,D4,tdedx,tDD[2+1][4+1],xkmn[200+1];
//variables from common /SPTT/ 
float sig[6+1][1252+1],stp[5+1],tsig[5+1],rM2[5+1],rim[1252+1]; 
//variables from common /EP12/ 
float ep[2][1252+1],dfdE[1252+1],BB[2+1];
//variables from common /ENER/ 
float Emin,Efin,Emax,gam,pkE;
//variables from common /const/
int ners;
//variables from common /EVA/  
int npm,nzch;
float PTM,bg,betasq=1.;
//variables from common/BDUMP/ 
int LEVDMP;   

//Other constants not in common block- I collected them here
const float  elm =  0.511004; //electron masss

float mydec = 1;

//*********************************************************************************************
//The following functions are fomr cov-lib.for
//*********************************************************************************************

void EPRED() {
// Function to read in heps.tab containing the tabulated result of 
// dielectric constant (epsilon). This is used for the cross section
// of small momentum transfer excitations.    
// values read from file go into global arrays ep, rim, and dfdE. etbl is thrown away. 

std::ifstream unit14;
int n2t=0,numt=0,jt=0;
float etbl;
std::string line;

unit14.open("bichdat/heps.tab");

std::getline(unit14, line);
std::istringstream iss(line);
iss >> n2t >> numt;
std::cout << "EPRED, n2t,numt= " << n2t << " " << numt << std::endl;
if (nume != numt) std::cout << " CAUTION: nume & numt differ " << std::endl;
if (N2 != n2t) std::cout << " CAUTION: N2 & n2t differ " << std::endl;

for (int j=1; j<=numt; j++) {
    std::getline(unit14, line);
    std::istringstream iss(line);
    iss >> jt >> etbl >> ep[1][j] >> ep[2][j] >> rim[j];
    dfdE[j] = rim[j] * 0.0092456 * E[j];
//    if (j<10) std::cout << jt << " " << etbl << " " <<  ep[1][j] << " " << ep[2][j] << " " << rim[j] << std::endl;
}
unit14.close();
return;
}

//*********************************************************************************************
void AERED() {
// Function to read in MACOM.TAB containing tabulated result of 
// Generalised Oscillation Strength (GOS) calculations for longitudinal 
// excitation (K+L shells) with large momentum tansfers.  
// values read from file go into global arrays sig. etbl is thrown away. 

std::ifstream unit15;
int n2t=0,numt=0,jt=0;
float etbl;
std::string line;

unit15.open("bichdat/macom.tab");

std::getline(unit15, line);
std::istringstream iss(line);
iss >> n2t >> numt;
std::cout << "AERED, n2t,numt= " << n2t << " " << numt << std::endl;

for (int j=1; j<=numt; j++) {
    std::getline(unit15, line);
    std::istringstream iss(line);
    iss >> jt >> etbl >> sig[6][j];
    if (j<10) std::cout << jt << " " << etbl << " " <<  sig[6][j] << std::endl;
}
unit15.close();
return;
}

//*********************************************************************************************

void EMRED() {
// Function to read in EMERC.TAB containing tabulated result of 
// Generalised Oscillation Strength (GOS) calculations for longitudinal 
// excitation (M shell) with large momentum tansfers.  
// values read from file go into global arrays sig and xkmn. etbl is thrown away. 

std::ifstream unit16;
int n2t=0,numt=0,jt=0;
float etbl;
std::string line;

unit16.open("bichdat/emerc.tab");

std::cout << "EMRED" << std::endl;
std::getline(unit16, line);
std::getline(unit16, line);
std::getline(unit16, line);
std::getline(unit16, line);

for (int j=1; j<=175; j++) {
    std::getline(unit16, line);
    std::istringstream iss(line);
    iss >> jt >> etbl >> sig[6][j] >> xkmn[j];
    if (j<10) std::cout << jt << " " << etbl << " " <<  sig[6][j] << " " <<  xkmn[j] << std::endl;
}
unit16.close();
return;
}

//*********************************************************************************************

void EVANS() {
// Initialization of kinematic parameters 
// Evans, p 891, Uehling Eq (4a).  Date : 15 June 1984
float etbl, xxx, pmom, ptE, beta, Emx;
etbl=0.0;
int jkm;

// open ouptut file in append mode
std::ofstream outfile;
outfile.open("COV.OPA", std::ios_base::app); 
 
std::cout << "Evans: particle mass (PTM) = " << PTM << " MeV" << std::endl;
std::cout << "Select input option [1=E_kin (MeV), 2=P_mom(MeV/c), 3=bet*gam]: ";
std::cin >> jkm;
std::cout << "Now enter the value: ";
std::cin >> xxx;

//Calculate the following differently depending on input option
// * gamma factor (gam)
// * particle kinetic energy (pkE)
// * bata.gamma (bg)
// * particle momentum (pmom)
// * particle total energy (ptE)

if (jkm==1) {
// xxx was kinetic energy
     pkE = xxx;
     gam = xxx/PTM + 1.0;
     bg = sqrt(gam*gam - 1.);
     pmom = PTM*bg;
     }

else if (jkm==2) {
// xxx was momentum
     pmom = xxx;
     bg = xxx / PTM;
     gam = sqrt(bg*bg + 1.0);
     pkE = PTM * (gam - 1.0);
     }

else {
// xxx was beta.gamma
     bg = xxx;                
     gam = sqrt(bg*bg + 1.0);
     pkE = PTM * (gam - 1.0);
     pmom = PTM*bg;
   }

// remaining parameters
   betasq = bg*bg / (1 + bg*bg);
   beta = bg / gam;
   ptE  = PTM * gam;

// Maximum energy transfer   Emax  (MeV) 
// Uehling, also Sternheimer & Peierls Eq.(53)
    Emax = PTM * (gam*gam - 1) / (PTM/(2.*elm) + (2.*elm)/PTM + gam);

// cannot distinguish i/o electrons
    if (npm==4) Emax = pkE / 2.;
        
    std::cout << "particle type = " << npm << "   Emax = " << Emax << " MeV" << std::endl;
    Emx  = (2.*elm) * bg*bg;

// write to output file. Need to properly pad the numbers as
// f11.4, f13.4, and f15.4
outfile << std::endl << "   beta*gamma=" << bg << "   momentum=" << pmom << "   E kinetic of incident particle=" << pkE;
// write to output file. Need to properly pad the numbers as
// f9.6, f12.5, e12.4, e12.4
outfile << "   beta^2=" << betasq << "   gamma=" << gam << "   Emax=" << Emax << " " << Emx << " MeV" << std::endl;


// Change Emax units to eV
    Emax = 1.e6 * Emax;

outfile.close();
return;
}

//*********************************************************************************************

void PREP() {
// Initialization routine for user parameter input & basic constants 
// Particle types   npm = 1    proton
//                      = 2    pion
//                      = 3    alpha
//                      = 4    electron (poistron) 
//                      = 5    kaon

float PMASS[6] = {0., 938.256, 139.578, 3727.328, 0.511004, 497.034};
float Emk=0.;

// open ouptut file in append mode
std::ofstream outfile;
outfile.open("COV.OPA", std::ios_base::app); 

std::cout << "PREP: Ry = " << Ry << " eV" << std::endl;
std::cout << "Enter particle type (1=P, 2=Pi, 3=Alpha, 4=e, 5=K) :";
std::cin >> npm;
std::cout << "Enter silicon thickness in microns :";
std::cin >> exth;

// convert microns into cm
        exth = exth / 10000;

// set particle mass values (MeV) according to input code number  
        PTM = PMASS[npm];

// set charge of incident particle
        zi = 1.;
        if (npm==3) zi =2.;

// properties of absober material (Silicon) 
//    ZA = atomic number 
//    AW = atomic weight
//   rho = density  (g/cm^3) 
//  atnu = # of atoms/cm^3 
        atnu = 6.0222e23 * rho / AW;

// write to output file
outfile << "PREP: particle mass= " << PTM << "MeV, charge=" << zi; 

// Initialization kinematic parameters
        EVANS();
        Efin = Emax;
        saxk = 153540. * zi*zi * rho / (betasq*AW);
//  Saxon Eq (3a) for k
        Emk  = saxk * Emax;

// write to output file
outfile << ", Z=" << ZA << ", A=" << AW << ", t=" << exth << "cm, K/Z=" << saxk << ", k*Emax/Z=" << Emk << std::endl;

outfile.close();
return;
}

//*********************************************************************************************

void PREPE() {
float etbl=0.0, u, ken, exs;

// open ouptut file in append mode
std::ofstream outfile;
outfile.open("COV.OPA", std::ios_base::app); 

// Definitions of energy scale (log) bin size 
// N2 = number of bins for each factor of 2 in energy 
        N2   = 64.;
        nume = 650;
        if (N2 == 64) nume = 1250;
        u    = log(2.) / N2; //check if log is log or ln
        um   = exp(u);
        ken  = log(1839. / 1.5) / u;
        Emin = 1839. / pow(2.,(ken/N2)); 
        //Emin = 1839/2; 
        E[1] = Emin; 
        exs  = 1.;

std::cout <<" PREPE " << N2 << ", " << ken << ", " << E[1] << ", " << Emin << std::endl;

// write to output file. Need to properly pad the numbers as
// I4, f8.3, f9.6, f10.6
outfile << std::endl << " PREPE:   N2=" << N2 << "      Emin=" << Emin << "   u=" << u << " e^u=" << um << std::endl;

        lemx = nume + 450;
        for (int L=1; L<=lemx; L++) {
          exs = exs * um;
          E[L+1] = E[L] * um;
          if ( (L/50)*50 == L) std::cout << " L,E=" << L << ", " << E[L] << ", " << exs << ", " << um << std::endl;
          DI[L]  = -log(1.0 - 1.0/exs) / u;
          dE[L]  = E[L+1] - E[L];
          if (L <= nume) h[L]   = 0.;
          if (E[L] <= Emax) LEH = L;
        }

        if (LEH > nume) LEH = nume;
        Etop = E[nume] * sqrt(um);
        if (Efin > Etop) Efin = Etop;

// write to output file. Need to properly pad the numbers as
// I4, f8.3, f9.6, f10.6
outfile <<  " PREPE: Efin, Etop, Emax=" << Efin << ", " << Etop << ", " << Emax << std::endl;
std::cout <<  " PREPE: Efin, Etop, Emax=" << Efin << ", " << Etop << ", " << Emax << std::endl;

outfile.close();
return;
}

//*********************************************************************************************

void HPART(float bbb, float fft, float sbb, float STPW, float SECM, std::ofstream& outfile) {
float rst, tdedx, secmv, rM2p, del2, TE;

        if (npm < 4) rst = bbb * log (Emax/Efin) + sbb * (1./Efin - 1./Emax) - betasq*(1. - Efin/Emax);

        if (npm >= 4) {
// see FSR-143 and Uehling Eq 9
          TE = pkE * 1.e6;

          std::cout << " ele " << Efin << " " << Emax << " " << TE << " " << gam << std::endl;

          rst= log(Emax/Efin) + log(Emax) - log(TE-Efin) - 1/(1.-Efin/TE);
          rst+= 2 + ((gam-1)*(gam-1) / (gam*gam)) * (1./8. - .5*(Efin*Efin/(TE*TE)));
          rst+=  ((2*gam - 1) / (gam*gam))  *  (log(Emax) - log(TE-Efin));

//        rst=alog(Emax/Efin)+log(Emax)-alog(TE-Efin)  - 1/(1.-Efin/TE)
//     1       + 2 +     ((gam-1)/gam)**2      *  (1./8.-.5*(Efin/TE)**2)
//     2       + ((2*gam - 1) / gam**2) * (alog(Emax) - alog(TE-Efin))

        }

std::cout << " residual dE/dx=" << rst << " " << rst*fft << " MeV/cm" << std::endl;
 
        tdedx = STPW/1.e6 + rst*fft;

outfile <<  " dE/dx=" << tdedx << " MeV/cm " << tdedx/2.329 << " MeV cm^2/g" << std::endl;
std::cout <<  " dE/dx=" << tdedx << " " << tdedx/rho << std::endl;

        secmv = SECM / 1.e6;
        rM2p  = Efin - 0.5 * betasq * Efin*Efin / Emax;
        del2  = secmv - fft * rM2p;

std::cout <<  " M2=" << secmv << "  M2''=" << fft*rM2p << " M2-M2''=" << del2 << std::endl;
outfile <<  " M2=" << secmv << "  M2''=" << fft*rM2p << " M2-M2''=" << del2 << std::endl;

return;
}

//*********************************************************************************************

void SPTS(std::ofstream& outfile) {

float SGM  = 0., STPW = 0., SECM = 0, sbb=720., bbb, fft, he2, eps, rm0;
int j, nlast, jpr = 20, ja  = 20;

        for (j=1; j<=nume; j++) { 
          if (E[j] > Emax) break;
          nlast= j;
          he2  = sig[5][j] * mydec;
          h[j] = he2 / (E[j]*E[j]);
          SGM  = SGM + h[j]*dE[j];
          STPW = STPW + h[j] * E[j] * dE[j];
          SECM = SECM + he2 * dE[j];
          eps  = STPW / SGM;
        }


// format (/1X,'SPTS: nlast  ',i5,' total cross section=',E15.5,3X,'dE/dx=',E15.5,', M2=',E15.5/2x,'see FSR-99 and CCS-9'/)
outfile << std::endl << " SPTS: nlast   " << nlast << " total cross section=" << SGM << "   dE/dx=" << STPW << ", M2=" << SECM << std::endl;
outfile << "  see FSR-99 and CCS-9" << std::endl;

outfile << " final E=" << Efin << ", Emax=" << Emax << ", he2=" << he2 << std::endl;

        bbb = 1. - sbb*betasq / Emax;
        fft = 14. * mydec / 1e6;

outfile << " sbb" << sbb << " eV,   bbb=" << bbb << "   fft=" << fft << std::endl;

        rm0 = bbb * ((1/Efin - 1/Emax) + 2 * (1/Efin*Efin - 1/Emax*Emax)) - betasq * log(Emax/Efin) / Emax;

outfile << " residual M0=" << rm0 << " " << rm0*fft << "/cm" << std::endl;
outfile << "  If residual M0 is large, look for error" << std::endl;

HPART(bbb,fft,sbb,STPW,SECM,outfile);

return;
}

//*********************************************************************************************

void SPECT() {
float fac, blg;
int L, jpr, jpd; 
float S0=0., S1=0., avI=0., avI1=0., pf, tmcb, uef, Q1, Qmin, epbe, sgg, thet, sgh, rmf;

// generate collision spectrum from ep-1,2 and ae

// open ouptut file in append mode
std::ofstream outfile;
outfile.open("COV.OPA", std::ios_base::app); 

     fac = 8. * pi * Ry*Ry * (0.529177e-8)*(0.529177e-8) / (elm * betasq);
     mydec = zi*zi * atnu * fac;
     blg = log ((2.*elm) * bg*bg) - betasq; 

//format f12.10, e12.4, f9.4
outfile << std::endl << "     SPECT F.307:  beta^2=" << betasq << "     atoms per cm^3=" << atnu << "   blg=" << blg << std::endl;

        jpr = 5;
        jpd = 5;
        for (L=1; L<=5; L++) { 
          rM2[L]  = 0;
          tsig[L] = 0;
          stp[L]  = 0;
        }
        bemx = betasq / Emax;
        pf = pkE * 1.e6;
        tmcb = 2. * elm * betasq;
 
// do loop for Fano Eq 47
        for (int j=1; j<=nume; j++) { 
          if (E[j] > Emax) break; //terminate do loop
          if (npm == 4) uef = 1 + E[j]*E[j]/(pf-E[j])*(pf-E[j]) + ((gam-1)*(gam-1)/(gam*gam)) * (E[j]*E[j]/(pf*pf)) - (2*gam-1)*E[j]/(gam*gam*(pf-E[j]));
          if (npm < 4) uef = 1 - E[j] * bemx;     
          if (j == 0) std::cout << " uef=" << uef << std::endl;
          S0  = S0   + dfdE[j] * dE[j];
          avI  += dfdE[j] * log(E[j]) * dE[j];
          avI1 += dfdE[j] * E[j] * log(E[j]) * dE[j];
          S1   += dfdE[j] * E[j] * dE[j];
          Q1 = Ry;
          if (E[j] < 100.)  Q1 = 0.025*0.025 * Ry;
          if (E[j] < 11.9)  Q1 = xkmn[j]*xkmn[j] * Ry;
          Qmin = E[j]*E[j] / tmcb;
          sig[1][j] = 0;
          if ( (E[j] >= 11.9) || (Q1>Qmin) ) sig[1][j] = E[j] * dfdE[j] * log(Q1 / Qmin); 
          epbe = 1 - betasq * ep[1][j];

//  Fano Eq 47
          if (epbe == 0) epbe = 1e-20;
          sgg = E[j] * dfdE[j]*(-.5)*log(epbe*epbe+(betasq*ep[2][j])*(betasq*ep[2][j]));
          thet = atan(ep[2][j] * betasq / epbe);
          if (thet < 0.) thet = thet + pi;         
//  plausible-otherwise I'd have a jump
//  Fano says [p 21]: 'arctan approaches pi for betasq*eps1 > 1'
          sgh = E[j]*E[j] * (betasq-ep[1][j]/(ep[1][j]*ep[1][j]+ep[2][j]*ep[2][j]))*thet;
          sgh = 0.0092456 * sgh;
          sig[3][j] = sgg + sgh;
          sig[4][j] = 2. * sig[6][j] * uef;


          sig[2][j] = 0; sig[5][j] = 0;
          for (L=1;L<=4;L++) {
                tsig[L]  = tsig[L]  + sig[L][j] * dE[j] / (E[j]*E[j]);
                stp[L]   = stp[L]   + sig[L][j] * dE[j] / (E[j]*E[j]);
                rM2[L]   = rM2[L]   + sig[L][j] *  dE[j];
                sig[5][j] = sig[5][j] + sig[L][j];
          }
          tsig[5] = tsig[5] + sig[5][j] * dE[j] / (E[j]*E[j]);
          stp[5]  = stp[5]  + sig[5][j] * dE[j] / E[j];  
          rM2[5]  = rM2[5]  + sig[5][j] * dE[j];

     } //close do loop for Fano Eq 47

std::cout << "  uef= " << uef << std::endl;
// formats 5F12.4/2(28x,5f12.3/) 
outfile << "         Integ. over sig =" << tsig << std::endl; 
outfile << "                            ";
for (L=1;L<=5;L++) outfile << stp[L] << " ";   
outfile << std::endl << "                            ";
for (L=1;L<=5;L++) outfile << rM2[L] << " ";
outfile << std::endl;

//formats (/9x,' S(0)=',f9.5,3x,'ln(I)=',f10.5,3x,S(1)=',f10.3,3x,'L(1)=',f10.3/)
outfile << "          S(0)=" << S0 << "   ln(I)=" << avI << "   S(1)=" << S1 << "   L(1)=" << avI1 << std::endl;

outfile << "  following data without density effect" << std::endl;
outfile << "  S(0)*blg=" << S0*blg << "   2*L(0)=" << 2*avI << std::endl;
outfile << "  S(1)*blg=" << S1*blg << "   2*L(1)=" << 2*avI1 << std::endl;

std::cout << " S(0)=" << S0 << "  L(0)=" << avI << std::endl;

        FSG  = tsig[5] * mydec;
        dEdx = stp[5] * (mydec/1.E6);
        rmf  = rM2[5] * (mydec/1.E6);

// format (/,10X,'Zeff=',F7.3,4X,'# coll/cm=',f11.3,4x,'dE/dx=',F9.4,' MeV/cm',3x,'M2=',f12.4,' keV**2/cm')
outfile << std::endl << "          Zeff=" << S0 << "    # coll/cm=" << FSG << "    dE/dx=" << dEdx << " MeV/cm";
outfile << "   M2=" << rmf << " keV^2/cm" << std::endl;
outfile << " mydec=" << mydec << "  # atoms/cm^3="<< atnu << "  fac=" << fac << std::endl;

// call SPTS routine
SPTS(outfile);

outfile.close();
return;
}

//*********************************************************************************************
//The following functions are fomr cov3.for
//*********************************************************************************************

void SHRINK(){
 
float S = 0.;
int N,L,lla;

        N = MIH - MIE;
        for (lla=1; lla<=LEH; lla++) {
          S += h[lla]*dE[lla+N];
          if (S > EX) break;
        }
        MIH += lla - 1;
        S  = 0.;

        for (L=LEH; L>0; L--) { //this took 4 lines in Fortran code
          S += h[L] * dE[L+N];
          if (S > EX) break;
        }

        LEH = L - (lla-1);

        for (L=1; L<=LEH; L++) h[L] = h[L+lla-1];
        for (L=(LEH+1); L<=nume; L++) h[L]=0;

  return;
}

//*********************************************************************************************

void NORMAL(std::ofstream& outfile1) {

float X=1.,Y,Z,S,T=1.,cmq=0.,EC;
int N,L,LE;
 
        Y = CM1 * 2.0;
        Z = CM2 * 2.0;

        outfile1 << " NORM "<<MIE<<", "<<MIH<<"   H0= "<<H0<<"   xn,Y,Z= "<<xn<<", "<<Y<<", "<<Z<<std::endl;
//     format(' NORM',2i5,'  H0=',1pe12.5,'  xn,y,z=',0p2f12.3,1pe12.5)

        CM0 = H0;
        CM1 = 0.;
        CMA=0;
        N = MIH - MIE;

        for (L=1; L<=LEH; L++) { 
          LE  = L + N;
          S   = h[L] * dE[LE];
          CM0 = CM0 + S;
          CM1 = CM1 + S*E[LE];
          CMA += S*E[LE]*E[LE];
        }

        if (CM0 != 0) cmq = CM1 / CM0;

        outfile1 << "       area= "<<CM0<<"   straight mean= "<<CM1<<"   CM1/CM0= "<<cmq<<std::endl;
//     format  (7x,'area=',1pe12.5,3x,'straight mean=',1pe12.5,3X,'CM1/CM0=',0pf14.4)

        if ((CM0 - H0) != 0.) T = (1. - H0) / (CM0 - H0);
        CM1 *= T;
        CM2 = 0.;
        CM3 = 0.;
        CM4 = 0.;
        for (L=1; L<=LEH; L++) h[L] *= T;

        std::cout<<" N,LEH= "<<N<<", "<<LEH<<std::endl;
        outfile1<<" N,LEH= "<<N<<", "<<LEH<<std::endl;

        for (L=1; L<=LEH; L++) {
           LE = L + N;
           EC = E[LE] - CM1;
           S  = h[L] * dE[LE];
           CM2 = CM2 + S*EC*EC;
           CM3 = CM3 + S*EC*EC*EC;
           CM4 = CM4 + S*EC*EC*EC*EC;
        }

        xn  *= CM0;
        if (Y != 0.) Y = CM1 / Y;
        if (Z != 0.) Z = CM2 / Z; 

        outfile1 << "       Precision control, "<<xn<<", "<<CM0;
        outfile1 <<", "<<CM1<<", "<<CM2<<", "<<CM3<<", "<<CM4<<std::endl;
        outfile1 <<"                        mean= "<<Y<<"    variance= "<<Z<<std::endl;
//     format  (7X,'Precision control, ',1P6E13.5,
//             /24x,'mean=',e12.5,'   variance=',e12.5)

  return;
}

//********************************************************************************************* 

void ZERO(std::ofstream& outfile1){
 
float xs = 0.;
int N,L;

        N = MIH - MIE;
        outfile1 << " ZERO "<<MIH<<", "<<MIE<<", "<<N<<", "<<LEH<<", "<<F0<<std::endl;

        for (L=1; L<=LEH; L++) xs = xs + h[L]*dE[L+N];
        xs = (1.-F0)*(1.-F0) / xs; 
        outfile1 <<(L+N)<<", "<<xs<<", "<<h[L]<<", "<<dE[L+N]<<std::endl;

        for (L=1; L<=LEH; L++) h[L] *= xs;
        N   = MIH - MIF;
        MIH = MIF;
        LEH = LEH + N;
        if (LEH > nume) LEH = nume;

        for (L=LEH; L>0; L--) h[L+N] = h[L]; //this took 5 lines in fortran looping up and calculating index
        for (L=1;L<=N;L++) h[L]=0;

        std::cout<<" ZERO: "<<N<<"   "<<LEF<<", "<<LEH<<", "<<(LEH+1)<<"   "<<MIE<<", "<<MIF<<", "<<MIH<<std::endl;
//         format  (' zero:',i3,2(3x,3i5))

        for (L=1;L<=LEF;L++) h[L] += 2.*F0*f[L];         
//   f[L] is the h of the previous convolution
//   This seems to be correct (18 July 1984)

  return;
}

//*********************************************************************************************

void RESET(std::ofstream& outfile1){
//  called from FOLD
//  LEF = LEH initially

int L,J;
float S;
 
        if (N2 < 128) {
           N2 = N2 * 2;
           U = log(2.) / (float)N2;
           outfile1 << " RESET: "<<MIE<<", "<<MIH;
           outfile1 << "  Coordinate change: doubling of point grid "<<N2;
           outfile1 << " points on a factor of 2"<<std::endl;
//        format  (1x,'RESET: ',2i4,'  Coordinate change: doubling of ',
//                 'point grid',i4,' points on a factor of 2',/)
//                  write (3,*) ' LEH,MIH,MIE=',LEH,MIH,MIE

           for (L=LEF; L>0; L--) f[2*L] = f[L];  // this is from top down. Took 5 lines in fortran.
                        
           for (L=4; L<=(2*LEF); L+=2) f[L-1] = (f[L] + f[L-2]) / 2.;

           LEF = 2 * LEF + 1;
           MIF = 2 * MIF;
           MIE = MIF;
           
           for (J=1;J<=LEH;J++) {
             S     = (float)(J+MIE);
             E[J]  = exp(S*U) * Emin;
             dE[J] = E[J]*U;
             S = (float)J;
             DI[J] = -log(1. - exp(-S*U)) / U;
           }
           MIE = MIF;   
        } //end of if N2<128       

             
//   from 4th Line above
        for (J=1;J<=LEH;J++) {
           S     = (float)(J+MIE);
           E[J]  = exp(S*U) * Emin;
           dE[J] = E[J]*U;
        }
        outfile1 << " LEH,MIH,MIE= "<<LEH<<", "<<MIH<<", "<<MIE<<std::endl;
  return;
}

//*********************************************************************************************

void FOLD(std::ofstream& outfile1){
// On 18 July 1984, I have some questions whether ST 62 is correct
        
int k,LE,LF,JH,JF,FLF,LFF;
int LH;
float S;
 
        for (LH=1; LH<=1250; LH++) {
          f[LH] = h[LH];
          h[LH] = 0.;
        }

        F0  = H0;
        LEF = LEH;
        MIF = MIH;
        if (LEF < 80) RESET(outfile1);
        H0  = F0*F0;
        MIH = MIF + N2;
        LEH = LEF;
        if (LEH > 1250) LEH = 1250;

        std::cout<< "  FOLD: MIH,LEF,LEH,nume= "<<MIH<<", "<<LEF<<", "<<LEH<<", "<<nume<<std::endl;
        outfile1<< "  FOLD: MIH,LEF,LEH,nume= "<<MIH<<", "<<LEF<<", "<<LEH<<", "<<nume<<std::endl;

        for (LH=1;LH<=LEH;LH++) {
          JH = LH + MIH;
          for (LF=1;LF<=LH;LF++){
            JF = LF + MIF;
             k = JH - JF;
             FLF = JH - MIF - DI[k] + 1.E-20;
             LFF = FLF;
             LE  = JF - MIE;
             S   = (float)(FLF - LFF);
             if (LFF == 0) LFF = 1;
             h[LH] = h[LH] + f[LF] * ((1.0-S)*f[LFF]+S*f[LFF+1]) * dE[LE];
           }
        h[LH] = h[LH] - f[LH]*f[LH] * 0.5*dE[LE];
        }
        for (LH=1;LH<=LEH;LH++) h[LH] *= 2.0;
        if (F0 > EX) ZERO(outfile1);
        SHRINK();
        NORMAL(outfile1);
  return;
}

//*********************************************************************************************

void OUTPUT(std::ofstream& outfile1, std::ofstream& outfile2) {               
  
        float B, C, D, S1, S2, S3, S4, X;   
        float ASP[1252],ASS[1252];
        float bax, bax1, dmp, dmp1,hhun;
        int hmax, N, kmax, NPP, NPM, nskip, k,J,k2,J2;
        int llow, lup;

        time_t t = time(0);   // get time now
        struct tm * now = localtime( & t );
        std::cout << (now->tm_year + 1900) << '-' << (now->tm_mon + 1) << '-' <<  now->tm_mday << "  " 
        << now->tm_hour << ":" << now->tm_min << std::endl;


        B = CM2 / (CM1*CM1);
        C = CM3 / sqrt(CM2*CM2*CM2);
        D = CM4 / (CM2*CM2);
        S1 = CN * D1;

        outfile1 << " OUT "<<B<<", "<<C<<", "<<D<<", "<<S1<<", "<<CN<<std::endl;

        S2 = CN * D2 / (S1*S1);
        S3 = D3 / sqrt(D2*D2*D2 * CN);
        S4 = D4 / (D2*D2 * CN) + 3.;
 
        outfile1 << "  OUTP: zero component = "<< H0 << std::endl << "                              mean";
        outfile1 << "        variance/mean^2    central3/var^(1.5)";
        outfile1 << "     central4/var^2" << std::endl;
        outfile1 << "   actual values= " << CM1<<", "<<B<<", "<<C<<", "<<D<<std::endl;
        outfile1 << "   theoret values " << S1<<", "<<S2<<", "<<S3<<", "<<S4<<std::endl;
//     format (2x,'OUTP: zero component =',1pe16.4/30x,'mean',
//          9x,'variance/mean**2',4x,'central3/var**1.5',
//          20H     central4/var**2/3x,'actual values=',4e20.4,/     //note old Hollerith 20H instead of quotes
//          3x,'theoret values',4e20.4,/)

        outfile1 << " OUTPUT:" << CM1 << ", " << S1 << "   ratio= " << (CM1/S1) << std::endl;
//     format (1x,'OUTPUT:',1p2e13.6,'  ratio=',0pf10.6)
        std::cout << " OUTPUT:" << CM1 << ", " << S1 << "   ratio= " << (CM1/S1) << std::endl;

        X = 1.;
        N = MIH - MIE;

      if (N1 > NU-ners) {

           outfile2 <<LEH<<", "<<N1<<", "<<N2<<", "<<N2P<<", "<<N;
           outfile2 <<", "<<H0<<", "<<thi<<", "<<xi<<", "<<rkap<<std::endl;

                if (N2 != N2P) {
  
                   outfile2 << lemx << std::endl;
                   for (k=1;k<=lemx;k++) {
                      outfile2 << E[k]; //  format (8f10.3) 
                      if ((k+1)%8) outfile2 << ", ";
                      else outfile2 << std::endl;
                   } 
                   outfile2 << std::endl;
                   N2P = N2;
                } // if N2 != N2P 

        ASP[1] = 0;
        ASS[1] = 0;
        hmax   = 0;


        for (k=1; k<=LEH; k++) {
           if (h[k]>= hmax) { 
             hmax = h[k];
             kmax = k;
           }
           J = k + N;
           ASP[k+1] = ASP[k] + h[k] * dE[J];
           ASS[k+1] = ASS[k] + h[k] * E[J] * dE[J];
         }

         std::cout<<" kmax+N= "<<kmax+N<<"   E= "<<E[kmax+N]<<"   h= "<<h[kmax]<<std::endl;
         outfile1<<" kmax+N= "<<kmax+N<<"   E= "<<E[kmax+N]<<"   h= "<<h[kmax]<<std::endl;
//     format (1x,'lmax+N=',i4,'  E=',f11.2,'  h=',f12.6)


         bax = (h[kmax] - h[kmax-1]) / (h[kmax] - h[kmax+1] );
         dmp = E[kmax+N] + 0.5 * dE[kmax+N] * (bax - 1.) / (bax + 1.);
         bax1 = (h[kmax-1] - h[kmax-2]) / (h[kmax-1] - h[kmax] );
         dmp1 = E[kmax+N-1] + 0.5 * dE[kmax+N-1]*(bax1 - 1.)/(bax1 + 1.);
         outfile1 <<" bax,dmp= "<<bax<<", "<<dmp<<"  lower:"<<bax1<<", "<<dmp1<<std::endl;

  	 hhun = 0.01 * h[kmax];
         for (k=1;k<=LEH;k++) {
           if (k < kmax && h[k] < hhun) llow = k;
           if (k >= kmax && h[k] >= hhun) lup  = k;
         }
         outfile1 <<" lower & upper cutoff: "<<llow<<", "<<E[llow+N]<<", "<<h[llow]<<", "<<lup;
         outfile1 <<", "<<E[lup+N]<<", "<<h[lup]<<std::endl;
         std::cout <<" lower & upper cutoff: "<<llow<<", "<<E[llow+N]<<", "<<h[llow]<<", "<<lup;
         std::cout <<", "<<E[lup+N]<<", "<<h[lup]<<std::endl;
//     format (' lower & upper cutoff:',2(i6,0pf9.1,1pe12.4))

         outfile2 << kmax+N<<", "<<E[kmax+N]<<", "<<h[kmax]<<"  "<<" kmax+N, E, h"<<std::endl;
         outfile2 << LEH<<", "<<dmp<<"  LEH, dmp"<<std::endl;

         for (k=1;k<=LEH;k++) {
            outfile2 << h[k]; // format (1p6e12.5)
            if ((k+1)%6) outfile2 << ", ";
            else outfile2 << std::endl;
         } 
         outfile2 << std::endl;
        
         NPP   = (lup + llow) / 2; 
         nskip = 1;
         NPM = (lup - llow) / 2;
         if (NPM > 70) nskip = 2;
         NPM+= nskip;
         std::cout<<" N,npp,npm,nskip= "<<N<<", "<<NPP<<", "<<NPM<<", "<<nskip;
  
         outfile1 << "       J    E/eV        phi(E)       dE/dx         M2   ";
         outfile1 << "       J    E/eV        phi(E)       dE/dx         M2   "<<std::endl;
//  format ( 2(7x,'j',4x,'E/eV',8x,'phi(E)',7x,'dE/dx',9x,'M2',3x))

         for (k=llow;k<=NPP;k+=nskip) {
            J  = k + N;
            J2 = J + NPM;
            k2 = k + NPM;
            outfile1 <<"   "<<J<<"   "<<E[J]<<"   "<<h[k]<<"   "<<ASP[k]<<"   "<<ASS[k];
            outfile1 <<"   "<<J2<<"   "<<E[J2]<<"   "<<h[k2]<<"   "<<ASP[k2]<<"   "<<ASS[k2];
//        format (3x,i5,0pF10.1,1p3E13.4,3x,i5,0pF10.1,1p3E13.4)
         }

     } // if N1 < NU-ners
  return;
}
//*********************************************************************************************

void CONV() {

int JXT, k;
float xmc,xx,S, STPP, DK, DQ, dmmpl, PB;

// open ouptut files in append mode
   std::ofstream outfile1, outfile2;
   outfile1.open("COV.OPA", std::ios_base::app); 
   outfile2.open("COV.SPE", std::ios_base::app);

   outfile1 << "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *" << std::endl;
   outfile1 << "   CONV   t=" << exth << "   FSG=" << FSG << std::endl;

   xmc = exth * FSG;
//number of collisions in thickness exth
   JXT = (int)(log(xmc)/log(2.0)) + 1;
   std::cout << "  JXT=" << JXT;
   xx = xmc / pow(2.0,JXT);
   CZ0 = xx / 1024.;
   NU  = JXT + 10 + 1;

//  format (3x,'F.386: NU=',2i3,4x,'t=',f10.7,' cm',
//          4X,'# coll=',F11.3,4X,'xx=',F9.5,4X,'CZ0=',1pe12.5)
   outfile1 << "   F.386: NU=" << JXT << ", " << NU << "    t= " << exth << "cm    #coll=" << xmc;
   outfile1 << "    xx=" << xx << "    CZ0=" << CZ0 << std::endl;

   outfile2 << CZ0 <<", "<< Emax <<", "<< betasq <<", "<< pkE <<", "<< bg <<", "<< PTM <<", "<< zi << std::endl;
   outfile2 << lemx <<", "<< NU <<", "<< U <<", "<< um <<", "<< Emin << std::endl;

   for (k=1; k<=lemx; k+=6) outfile2 <<E[k]<<", "<<E[k+1]<<", "<<E[k+2]<<", "<< E[k+3]<<", "<<E[k+4]<<", "<<E[k+5]<<std::endl;
// format (1p6e13.6)

        H0  = 0.;
        CM1 = 1.0;
        CM2 = 1.0;
        xn  = 1.;
        EX  = 1.e-15;              
        MIE = 0;                 
        MIF = 0;
        MIH = 0;
        NORMAL(outfile1);

        xn = 1.;
        D1 = CM1;
        D2 = CM2 + CM1*CM1;
        D3 = CM3 + 3.0*CM2*CM1 + CM1*CM1*CM1;
        D4 = CM4 + 4.0*CM3*CM1 + 6.0*CM1*CM1*CM2 + CM1*CM1*CM1*CM1;
        S  = D2 / D1;

        outfile1 << "  conv F.612:    Initial distribution" << std::endl;
        outfile1 << "  delta 1= " << D1 << "   delta 2=" << S << "   exp(u)= " << um << std::endl;
// format (2x,'conv  F.612:',3X,'Initial distribution',/,
//         '  delta 1=',1PE12.4,'  delta 2 =',E12.4,'  exp(u)=',E12.5)

        STPP = 2. * D1 * saxk * CM0 / ZA;
        std::cout << " Eav " << D1 <<", "<< CM0 <<", "<< STPP <<", "<< CMA << std::endl;
        DK   = ((2. * dEdx / STPP) / ZA) * CMA;
        DQ   = DK / Emax - 1.;
        SHRINK();
 
        H0 = 1. - CZ0;
        for (k=1; k<=LEH; k++) h[k] *= CZ0;
//   first convolution: H0 + H(E)

        CN  = CZ0;
        CM1 = CZ0 * D1;
        CM2 = CZ0 * D2;
        thi = exth / pow(2,JXT+10);
        xi  = saxk * thi * ZA;
        rkap = xi / Emax;
        outfile1 << std::endl << "  F 401: "<<NU<<", "<<LEH<<", "<<Emin<<", "<<D1<<", "<<DQ<<", "<<Emax<<std::endl;
//   format(/2x,'F 401: ',2i5,1p5e14.7)

        for (N1=1; N1<=NU; N1++) {
           std::cout << "convol: N1 " << N1;
           thi  = 2. * thi;
           xi   = 2. * xi;
           rkap = 2. * rkap;
           CN   = 2. * CN;

    outfile1 << "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *" << std::endl;
    outfile1 << " CONVOL NUMBER= "<<N1<<" LEH= "<<LEH<<"  MIE, MIH= "<<MIE<<", "<<MIH;
    outfile1 << "    mean collision number= "<< CN << "   CZ0= "<< CZ0 << std::endl; 
//    format (/,' CONVOL NUMBER=',i3,' leh=',i4,2x,'MIE,MIH=',
//            2i5,4x, 'mean collision number=',1PE12.4,'  CZ0=',f9.6,/)

           N2P = N2;
           FOLD(outfile1);

//  Landau-Vavilov parameters
           PB = -.4227843351 - log(rkap) - betasq;
           dmmpl = xi * (PB + .225);

     outfile1 << " t = "<<thi<<"  xi= "<<xi<<"  kappa= "<<rkap;
     outfile1 << " <lam>= "<<PB<<"   Landau theory <del>-dmp= "<<dmmpl<<std::endl;


//    format (/1x,'t =',1pe12.5,2x,'xi=',e12.5,2x,'kappa=',e12.5,
//            2x,'<lam>=',0pf12.6/3x,'Landau theory <del>-dmp=',f12.3)

// **** This line was not doing anything. N not used in this function. Could it be a typo? Should be N1?
//         N = MIH - MIE;
// *****2017*******************************************************************************************

     outfile1 << " CONV, F.408 "<<N1<<" ,"<<LEH<<" ,"<<CN<<" ,"<<xi<<" ,"<<rkap<<" ,"<<PB<<std::endl;
//    format (' CONV, F.408',2i5,f9.3,1p3e11.4/)
 
           OUTPUT(outfile1, outfile2);
        } // end for N1=1 to NU loop

  outfile1.close();
  outfile2.close(); 
  return;
}

//*********************************************************************************************
// Finally the main program from cov.for
//*********************************************************************************************

int mymain() {

std::ofstream outfile;

outfile.open("COV.OPA");
time_t t = time(0);   // get time now
struct tm * now = localtime( & t );
outfile << (now->tm_year + 1900) << '-' << (now->tm_mon + 1) << '-' <<  now->tm_mday; 
outfile << "  " << now->tm_hour << ":" << now->tm_min << std::endl;
outfile.close();

PREP();  //initializes variables in memory
PREPE(); //initializes variables in memory
EPRED(); //reads in data from file
AERED(); //reads in data from file
EMRED(); //reads in data from file
SPECT(); //does something

outfile.open("COV.SPE");
outfile << ners <<", "<< betasq <<", "<< tdedx <<std::endl;
outfile << (now->tm_year + 1900) << '-' << (now->tm_mon + 1) << '-' <<  now->tm_mday; 
outfile << "  " << now->tm_hour << ":" << now->tm_min << std::endl;
outfile << " COV.SPE , t= " << exth << std::endl;
outfile.close();

CONV(); //does most of the work. Calls OUTPUT, FOLD, RESET, ZERO, NORMAL, SHRINK

return(0);
}
