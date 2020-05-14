
// TODO: vmax with and without integral discriminant? yes to show it is integral independent?

#include <Riostream.h>
#include <stdlib.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TRandom.h>
#include "TNtuple.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TMath.h"
#include <vector>

#ifndef LABNUC
#define LABNUC
struct bragg_signal {
  short int s[128];
};
#endif

int AnaBragg(const char *filename, int intto=128, int nsig=0) {

  // APERTURA FILE E LETTURA CONTENUTI

  bragg_signal signal;

  TFile *fin=new TFile(filename);
  if (!fin->IsOpen()) {
    std::cout << "file not found! " << std::endl;
    return -1;
  }

  TTree *tree = (TTree*)fin->Get("bragg");
  if (!tree) {
    std::cout << "Bragg tree not found! " << std::endl;
    return -2;
  }

  TBranch *br = tree->GetBranch("signals");
  if (!br) {
    std::cout << "Signal branch not found! " << std::endl;
    return -3;
  }

  br->SetAddress(&signal);
  int nev = br->GetEntries(); //numero eventi
  std::cout << "Number of events in file : " << nev << std::endl;

  // ANALIZZA EVENTO x EVENTO

  // altri parametri iniziali DA VERIFICARE ED EVENTUALMENTE MODIFICARE
  float thr_frac = 0.3; // soglia rispetto al vmax per il calcolo della larghezza
  int intfrom = 0;// regione di integrazione da 0 a intto
  if (intto>128) intto=128;
  int blfrom = 85, blto = 128; // regione per il calcolo della baseline (verificato)

  // valori da salvare nella ntupla
  float bl; // baseline evento x evento
  float integral; // integrale di carica
  float vmax; // massimo relativo alla bl
  float width; // larghezza temporale dei segnali

  // file di output
  char outfilename[200];
  strcpy(outfilename,"anabragg_");
  const char *cc=strrchr(filename,'/');
  if (cc) {cc++; strcat(outfilename,cc);}
  else strcat(outfilename,filename);

  TFile *fout=new TFile(outfilename,"RECREATE"); // output file

  TNtuple *nt=new TNtuple("nt","","ev:vmax:integral:width:baseline");

  int maxev=nev;
  if (nsig && nsig<nev) maxev=nsig; // controllo sul massimo di eventi

  // PRIMO LOOP PER CALCOLO BASELINE

  // array per salvare valori di bl da scrivere poi su il tree
  float* bl_arr = new float[maxev];

  // variabili per ricavare una stima unica del valore della baseline per questo dataset
  float sum_bl = 0;
  float sum_err_bl = 0;

  for (int j=0; j<maxev; j++){
      br->GetEntry(j); // recupero evento

      //calcolo della baseline
      float bl_sum=0;
      float bl_sum2=0;

      for(int i=blfrom; i<blto; i++){
          bl_sum+=signal.s[i];
          bl_sum2+=signal.s[i]*signal.s[i];

      }

      bl=1.0*bl_sum/(blto-blfrom); // baseline dell'evento j-esimo
      float bl_sigma=TMath::Sqrt((bl_sum2-1.0*bl*bl*(blto-blfrom))/(blto-blfrom-1)); // j-dev std
      float bl_err=bl_sigma/TMath::Sqrt(blto-blfrom);

      bl_arr[j]=bl; //salvo valori di ciascun evento nell' array

      //aggiorno somme per calcolo della bl_avg:
      sum_bl+= bl/(bl_err*bl_err);
      sum_err_bl+= 1/(bl_err*bl_err);

  }

  // calcolo del valore medio (media pesata)

  float bl_avg = sum_bl/sum_err_bl;
  float bl_avg_err = TMath::Sqrt(1/sum_err_bl);

  std::cout << "baseline media = " << bl_avg << "+-" << bl_avg_err << std::endl;

  
  // LOOP SUGLI EVENTI
  for (int i=0; i<maxev; i++) {

    // recupera l'evento
    br->GetEntry(i);

    // inizializza a zero
    bl=bl_arr[i]; //carico il valore calcolato in precedenza per salvarlo nella ntupla
    integral=0;
    vmax=0;				     
    width=0;

    int pos_max=0; // posizione temporale del massimo
    
    // calcolo integrali e vmax
    for (int j=intfrom; j<intto; j++) {
      integral += (signal.s[j] - bl_avg); // somma eventi
      if ( (signal.s[j] - bl_avg) > vmax ){
          vmax = (signal.s[j] - bl_avg); // calcolo vmax
          pos_max = j;
      }
    }
    integral += gRandom->Rndm(); // serve per conversione a float da int, contrario del troncamento

    // calcolo width:

    int t_start=0;
    int t_end=0;
    float v_target=vmax*thr_frac;

    for (int j=intfrom; j<pos_max; j++){
        if((signal.s[j]-bl_avg) <= v_target) {
            t_start=j;
        }
    }

    t_start+=0.5; //media con punto successivo che potrebbe essere più vicino a v_target

    for (int j=pos_max; j<intto; j++){
        if(((signal.s[j]-bl_avg) >= v_target)) {
            t_end=j;
        }
    }

    t_end+=0.5; //media con punto successivo che potrebbe essere più vicino a v_target

    width=t_end-t_start;
    width/=10; // conversione in us

    nt->Fill(i,vmax,integral,width,bl);
  }
  std::cout << maxev << " events analyzed..." << std::endl;

  fout->Write();
  fout->Close();

  fin->Close();

  new TFile(outfilename); // riapre il file dei risultati

  return 0;
}
