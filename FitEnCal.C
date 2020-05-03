
#include <string.h>
#include "TGraph.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TF1.h"

double fitf(double *x, double *par) {
    double fitval = par[0] + par[1]*x[0];
    return fitval;
}

void FitEnCal(const char* fname) {

    char title[200];
    strcpy(title, fname);

    char type[200] = "Energy Calibration ";

    // crea l'oggetto GRAFICO CON ERRORI leggendo da file di testo
    TGraph *g = new TGraph(fname);
    g->SetName("EnCal");
    // formattazione del grafico con comandi testuali
    g->SetMarkerColor(4);
    g->SetMarkerStyle(20);
    g->SetTitle(strcat(type, title));
    g->GetXaxis()->SetTitle("Integrated charge [a.u.]");
    g->GetXaxis()->SetTitleSize(0.04);
    g->GetXaxis()->CenterTitle();
    g->GetYaxis()->SetTitle("Alpha energy [keV]");
    g->GetYaxis()->SetTitleSize(0.04);
    g->GetYaxis()->CenterTitle();

    TCanvas *c1 = new TCanvas("c1"); // apre la finestra grafica

    gStyle->SetStatY(0.35); // y-pos of statpad

    g->Draw("ap"); // disegna il grafico a punti (p) in un nuovo frame (a)

    g->GetYaxis()->SetRangeUser(5000, 6000); // setta il range dell'asse Y

    TF1 *f1 = new TF1("f1", fitf, 0., 15000., 2); // formula definita da funzione

    // setta i valori iniziali dei parametri
    f1->SetParameters(0, 1);
    f1->SetParNames("q", "m");

    TFitResultPtr r = g->Fit(f1, "RS"); // R = Region   S = StoreResults

    gStyle->SetOptFit(1111); // cliccare sulla pad dei risultati con l'editor e provare a cambiare qualcosa...

}


