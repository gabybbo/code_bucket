
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"


void Plot_Range_P(const char* fname) {
    TGraph *g = new TGraph(fname);
    g->SetName("Range");
    g->SetMarkerColor(4);
    g->SetMarkerStyle(20);
    g->SetTitle("Plot Range vs P^{-1}");
    g->GetXaxis()->SetTitle("P^{-1} [mbar^{-1}]");
    g->GetXaxis()->SetTitleSize(0.04);
    g->GetXaxis()->CenterTitle();
    g->GetYaxis()->SetTitle("Range [cm]");
    g->GetYaxis()->SetTitleSize(0.04);
    g->GetYaxis()->CenterTitle();

    TCanvas *c1 = new TCanvas("c1"); // apre la finestra grafica

    g->Draw("ap");
}
