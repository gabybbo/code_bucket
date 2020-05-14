#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TLegend.h"

void PlotIntegral(char* file_Pu, char* file_Am, char* file_Cm){

    TGraphErrors* g1 = new TGraphErrors(file_Pu);
    TGraphErrors* g2 = new TGraphErrors(file_Am);
    TGraphErrors* g3 = new TGraphErrors(file_Cm);

    TCanvas* c1 = new TCanvas();
    c1->SetGrid();

    g1->GetXaxis()->SetTitle("Pressure [mbar]");
    g1->GetYaxis()->SetTitle("Integral [a.u.]");
    g1->GetYaxis()->SetRangeUser(4050,4650);

    g1->SetTitle("Plot Integral vs Pressure");
    g1->SetMarkerColor(kBlue);
    g1->SetLineColor(kBlue);
    g1->SetLineWidth(2);
    g2->SetMarkerColor(kRed);
    g2->SetLineColor(kRed);
    g2->SetLineWidth(2);
    g3->SetMarkerColor(kGreen+3);
    g3->SetLineColor(kGreen+3);
    g3->SetLineWidth(2);

    TLegend* l = new TLegend();

    l->AddEntry(g1, "^{239}Pu");
    l->AddEntry(g2, "^{241}Am");
    l->AddEntry(g3, "^{244}Cm");

    g1->Draw("ap");
    g2->Draw("samep");
    g3->Draw("samep");

    l->Draw();


}
