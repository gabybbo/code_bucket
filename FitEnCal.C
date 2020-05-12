
#include <string.h>
#include "TMatrixDSym.h"
#include "TLine.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TF1.h"

double fitf(double *x, double *par) {
    double fitval = par[0] + par[1]*x[0];
    return fitval;
}

void FitEnCal(const char* fname) {

    // crea l'oggetto GRAFICO CON ERRORI leggendo da file di testo
    TGraphErrors *g = new TGraphErrors(fname);
    g->SetName("EnCal");
    // formattazione del grafico con comandi testuali
    g->SetMarkerColor(4);
    g->SetMarkerStyle(20);
    g->SetTitle("Energy Calibration");
    g->GetXaxis()->SetTitle("Alpha energy [keV]");
    g->GetXaxis()->SetTitleSize(0.04);
    //g->GetXaxis()->CenterTitle();
    g->GetYaxis()->SetTitle("Integrated charge [a.u.]");
    g->GetYaxis()->SetTitleSize(0.04);
    //g->GetYaxis()->CenterTitle();

    TCanvas *c1 = new TCanvas("c1"); // apre la finestra grafica

    gStyle->SetStatY(0.35); // y-pos of statpad

    g->Draw("ap"); // disegna il grafico a punti (p) in un nuovo frame (a)

    TF1 *f1 = new TF1("f1", fitf, 0., 15000., 2); // formula definita da funzione

    // setta i valori iniziali dei parametri
    f1->SetParameters(0, 1);
    f1->SetParNames("q", "m");
    f1->SetLineWidth(1);

    TFitResultPtr r = g->Fit(f1, "EX0RS"); // R = Region   S = StoreResults
    f1->SetLineWidth(1);
    f1->Draw("same");

    gStyle->SetOptFit(1111); //

    TMatrixDSym cov = r->GetCovarianceMatrix();  //  to access the covariance matrix
    double chi2   = r->Chi2(); // to retrieve the fit chi2

    std::cout << "\nCOVARIANCE MATRIX :\n ";
    cov.Print();

    std::cout << "chi2 : " << chi2 << std::endl;
    std::cout << "chi2/dof :" << chi2/(r->Ndf()) << std::endl;


    TGraphErrors *gr = new TGraphErrors(fname); // parte da una copia del grafico originale
    for (int i=0; i<g->GetN(); i++) {
        double res = g->GetY()[i] - f1->Eval(g->GetX()[i]); // residuo
        gr->SetPoint(i,g->GetX()[i],res);
        double eresy = g->GetEY()[i]; // contributo delle Yi
        double eresx = f1->Derivative(g->GetX()[i])*g->GetEX()[i]; // contrib. Xi (approx. 1 ordine)
        double eres = TMath::Sqrt(eresy*eresy+eresx*eresx);
        gr->SetPointError(i,0,eres);
        std::cout << res << " " << eres << std::endl;
    }
    gr->SetName("gr");
    // formattazione del grafico con comandi testuali
    gr->SetMarkerColor(4);
    gr->SetMarkerStyle(20);
    gr->SetTitle("Energy Calibration, residui");
    gr->GetXaxis()->SetTitle("Alpha energy [keV]");
    gr->GetXaxis()->SetTitleSize(0.04);
    //gr->GetXaxis()->CenterTitle();
    gr->GetYaxis()->SetTitle("Integrated charge [a.u.]");
    gr->GetYaxis()->SetTitleSize(0.04);
    //gr->GetYaxis()->CenterTitle();

    TCanvas *c2 = new TCanvas("c2");
    c2->SetGrid();
    //gr->GetXaxis()->SetRangeUser(-1.2,1.2);
    gr->Draw("ap");
    //gr->GetYaxis()->SetRangeUser(-0.8,0.8); // setta il range dell'asse Y

    TGraph *grn = new TGraph(fname); // parte da una copia del grafico originale
    for (int i=0; i<gr->GetN(); i++) {
        double resnorm = gr->GetY()[i]/gr->GetEY()[i]; // residuo
        grn->SetPoint(i,gr->GetX()[i],resnorm);
    }
    grn->SetName("grn");
    // formattazione del grafico con comandi testuali
    grn->SetMarkerColor(4);
    grn->SetMarkerStyle(20);
    grn->SetTitle("Energy Calibration, residui normalizzati");
    grn->GetXaxis()->SetTitle("Alpha energy [keV]");
    grn->GetXaxis()->SetTitleSize(0.04);
    //grn->GetXaxis()->CenterTitle();
    grn->GetYaxis()->SetTitle("res/#sigma_{res} [V/V]");
    grn->GetYaxis()->SetTitleSize(0.04);
    //grn->GetYaxis()->CenterTitle();


    TCanvas *c3 = new TCanvas("c3");
    c3->SetGrid();
    //grn->GetXaxis()->SetRangeUser(-1.2,1.2);
    grn->Draw("ap");
    //grn->GetYaxis()->SetRangeUser(-1,1); // setta il range dell'asse Y
    for (int i=0; i<grn->GetN(); i++){
        //if(grn->GetX()[i]<-1.2 || grn->GetX()[i]>1.2) continue;
        TLine* l = new TLine(grn->GetX()[i],grn->GetY()[i],grn->GetX()[i],0);
        l->SetLineStyle(2);
        l->SetLineWidth(2);
        l->Draw();
    }

}


