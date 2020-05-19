
#include <string.h>
#include "TMatrixDSym.h"
#include "TLine.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TF1.h"
#include "TLatex.h"
#include "TPaveStats.h"

double fitf(double *x, double *par) {
    double fitval = par[0] + par[1]*x[0];
    return fitval;
}

void ParametersPlotting(TCanvas* c, TFitResultPtr r){ //can be improvedS
    double xmax = c->GetUxmax();
    double xmin = c->GetUxmin();
    double ymin = c->GetUymin();
    double ymax = c->GetUymax();

    double xrange = xmax-xmin;
    double yrange = (ymin-ymax);

    TMatrixDSym cov = r->GetCovarianceMatrix();  //  to access the covariance matrix
    double chi2  = r->Chi2(); // to retrieve the fit chi2

    double m = r->Parameter(1);
    double erm = TMath::Sqrt(cov(1,1));
    double q = r->Parameter(0);
    double erq = TMath::Sqrt(cov(0,0));
    double cov_mq = cov(1,0);
    double prob = r->Prob();


    TString s1 = Form("m = %3.3f #pm %3.3f [#mus bar]",m,erm);
    TLatex* t1 = new TLatex(2,3.5,s1);
    t1->SetTextFont(42);
    t1->SetLineWidth(2);
    t1->Draw();
    TString s2 = Form("q = %3.3f #pm %3.3f [#mus]",q,erq);
    TLatex* t2 = new TLatex(2,3.7,s2);
    t2->SetTextFont(42);
    t2->SetLineWidth(2);
    t2->Draw();
    TString s3 = Form("#rho(m,q) = %3.3f",cov_mq);
    TLatex* t3 = new TLatex(2,3.9,s3);
    t3->SetTextFont(42);
    t3->SetLineWidth(2);
    t3->Draw();
    TString s4 = Form("#chi^{2} / ndf = %3.3f / %d",chi2,r->Ndf());
    TLatex* t4 = new TLatex(2,4,s4);
    t4->SetTextFont(42);
    t4->SetLineWidth(2);
    t4->Draw();
    TString s5 = Form("Prob = %3.3f",prob);
    TLatex* t5 = new TLatex(2,4.1,s5);
    t5->SetTextFont(42);
    t5->SetLineWidth(2);
    t5->Draw();

    TPaveStats *ptstats = new TPaveStats(0.6020619,0.6268293,0.8997938,0.8993902,"brNDC");
    ptstats->SetName("stats");
    ptstats->SetBorderSize(1);
    ptstats->SetFillColor(0);
    ptstats->SetFillStyle(0);
    ptstats->SetTextAlign(12);
    ptstats->SetTextFont(42);
    ptstats->SetOptStat(0);
    ptstats->SetOptFit(0);
    ptstats->Draw();

}

void FitWidth(const char* fname) {

    // crea l'oggetto GRAFICO CON ERRORI leggendo da file di testo
    TGraphErrors *g = new TGraphErrors(fname);
    g->SetName("FitWidth");
    // formattazione del grafico con comandi testuali
    g->SetMarkerColor(4);
    g->SetMarkerStyle(20);
    g->SetTitle("Plot Width vs P^{-1} ^{244}Cm");
    g->GetXaxis()->SetTitle("P^{-1} [bar^{-1}]");
    g->GetXaxis()->SetTitleSize(0.04);
    //g->GetXaxis()->CenterTitle();
    g->GetYaxis()->SetTitle("Width [#mus]");
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

    gStyle->SetOptFit(0); //

    ParametersPlotting(c1,r);

    TMatrixDSym cov = r->GetCovarianceMatrix();  //  to access the covariance matrix
    double chi2   = r->Chi2(); // to retrieve the fit chi2

    std::cout << "\nCOVARIANCE MATRIX :\n ";
    cov.Print();

    std::cout << "chi2 : " << chi2 << std::endl;
    std::cout << "chi2/dof :" << chi2/(r->Ndf()) << std::endl;

    float sigmay_post = 0;
    for(int i=0; i<g->GetN(); i++) {
        sigmay_post+= (g->GetY()[i] - f1->Eval(g->GetX()[i]))*(g->GetY()[i] - f1->Eval(g->GetX()[i]));
    }
    sigmay_post = TMath::Sqrt(sigmay_post/(g->GetN()-2));
    std::cout << "SIGMAY " << sigmay_post << std::endl;


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
    gr->SetTitle("Plot Width vs P^{-1} ^{244}Cm, residui");
    gr->GetXaxis()->SetTitle("P^{-1} [bar^{-1}]");
    gr->GetXaxis()->SetTitleSize(0.04);
    //gr->GetXaxis()->CenterTitle();
    gr->GetYaxis()->SetTitle("Width [#mus]");
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
    grn->SetTitle("Plot Width vs P^{-1} ^{244}Cm, residui normalizzati");
    grn->GetXaxis()->SetTitle("P^{-1} [bar^{-1}]");
    grn->GetXaxis()->SetTitleSize(0.04);
    //grn->GetXaxis()->CenterTitle();
    grn->GetYaxis()->SetTitle("res/#sigma_{res} [#mus/#mus]");
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



