#include "TNtuple.h"
#include "TMath.h"
#include "Riostream.h"


void WidthCalc(TNtuple* nt, float int_l, float int_u, float vmax_l, float vmax_u){

    // variables to read
    float integral;
    float vmax;
    float width;

    // setting addresses to store values into variables
    nt->SetBranchAddress("vmax",&vmax);
    nt->SetBranchAddress("integral", &integral);
    nt->SetBranchAddress("width", &width);

    float sum=0;
    float sum2=0;
    int selected_ev=0;

    for(int i=0; i<nt->GetEntries(); i++){
        nt->GetEntry(i);
        if(integral>int_l && integral<int_u && vmax>vmax_l && vmax<vmax_u){
            sum+=width;
            sum2+=width*width;
            selected_ev+=1;
        }

    }

    float width_avg = sum/selected_ev;
    float dev_std = TMath::Sqrt((sum2-selected_ev*width_avg*width_avg)/(selected_ev-1));
    float err_avg = dev_std/(TMath::Sqrt(selected_ev));

    std::cout << "width = " << width_avg << " #pm " << err_avg << std::endl;
    std::cout << "dev_std = " << dev_std << " selected ev: " << selected_ev << std::endl;


}

