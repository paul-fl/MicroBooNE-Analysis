#include "TFile.h"
#include "TH1D.h"
#include <string> 




void modify_histograms(std::string run){
    char input_file[180];
    sprintf(input_file, "~/Desktop/dm_sets/dark_tridents_analysis/%s_samples/dark_trident_hist_for_collie_%s_CNN.root",run.c_str(),run.c_str());
    std::cout << input_file << std::endl;

    TFile *f = new TFile(input_file,"UPDATE");

    // CV histograms 
    std::vector<TH1D*> CV_vec; 
    TH1D *h_overlay = (TH1D*) f->Get("bkg_overlay");
    TH1D *h_dirt = (TH1D*) f->Get("bkg_dirt");
    TH1D *h_EXT = (TH1D*) f->Get("bkg_EXT");
    TH1D *h_data = (TH1D*) f->Get("data");
    TH1D *h_signal = (TH1D*) f->Get("signal");
    CV_vec.push_back(h_overlay);
    CV_vec.push_back(h_dirt);
    CV_vec.push_back(h_EXT);
    CV_vec.push_back(h_signal);


    // Uncertainties 
    std::vector<TH1D*> uncert_vec;
    TH1D *h_overlay_stat = (TH1D*) f->Get("overlay_stat_uncertainty");
    TH1D *h_dirt_stat = (TH1D*) f->Get("dirt_stat_uncertainty");
    TH1D *h_ext_stat = (TH1D*) f->Get("EXT_stat_uncertainty");
    TH1D *h_signal_stat = (TH1D*) f->Get("signal_stat_uncertainty");
    uncert_vec.push_back(h_overlay_stat);
    uncert_vec.push_back(h_dirt_stat);
    uncert_vec.push_back(h_ext_stat);
    uncert_vec.push_back(h_signal_stat);


    std::vector<TH1D*> multisim_frac;
    TH1D *h_ppfx_fraction = (TH1D*) f->Get("ppfx_uncertainty_frac");
    TH1D *h_genie_fraction = (TH1D*) f->Get("Genie_uncertainty_frac");
    TH1D *h_reint_fraction = (TH1D*) f->Get("Reinteraction_uncertainty_frac");
    multisim_frac.push_back(h_ppfx_fraction);
    multisim_frac.push_back(h_genie_fraction);
    multisim_frac.push_back(h_reint_fraction);


    std::vector<TH1D*> up_vec;
    std::vector<TH1D*> down_vec;


    TH1D *h_ppfx_up = (TH1D*) h_overlay->Clone("ppfx_up");
    up_vec.push_back(h_ppfx_up);
    TH1D *h_ppfx_down = (TH1D*) h_overlay->Clone("ppfx_down");
    down_vec.push_back(h_ppfx_down);
    TH1D *h_genie_up = (TH1D*) h_overlay->Clone("Genie_up");
    up_vec.push_back(h_genie_up);
    TH1D *h_genie_down = (TH1D*) h_overlay->Clone("Genie_down");
    down_vec.push_back(h_genie_down);
    TH1D *h_reint_up = (TH1D*) h_overlay->Clone("Reint_up");
    up_vec.push_back(h_reint_up);
    TH1D *h_reint_down = (TH1D*) h_overlay->Clone("Reint_down");
    down_vec.push_back(h_reint_down);

    
    // Up multisim variations
    for(int j{0}; j<3; j++){
        std::cout << up_vec[j]->GetName() << std::endl;
        for(int i{0}; i < up_vec[j]->GetNbinsX(); i++){
            double new_content = up_vec[j]->GetBinContent(i+1) + up_vec[j]->GetBinContent(i+1)*multisim_frac[j]->GetBinContent(i+1);

            /*
            std::cout << "Old: " << up_vec[j]->GetBinContent(i+1) << std::endl;;
            std::cout << "New: " << new_content << std::endl;
            */
           up_vec[j]->SetBinContent(i+1,new_content);
        }

    }


    // Down multisim variations 
    for(int j{0}; j<3; j++){
        std::cout << down_vec[j]->GetName() << std::endl;
        for(int i{0}; i < down_vec[j]->GetNbinsX(); i++){
            double new_content = down_vec[j]->GetBinContent(i+1) - down_vec[j]->GetBinContent(i+1)*multisim_frac[j]->GetBinContent(i+1);

            /*
            std::cout << "Old: " << up_vec[j]->GetBinContent(i+1) << std::endl;;
            std::cout << "New: " << new_content << std::endl;
            */
           down_vec[j]->SetBinContent(i+1,new_content);
        }

    }


    up_vec[0]->Write("ppfx_up", TObject::kOverwrite);
    up_vec[1]->Write("Genie_up", TObject::kOverwrite);
    up_vec[2]->Write("Reinteraction_up", TObject::kOverwrite);

    down_vec[0]->Write("ppfx_down", TObject::kOverwrite);
    down_vec[1]->Write("Genie_down", TObject::kOverwrite);
    down_vec[2]->Write("Reinteraction_down", TObject::kOverwrite);



    // Setting statistical error for MC samples 
    for(int j{0}; j < 4; j++){
        std::cout << CV_vec[j]->GetName() << std::endl;
        for(int i{0}; i < CV_vec[j]->GetNbinsX(); i++){
            CV_vec[j]->SetBinError(i+1,uncert_vec[j]->GetBinContent(i+1));
        }
    }

    // Setting stats error for data
    for(int i{0}; i < h_data->GetNbinsX(); i++){
        h_data->SetBinError(i+1, TMath::Sqrt(h_data->GetBinContent(i+1)));
    }

    CV_vec[0]->Write("bkg_overlay",TObject::kOverwrite);
    CV_vec[1]->Write("bkg_dirt",TObject::kOverwrite);
    CV_vec[2]->Write("bkg_EXT",TObject::kOverwrite);
    CV_vec[3]->Write("signal",TObject::kOverwrite);
    h_data->Write("data", TObject::kOverwrite);





    f->Close();


    
}