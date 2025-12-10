#include "TFile.h"
#include "TH1F.h"
#include "TVector3.h"
#include "TMath.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TStyle.h"
#include "TString.h"
#include "TLegend.h"


/* 

Script to produce histograms for old dark trident samples 


*/ 



void MacroPlot(std::string mass_point){

    std::map<std::string,double> dt_pots = {
        {"0.01",1.68e+23}, {"0.02",2.75e+23},{"0.03",3.95e+23},{"0.04",6.65e+23},{"0.05",1.24e+24},
        {"0.06",2.60e+24}, {"0.07",6.60e+24},{"0.08",2.24e+25},{"0.09",1.38e+26},{"0.1",2.78e+27}
        
    };

    std::map<std::string,double> xsec_correction = {
        {"0.01",28.13}, {"0.02",4.3},{"0.03",1.22},{"0.04",0.5},{"0.05",0.248},
        {"0.06",0.16}, {"0.07",0.14},{"0.08",0.20},{"0.09",0.59},{"0.1",6.3}

    };

    std::map<std::string,double> fermion_scalings ={
        {"0.01",12.7}, {"0.02",14.02},{"0.03",16.9},{"0.04",18.58},{"0.05",11.19},
        {"0.06",20.39}, {"0.07",18.76},{"0.08",15.71},{"0.09",13.3},{"0.1",12.3}

    };



    // Variable to plot and cuts
    std::cout << "Processing mass point : " << mass_point << std::endl;
    // cuts used in the first iteration 
    std::string cut_vars = "n_vertex == 1 && n_showers >= 1 && shower_dir >= 0.99 && signal_score >= 0.9";
    std::string plot_var = "E_total"; 
    std::string parent_directory = "/home/lmlepin/Desktop/dm_sets/dark_tridents_analysis/";
    std::string out_image_string = parent_directory + "old_samples/energy_distribution_" + mass_point + ".png";
    std::string output_root_files= parent_directory + "old_samples/dt_hists_mass_" + mass_point+".root";
    std::cout<<"Test output image :" << out_image_string << std::endl;
    const char *out_img_name = out_image_string.c_str();
    bool is_scalar = false; 
    bool in_progress_option=false;
    bool log_option=true;
    bool save_hist_option=true; 
    bool save_plot=true;
    char y_label[128] = "Number of events after 2e21 POT";
    char x_label[128] = "Energy [MeV]";
    // Get this from function, for signal provide this in config file (?)
    double pot_standard = 6.15622e+21; // This corresponds to pi0 
    double pot_dirt = 3.8386e+23; // corresponds to etas 
    double pot_ext = 4.35574e+23; // corresponds to gammas 
    double pot_dt; 
    double expected_pot = 2.e21; // 5 runs POT 

    if(is_scalar){
        pot_dt = dt_pots[mass_point]/(xsec_correction[mass_point]*4); // Includes fermion scaling and fixing factor for the geneator
    }
    else{
        pot_dt = dt_pots[mass_point]/(fermion_scalings[mass_point]*xsec_correction[mass_point]*4); // Includes fermion scaling and fixing factor for the geneator
    }

    // Set gStyle parameters and canvas dimensions:
    gROOT->ForceStyle(1);
    gStyle->SetOptStat(0); // Turn off statistics chart
    gStyle->SetPadTopMargin(0.07); gStyle->SetPadBottomMargin(0.16);
    gStyle->SetPadRightMargin(0.05); gStyle->SetPadLeftMargin(0.15);
    gStyle->SetHistLineColor(kBlack); gStyle->SetHistLineWidth(1);
    gStyle->SetLegendTextSize(0.035); 


    int canv_x=1000;
    int canv_y=700;
    double y_min=0.1;
    double y_max=10000.;
    double x_min=0.;
    double x_max=1500.; 
    int nbins=5;
    double min_val = 0.;
    double max_val = 1500.;

    auto fileName_neutrino = parent_directory + "NCPi0_info_mpid_scores_8441_steps_EXT.csv";
    auto fileName_b = parent_directory + "NCeta_info_mpid_scores_8441_steps_EXT.csv";
    auto fileName_c = parent_directory + "NCgamma_info_mpid_scores_8441_steps_EXT.csv";
    auto fileName_signal = parent_directory + "dt_overlay_"+ mass_point +"_CNN_scores_8441_steps.csv";
    auto df_standard = ROOT::RDF::MakeCsvDataFrame(fileName_neutrino);
    auto df_dirt = ROOT::RDF::MakeCsvDataFrame(fileName_b);
    auto df_ext = ROOT::RDF::MakeCsvDataFrame(fileName_c);
    auto df_dt = ROOT::RDF::MakeCsvDataFrame(fileName_signal);


    auto total_dt = *df_dt.Count()*(expected_pot/(pot_dt));
    auto total_standard = *df_standard.Count()*(expected_pot/pot_standard);
    auto total_dirt = *df_dirt.Count()*(expected_pot/pot_dirt);
    auto total_ext = *df_ext.Count()*(expected_pot/pot_ext);

    std::cout << "Total simulated DT: " << total_dt << std::endl;
    std::cout << "Total simulated Pi0: " << total_standard << std::endl;
    std::cout << "Total simulated eta: " << total_dirt << std::endl;
    std::cout << "Total simulated gamma: " << total_ext << std::endl;

    std::cout << "===============================" << std::endl;

    auto filteredEvents_standard = df_standard.Filter(cut_vars);
    auto filteredEvents_dirt = df_dirt.Filter(cut_vars);
    auto filteredEvents_ext = df_ext.Filter(cut_vars);
    auto filteredEvents_dt = df_dt.Filter(cut_vars);

    auto cut_dt = *filteredEvents_dt.Count()*(expected_pot/pot_dt);
    auto cut_ncpi0 = *filteredEvents_standard.Count()*(expected_pot/pot_standard);
    auto cut_nceta = *filteredEvents_dirt.Count()*(expected_pot/pot_dirt);
    auto cut_ncgamma = *filteredEvents_ext.Count()*(expected_pot/pot_ext);

    std::cout << "Total a/c DT: " << cut_dt << std::endl;
    std::cout << "Total a/c Pi0: " << cut_ncpi0 << std::endl;
    std::cout << "Total a/c eta: " << cut_nceta << std::endl;
    std::cout << "Total a/c gamma: " << cut_ncgamma << std::endl;

    std::cout << "===============================" << std::endl;


    auto intCol = filteredEvents_standard.Take<double>(plot_var);
    auto intColB = filteredEvents_dirt.Take<double>(plot_var);
    auto intColC = filteredEvents_ext.Take<double>(plot_var);
    auto intColD = filteredEvents_dt.Take<double>(plot_var);

    // Create histograms
    TH1D* h1 = new TH1D("ncpi0"," ", nbins, min_val, max_val);
    TH1D* h2 = new TH1D("nceta"," ", nbins, min_val, max_val);
    TH1D* h3 = new TH1D("ncgamma"," ", nbins, min_val, max_val);
    TH1D* h4 = new TH1D("signal"," ", nbins, min_val, max_val);
    TH1D* background_histo = new TH1D("background","background",nbins,min_val,max_val);


    // Fill histograms
    for (double n : intCol) {
        h1->Fill(n);
    }

    
    for (double n: intColB){
        h2->Fill(n);
    }

    
    for (double n: intColC){
        h3->Fill(n);
    }

    for (double n: intColD){
        h4->Fill(n);
    }


    auto c = new TCanvas();
    c->SetCanvasSize(canv_x,canv_y);

    THStack * hs1 = new THStack("hs1","Distribution");
    TPaveText* ptext = new TPaveText(0.4, 0.83, 0.5, 0.87, "NDC");
    TText *t1 = ptext->AddText("MicroBooNE simulation, preliminary");
    ptext->SetFillColor(kWhite);
    ptext->SetTextFont(53); 
    ptext->SetTextSize(40);
    ptext->SetBorderSize(0);

    // Define color palette 
    Int_t palette[4];
    palette[0] = 920;
    palette[1] = 880 - 4;
    palette[2] = 860 -7;
    palette[3] = 860 + 6;
    gStyle->SetPalette(4,palette);


    h1->Scale(expected_pot/pot_standard);
    h2->Scale(expected_pot/pot_dirt);
    h3->Scale(expected_pot/pot_ext);
    h4->Scale(expected_pot/pot_dt);


    background_histo->Add(h1);
    background_histo->Add(h2);
    background_histo->Add(h3);

    hs1->Add(h4);
    hs1->Add(h3);
    hs1->Add(h2);
    hs1->Add(h1);


    std::cout<<"On histogram: " << std::endl;
    std::cout<<"Number of DM events: " << h4->Integral() << std::endl;
    std::cout<<"Number of NCpi0: " << h1->Integral() << std::endl;
    std::cout<<"Number of NCeta: " << h2->Integral() << std::endl;
    std::cout<<"Number of NCgamma: " << h3->Integral() << std::endl;

    float n_dt = std::floor(h4->Integral());
    float n_standard = std::floor(h1->Integral());
    float n_dirt  = std::floor(h2->Integral());
    float n_ext = h3->Integral(); 
    
    std::stringstream stream_dt;
    stream_dt<< std::fixed << std::setprecision(2) << n_dt;
    std::string dt_string = stream_dt.str();

    std::stringstream stream_standard;
    stream_standard<< std::fixed << std::setprecision(2) << n_standard;
    std::string standard_string = stream_standard.str();

    std::stringstream stream_dirt;
    stream_dirt << std::fixed << std::setprecision(2) << n_dirt;
    std::string dirt_string = stream_dirt.str();

    std::stringstream stream_ext;
    stream_ext << std::fixed << std::setprecision(2) << n_ext;
    std::string ext_string = stream_ext.str(); 

    std::string n_dt_string = "Dark trident: "+ dt_string;
    std::string n_standard_string =  "Pi0: "+ standard_string; 
    std::string n_dirt_string = "eta: "+ dirt_string;
    std::string n_ext_string = "gamma: "+ ext_string;
    const char *dt_char = n_dt_string.c_str();
    const char *standard_char = n_standard_string.c_str();
    const char *dirt_char = n_dirt_string.c_str();
    const char *ext_char = n_ext_string.c_str();
    
   /* else{
        const char *dt_char = "Dark trident";
        const char *standard_char = "Neutrino";
        const char *dirt_char = "Dirt";
        //const char *ext_char = "Cosmics";
    }*/ 

    hs1->SetMinimum(y_min);
    hs1->SetMaximum(y_max);
    hs1->SetTitle("");
    hs1->Draw("pfc HIST");
    auto legend = new TLegend(0.15,0.94,0.98,0.98); // Legend on top of the plot
    legend->SetNColumns(4);
    // legend->SetHeader("Datasets","C"); // option "C" allows to center the header
    // legend->AddEntry(hsum,"Total");
    // legend->AddEntry(hsum_back,"Total background");
    legend->AddEntry(h1,standard_char,"f");
    legend->AddEntry(h2,dirt_char ,"f");
    legend->AddEntry(h3,ext_char,"f");
    legend->AddEntry(h4,dt_char,"f");
    legend->Draw("SAME");

    if(in_progress_option){
       ptext->Draw("SAME");

    }
    if(log_option){
        gPad->SetLogy();

    } 
    
    hs1->GetXaxis()->SetLabelSize(0.04);
    hs1->GetXaxis()->SetTitleOffset(0.9);
    hs1->GetXaxis()->SetTitleSize(0.045);
    hs1->GetYaxis()->SetLabelSize(0.04);
    hs1->GetYaxis()->SetTitleOffset(1.12);
    hs1->GetYaxis()->SetTitleSize(0.045);
    hs1->GetXaxis()->SetTitle(x_label);
    hs1->GetYaxis()->SetTitle(y_label);
    gPad->Modified(); 



    if(save_plot){
        c->SaveAs(out_img_name);
    }

    if(save_hist_option){
        // Create ROOT file to store histograms 
        TFile* out_file = new TFile(output_root_files.c_str(),"RECREATE");
        h4->GetXaxis()->SetTitle(x_label);
        h4->GetYaxis()->SetTitle(y_label);
        background_histo->GetXaxis()->SetTitle(x_label);
        background_histo->GetYaxis()->SetTitle(y_label);
        h4->Write();
        background_histo->Write();
        out_file->Close();
    }

}


