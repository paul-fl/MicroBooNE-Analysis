#include "TFile.h"
#include "TH1F.h"
#include "TVector3.h"
#include "TMath.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TStyle.h"
#include "TString.h"
#include "TLegend.h"

using namespace std; 

void Test(){
    gStyle->SetOptStat(0);

    string base_dir = "/home/lmlepin/Desktop/dm_sets/";
    string file = "run1_standard_NuMI_weights_sp.root";
    string input_file = base_dir + file; 
    TFile* f = TFile::Open(input_file.c_str());
    TTree* vertex_tree = new TTree;
    TTree* weight_tree = new TTree;
    vertex_tree = (TTree *) f->Get("singlephotonana/vertex_tree");
    weight_tree = (TTree *) f->Get("singlephotonana/EventWeights");



    cout << "Number of entries: " << vertex_tree->GetEntries() << endl; 

    int n_vertex, n_showers, n_tracks; 
    vector<double> *shower_energies = nullptr; 
    vector<double> *track_energies = nullptr;

    vertex_tree->SetBranchAddress("reco_vertex_size",&n_vertex);
    vertex_tree->SetBranchAddress("reco_asso_showers",&n_showers);
    vertex_tree->SetBranchAddress("reco_asso_tracks", &n_tracks);
    vertex_tree->SetBranchAddress("reco_shower_energy_max", &shower_energies);
    vertex_tree->SetBranchAddress("reco_track_calo_energy_max", &track_energies);

    vector<short> *ppfx_weights = nullptr;
    weight_tree->SetBranchAddress("weightsReint", &ppfx_weights);




    TH1D *energy_dist = new TH1D("E", "Distribution; E [MeV]; events", 10, 0., 1500);
    TH2D *density_hists = new TH2D("Den", "Reinteraction weights; E[MeV]; events", 10, 0., 1500., 30, 0., 800.);

    /*
    TH1D *reweighted_0 =  new TH1D("Re0", "Distribution; E [MeV]; events", 30, 0., 3000);
    TH1D *reweighted_1 =  new TH1D("Re1", "Distribution; E [MeV]; events", 30, 0., 3000);

    reweighted_0->SetLineColor(kRed);
    reweighted_1->SetLineColor(kMagenta);
    */ 

   TH1D *reweighted_hists[600];
   char name[20];
   char title[100];
   for (Int_t i=0;i<600;i++) {
      sprintf(name,"Re%d",i);
      sprintf(title,"Distribution h%d",i);
      reweighted_hists[i] = new TH1D(name,title,30,0.,3000);
   }


    for(int i = 0; i < vertex_tree->GetEntries() ; i++){
        double total_E = 0.; 
        vertex_tree->GetEntry(i);
        weight_tree->GetEntry(i); 

        /*
        cout << "Number of reco vertex: " << n_vertex << endl;
        cout << "Number of reco showers: " << n_showers << endl; 
        cout << "Number of reco tracks: " << n_tracks << endl;
        */ 


        
        if(n_showers >= 1){
            for(int j = 0; j < n_showers; j++){
                total_E+=shower_energies->at(j);
            }
        }

        if(n_tracks >= 1){
            for(int j = 0; j < n_tracks; j++){
                total_E+=track_energies->at(j);
            }

        }

        if(n_vertex ==1 && (n_showers == 1 && n_tracks == 0)){
            energy_dist->Fill(total_E);
            for(int i = 0; i < 600; i++){
                reweighted_hists[i]->Fill(total_E, ppfx_weights->at(i)/1000.);
            }
            //cout << "Weight vector" << ppfx_weights->size() << endl;
            //reweighted_1->Fill(total_E, ppfx_weights->at(1)/1000.);

        }
        //cout << "Total shower energy: " << total_E << endl;

    }




    for(int i = 0; i < 600; i++){
        for(int j = 0; j < reweighted_hists[i]->GetNbinsX() ; j++ ){
            density_hists->Fill(reweighted_hists[i]->GetBinCenter(j+1), reweighted_hists[i]->GetBinContent(j+1));
        }
    }


    density_hists->Draw("COLZ");
    energy_dist->SetLineColor(kRed);
    energy_dist->SetLineWidth(5);
    energy_dist->Draw("HIST SAME");

    /*
    reweighted_hists[0]->Draw("HIST SAME");
    reweighted_hists[1]->Draw("HIST SAME");*/ 



}