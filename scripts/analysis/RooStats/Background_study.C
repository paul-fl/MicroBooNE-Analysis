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



void StudyBackgrounds(){
    string base_dir = "/home/lmlepin/Desktop/dm_sets/";
    string file = "run1_NuMI_CV_weights_sp.root";
    string input_file = base_dir + file; 
    TFile* f = TFile::Open(input_file.c_str());
    TTree* vertex_tree = new TTree;

    vertex_tree = (TTree *) f->Get("singlephotonana/vertex_tree");
    cout << "Number of entries: " << vertex_tree->GetEntries() << endl; 

    int n_vertex, n_showers, n_tracks, mctruth_nu_pdg, mctruth_lepton_pdg,
    mctruth_cc_or_nc, mctruth_num_exiting_pi0; 
    vector<double> *shower_energies = nullptr; 
    vector<double> *track_energies = nullptr;
    vector<int>  *mctruth_daughters_pdg = nullptr;

    vertex_tree->SetBranchAddress("reco_vertex_size",&n_vertex);
    vertex_tree->SetBranchAddress("reco_asso_showers",&n_showers);
    vertex_tree->SetBranchAddress("reco_asso_tracks", &n_tracks);
    vertex_tree->SetBranchAddress("reco_shower_energy_max", &shower_energies);
    vertex_tree->SetBranchAddress("reco_track_calo_energy_max", &track_energies);
    vertex_tree->SetBranchAddress("mctruth_nu_pdg", &mctruth_nu_pdg);
    vertex_tree->SetBranchAddress("mctruth_lepton_pdg", &mctruth_lepton_pdg);
    vertex_tree->SetBranchAddress("mctruth_cc_or_nc", &mctruth_cc_or_nc); // Zero is CC and One is NC (?) 
    vertex_tree->SetBranchAddress("mctruth_daughters_pdg", &mctruth_daughters_pdg);
    vertex_tree->SetBranchAddress("mctruth_num_exiting_pi0", &mctruth_num_exiting_pi0);


    double total_cc = 0.;
    double total_nc = 0.;
    double cc_numus = 0.;
    double nc_numus = 0.;
    double n_nues = 0.;
    double n_numus = 0.; 
    double cc_nues = 0.;
    double nc_nues = 0.;
    double total_evts = 0.; 

    double total_pi0 = 0.;
    double total_etas = 0.;
    double total_gammas = 0.; 
    

    for(int i = 0; i < vertex_tree->GetEntries() ; i++){
        vertex_tree->GetEntry(i);

        /*
        cout << "Entry number: " << i << endl; 
        cout << "Nu PDG: " << mctruth_nu_pdg << endl; 
        cout << "Lepton PDGL " << mctruth_lepton_pdg << endl; 
        cout << "CC or NC: " << mctruth_cc_or_nc << endl;
        cout << "Number of daughters: " << mctruth_daughters_pdg->size() << endl;

        */ 

        if(n_showers == 1 && n_tracks == 0){

        if(mctruth_nu_pdg == 14 || mctruth_nu_pdg == -14){
            n_numus+=1.;

            if(mctruth_cc_or_nc == 0){
                total_cc += 1.;
                cc_numus += 1.; 
            }

            else{
                total_nc += 1.;
                nc_numus += 1.; 
            }
        }

        else if(mctruth_nu_pdg == 12 || mctruth_nu_pdg == -12){
            n_nues+=1.;

            if(mctruth_cc_or_nc == 0){
                total_cc += 1.;
                cc_nues += 1.; 
            }

            else{
                total_nc += 1.;
                nc_nues += 1.; 
            }



        }

        if(mctruth_num_exiting_pi0 >= 1){
            total_pi0 += 1.;
        }

        

        total_evts += 1.;
        }


    }

    cout << "\n" << endl;


    cout << "Total number of events after cut: " << total_evts << endl; 
    cout << "Fraction of events that pass the cut: " << total_evts/vertex_tree->GetEntries() << endl; 

    cout << "\n" << endl; 


    cout << "Total CC events: " << total_cc << " total NC events: " << total_nc << endl;
    cout << "Fraction of CC events: " << total_cc/total_evts << " fraction of NC events: " << total_nc/total_evts << endl;



    cout << "\n" << endl; 
    cout << "Total number of numus: " << n_numus << ", Fraction from total events: " << n_numus/total_evts << endl;
    cout << "Total CC numus: " << cc_numus << " total NC numus: " << nc_numus << endl;
    cout << "Fraction of CC events: " << cc_numus/n_numus << " fraction of NC events: " << nc_numus/n_numus << endl;

    cout << "\n" << endl;

    cout <<"Total number of nues: " << n_nues << ", Fraction from total events: " << n_nues/total_evts << endl;
    cout << "Total CC nues: " << cc_nues << " total NC nues: " << nc_nues << endl;
    cout << "Fraction of CC events: " << cc_nues/n_nues << " fraction of NC events: " << nc_nues/n_nues << endl;


    cout << "\n" << endl;
    cout << "Events with one pi zero: " << total_pi0 << endl;


} 


void PlotVariables(){

    string base_dir = "/home/lmlepin/Desktop/dm_sets/";
    string file = "run1_NuMI_CV_weights_sp.root";
    string input_file = base_dir + file; 
    TFile* f = TFile::Open(input_file.c_str());
    TTree* vertex_tree = new TTree;

    vertex_tree = (TTree *) f->Get("singlephotonana/vertex_tree");
    cout << "Number of entries: " << vertex_tree->GetEntries() << endl; 

    TCut c1 = "reco_asso_showers >= 0 && reco_asso_tracks >= 0";

    TCut numu_cc = "reco_vertex_size == 1 && mctruth_cc_or_nc == 0 && mctruth_nu_pdg == 14";
    TCut antinumu_cc = "reco_vertex_size == 1 && mctruth_cc_or_nc == 0 && mctruth_nu_pdg == -14";
    TCut mu_cc = (numu_cc||antinumu_cc)&&c1;

    TCut numu_nc = "reco_vertex_size == 1 && mctruth_cc_or_nc == 1 && mctruth_nu_pdg == 14";
    TCut antinumu_nc = "reco_vertex_size == 1 && mctruth_cc_or_nc == 1 && mctruth_nu_pdg == -14";
    TCut mu_nc = (numu_nc||antinumu_nc)&&c1;


    TCut nue_cc = "reco_vertex_size == 1 && mctruth_cc_or_nc == 0 && mctruth_nu_pdg == 12";
    TCut antinue_cc = "reco_vertex_size == 1 && mctruth_cc_or_nc == 0 && mctruth_nu_pdg == -12";
    TCut e_cc = (nue_cc||antinue_cc)&&c1;

    TCut nue_nc = "reco_vertex_size == 1 && mctruth_cc_or_nc == 1 && mctruth_nu_pdg == 12";
    TCut antinue_nc = "reco_vertex_size == 1 && mctruth_cc_or_nc == 1 && mctruth_nu_pdg == -12";
    TCut e_nc = (nue_nc||antinue_nc)&&c1;

    TCanvas *c = new TCanvas();

    vertex_tree->Draw("reco_shower_energy_max >> E_dist(80, 0., 1500)",c1);
    TH1F *energy_dist = (TH1F*) gDirectory->Get("E_dist");
    energy_dist->Draw("HIST"); 

    TCanvas *c2 = new TCanvas();
    gStyle->SetPalette(kBlueRedYellow);

    vertex_tree->Draw("reco_shower_energy_max >> E_dist_mu_cc(50, 0., 1500)",mu_cc);
    vertex_tree->Draw("reco_shower_energy_max >> E_dist_e_cc(50, 0., 1500)",e_cc);
    vertex_tree->Draw("reco_shower_energy_max >> E_dist_mu_nc(50, 0., 1500)",mu_nc);
    vertex_tree->Draw("reco_shower_energy_max >> E_dist_e_nc(50, 0., 1500)",e_nc);


    TLegend *lg = new TLegend();
    THStack *hs = new THStack("hs","Run 1; [MeV]; A.U");

    TH1F *h1 = (TH1F*) gDirectory->Get("E_dist_mu_cc");
    TH1F *h2 = (TH1F*) gDirectory->Get("E_dist_e_cc");
    TH1F *h3 = (TH1F*) gDirectory->Get("E_dist_mu_nc");
    TH1F *h4 = (TH1F*) gDirectory->Get("E_dist_e_nc");

    h1->SetFillColor(kRed);
    lg->AddEntry(h1,"#nu_{#mu} CC");
    hs->Add(h1);
    h3->SetFillColor(kGreen);
    lg->AddEntry(h3,"#nu_{#mu} NC");
    hs->Add(h3);
    h2->SetFillColor(kBlue);
    lg->AddEntry(h2,"#nu_{e} CC");
    hs->Add(h2);
    h4->SetFillColor(kMagenta);
    lg->AddEntry(h4,"#nu_{e} NC");
    hs->Add(h4);
    hs->Draw("pfc"); 
    lg->Draw("SAME");
}



void PlotTopologies(){
    string base_dir = "/home/lmlepin/Desktop/dm_sets/";
    string file = "run1_NuMI_CV_weights_sp.root";
    string input_file = base_dir + file; 
    TFile* f = TFile::Open(input_file.c_str());
    TTree* vertex_tree = new TTree;

    vertex_tree = (TTree *) f->Get("singlephotonana/vertex_tree");
    cout << "Number of entries: " << vertex_tree->GetEntries() << endl; 

    gStyle->SetOptStat(0);

    TCanvas *c = new TCanvas();

    vertex_tree->Draw("reco_asso_showers:reco_asso_tracks >> top_dist(7, 0., 7.,7,0.,7.)");
    TH2F *topology_dist = (TH2F*) gDirectory->Get("top_dist");
    topology_dist->SetTitle("Run1 neutrino background topologies");
    topology_dist->Scale(1./vertex_tree->GetEntries());
    topology_dist->GetXaxis()->SetTitle("Number of tracks");
    topology_dist->GetYaxis()->SetTitle("Number of showers");
    topology_dist->Draw("COLZ TEXT");
}




void PlotCut(){
    string base_dir = "/home/lmlepin/Desktop/dm_sets/";
    string file = "run1_dirt_NuMI_sp_weights.root";
    string input_file = base_dir + file; 
    TFile* f = TFile::Open(input_file.c_str());
    TTree* vertex_tree = new TTree;

    vertex_tree = (TTree *) f->Get("singlephotonana/vertex_tree");
    cout << "Number of entries: " << vertex_tree->GetEntries() << endl; 

    TCut c1 = "reco_asso_showers == 1 && reco_asso_tracks == 0";

    TCanvas *c = new TCanvas();

    vertex_tree->Draw("reco_shower_energy_max >> E_dist(20, 0., 1500)",c1);
    TH1F *energy_dist = (TH1F*) gDirectory->Get("E_dist");
    energy_dist->Draw("HIST"); 

}

