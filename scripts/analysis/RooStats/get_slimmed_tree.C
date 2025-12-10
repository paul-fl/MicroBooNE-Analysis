#include "TTree.h"
#include <iostream>
#include <string>



auto list_branches = {"run_number",
    "subrun_number",
    "event_number",
    "reco_asso_showers",
    "reco_asso_tracks",
    "reco_vertex_x",
    "reco_vertex_y",
    "reco_vertex_z",
    "reco_slice_objects",
    "reco_slice_nuscore",
    "sss_candidate_min_dist",
    "reco_shower_length",
    "reco_shower_opening_angle",
    "reco_shower_startx",
    "reco_shower_starty",
    "reco_shower_startz",
    "reco_shower_start_dist_to_SCB",
    "reco_shower_end_dist_to_SCB",
    "reco_shower_conversion_distance",
    "reco_shower_impact_parameter",
    "reco_shower_energy_max",
    "reco_shower_nuscore",
    "reco_shower_trackscore",
    "reco_shower_theta_yz",
    "reco_shower_phi_yx",
    "reco_shower_opening_angle",
    "reco_shower_dirx",
    "reco_shower_diry",
    "reco_shower_dirz",
    "reco_shower_implied_dirx",
    "reco_shower_implied_diry",
    "reco_shower_implied_dirz",
    "reco_shower_dEdx_plane0_max",
    "reco_shower_dEdx_plane1_max",
    "reco_shower_dEdx_plane2_max",
    "reco_track_dirx",
    "reco_track_diry",
    "reco_track_dirz",
    "reco_track_startx",
    "reco_track_starty",
    "reco_track_startz",
    "reco_track_endx",
    "reco_track_endy",
    "reco_track_endz",
    "reco_track_start_dist_to_SCB",
    "reco_track_end_dist_to_SCB",
    "reco_track_theta_yz",
    "reco_track_phi_yx",
    "reco_track_calo_energy_max",
    "reco_track_proton_kinetic_energy",
    "reco_track_dEdx_plane0",
    "reco_track_dEdx_plane1",
    "reco_track_dEdx_plane2",
    "reco_track_nuscore",
    "reco_track_isclearcosmic",
    "reco_track_is_nuslice" };





void CreateSlimmedTree(TTree *input_tree, std::string output_file){
    std::string output_dir = "/uboone/data/users/lmoralep/sp_files/";
    input_tree->SetBranchStatus("*",0);
    for(auto activeBranchName: list_branches){
        input_tree->SetBranchStatus(activeBranchName,1);
    }
    std::cout << "Creating slimmed TTree " << std::endl;
    TFile newfile((output_dir + output_file).c_str(), "recreate");
    TTree *newtree = input_tree->CopyTree("reco_asso_showers == 1 && reco_asso_tracks == 1");
    //newtree->Print();
    newfile.Write();
    newfile.Close();
    std::cout << "New tree created successfully " << std::endl;
}




int main(){



    std::string base_dir = "/pnfs/uboone/persistent/users/lmoralep/single_photon_files/";
    //std::string base_dir = "/uboone/data/users/lmoralep/sp_files/";
    std::string nu_overlay_file = "dt_overlay_0.05_sp.root";
    TFile *base_file = new TFile((base_dir + nu_overlay_file).c_str());
    TTree *base_tree = (TTree*) base_file->Get("singlephotonana/vertex_tree");
    std::string output_file = "dt_overlay_0.05_slimmed_1shower_1track.root";
    CreateSlimmedTree(base_tree, output_file);

    
    return 0;
}



