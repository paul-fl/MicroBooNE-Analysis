//Root Headers
#include <TMinuit.h>
#include <TGraph2DErrors.h>
#include <TF2.h>
#include <TTree.h>
#include <TFile.h>
#include <TSystem.h>
#include "TCanvas.h"

#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooRandom.h"

#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooMinuit.h"

#include "RooStats/ModelConfig.h"

// For HistoFactory
#include "RooStats/HistFactory/Measurement.h"
#include "RooStats/HistFactory/MakeModelAndMeasurementsFast.h"
#include "RooStats/HistFactory/FlexibleInterpVar.h"



using namespace RooFit;
using namespace RooStats;
using namespace RooStats::HistFactory;
using namespace std;

void CreateStats(string mass_point){

    string parent_directory = "/home/lmlepin/Desktop/dm_sets/dark_tridents_analysis/Official/";
    string input_file = parent_directory + "dt_hists_mass_"+ mass_point+".root";
    TFile* f = TFile::Open(input_file.c_str());
    TH1D* signal = (TH1D*) f->Get("signal");
    TH1D* bkg = (TH1D* ) f->Get("background");

    int nbins = signal->GetNbinsX(); 
    signal->SetDirectory(0);
    bkg->SetDirectory(0);
    f->Close();

    signal->Sumw2();
    bkg->Sumw2();

    std::cout<<"Number of signal events: " << signal->Integral() << std::endl;
    std::cout<<"Number of bkg events: " << bkg->Integral() << std::endl;

    float Nsig_0 = signal->Integral();
    float Nbkg = bkg->Integral();

    // Initialize workspace
    string workspace_name = "w_" + mass_point;
    RooWorkspace w(workspace_name.c_str());

    // Declare observable
    RooRealVar score("score","Signal score",0.,1500.);
    // Create binned dataset 
    RooDataHist data("data","data",score,bkg); //setting it to background for the moment...
    // Import to workspace
    w.import(data);

    // Build model from high stats MC: DM and Bkg
    RooDataHist hsignal("signal","signal",score,signal);
    RooDataHist hbkg("bkg","bkg",score,bkg);
    // Represent signal as pdf over the signal score
    RooHistPdf sigpdf("sigpdf","sigpdf",score,hsignal);
    RooHistPdf bkgpdf("bkgpdf","bkgpdf",score,hbkg);
    w.import(sigpdf);
    w.import(bkgpdf);

    // Number of dark trident events for nominal coupling to SM (epsilon = 0.001)
    RooRealVar nsig_0("nsig_0","nsig_0", Nsig_0);
    nsig_0.setConstant(true);
    w.import(nsig_0);


    // Put a constraint for BG
    w.factory("Gaussian::constraint(b0[0,1000],nbkg[0,1000],sigmab[1])");
    w.factory("SUM:pdf(nsig[0,1000]*sigpdf, nbkg*bkgpdf)");
    w.factory("PROD:model(pdf,constraint)");

    w.var("b0")->setVal(Nbkg);
    w.var("b0")->setConstant(true); 
    w.var("sigmab")->setVal(0.1*Nbkg); // Trying 10% syst uncerainty 

    RooAbsPdf *pdf = w.pdf("model");

    // Set the desired value of signal and background events
    w.var("nsig")->setVal(Nsig_0);
    w.var("nbkg")->setVal(Nbkg);

    RooFitResult *r = pdf->fitTo(data,Save(),SumW2Error(kTRUE));
    r->Print();

    RooPlot* myFrame = score.frame(Title("Signal score distibution"),Bins(nbins));

    
    data.plotOn(myFrame);
    pdf->plotOn(myFrame);
    //draw the two separate pdf's
    pdf->plotOn(myFrame, RooFit::Components("bkgpdf"), RooFit::LineStyle(kDashed) );
    pdf->plotOn(myFrame, RooFit::Components("sigpdf"), RooFit::LineColor(kRed), RooFit::LineStyle(kDashed) );

    string mc_name = "ModelConfig_" + mass_point;
    RooStats::ModelConfig mc(mc_name.c_str(),&w);
    mc.SetPdf(*pdf);
    mc.SetParametersOfInterest(*w.var("nsig"));
    mc.SetObservables(*w.var("score"));
    // Define set of nuisance parameters
    w.defineSet("nuisParams","nbkg");
    mc.SetNuisanceParameters(*w.set("nuisParams"));

    // Needed for hypothesis testing
    mc.SetSnapshot(*w.var("nsig"));
    mc.SetGlobalObservables(*w.var("b0"));

    mc.Print();
    w.import(mc);


    //Write the workspace in root file
    string output_file = parent_directory + "stats_space_energy_"+mass_point+".root";
    w.writeToFile(output_file.c_str(),"RECREATE");



    // Make canvas and draw RooPlots
    TCanvas *c = new TCanvas("c","c",500, 500);
    c->SetLogy();
    myFrame->Draw();
    string out_img_name = parent_directory + "fitted_model_mass_ " + mass_point + ".png";
    c->SaveAs(out_img_name.c_str());
}

void CreateSpaces(){
    vector<string> mass_points = {"0.01", "0.02", "0.03", "0.04", "0.05", "0.06", "0.07", "0.08", "0.09", "0.1"};
    for(string m: mass_points){
        cout << "Proccesing mass point: " << m << endl;
        CreateStats(m);
    }
}