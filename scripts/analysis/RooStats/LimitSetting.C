//Root Headers
#include <TCanvas.h>
#include <TMinuit.h>
#include <TGraph2DErrors.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TF2.h>
#include <TTree.h>
#include <TFile.h>
#include <TSystem.h>
#include <TGraphAsymmErrors.h>

#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooRandom.h"

#include "RooStats/ModelConfig.h"


#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/HybridCalculator.h"
#include "RooStats/FrequentistCalculator.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/HypoTestPlot.h"

#include "RooStats/NumEventsTestStat.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/SimpleLikelihoodRatioTestStat.h"
#include "RooStats/RatioOfProfiledLikelihoodsTestStat.h"
#include "RooStats/MaxLikelihoodEstimateTestStat.h"
#include "RooStats/NumEventsTestStat.h"

#include "RooStats/HypoTestInverter.h"
#include "RooStats/HypoTestInverterResult.h"
#include "RooStats/HypoTestInverterPlot.h"

#include <stdlib.h>

/* 

   Macro to set limits on the dark trident parameter space
   Based on RooStats tutorials and Anyssa's code for Dark Side  

*/ 



using namespace RooStats;
using namespace RooFit;
using namespace std;




HypoTestInverterResult* SimpleHypoTestInv( RooWorkspace* w,
                     const char* modelConfigName,
                     const char* dataName );




void PLLTestLimits(){

    // Set gStyle parameters and canvas dimensions:
    gROOT->ForceStyle(1);
    gStyle->SetOptStat(0); // Turn off statistics chart
    gStyle->SetPadTopMargin(0.07); gStyle->SetPadBottomMargin(0.16);
    gStyle->SetPadRightMargin(0.05); gStyle->SetPadLeftMargin(0.15);
    //gStyle->SetHistLineColor(kBlack); gStyle->SetHistLineWidth(1);
    gStyle->SetLegendTextSize(0.035); 


    string parent_directory = "/home/lmlepin/Desktop/dm_sets/dark_tridents_analysis/Official/";
    string old_limits_file = parent_directory + "old_limits.csv";

    // Open limits obtained using single bin counting experiment 
    auto df_old_limits = ROOT::RDF::MakeCsvDataFrame(old_limits_file.c_str());
    auto a01_col = df_old_limits.Take<double>("epsilon_a01");
    auto a1_col = df_old_limits.Take<double>("epsilon_a1");
    vector<string>  masses= {"0.01", "0.02", "0.03", "0.04", "0.05", "0.06", "0.07", "0.08", "0.09", "0.1"};


    TFile *file = NULL;
    TGraphErrors* gobs = new TGraphErrors;
    TGraph * g0 = new TGraph;
    TGraph * g_a01 = new TGraph; 
    TGraphAsymmErrors* g1 = new TGraphAsymmErrors;
    TGraphAsymmErrors* g2 = new TGraphAsymmErrors;
    RooWorkspace* w;
    RooRealVar* nsig_0;
    HypoTestInverterResult *r;
    double mass_point; 
    int counter = 0; 
    for(string m: masses){
      cout << "Processing mass point: " << m << std::endl; 

      string file_name = parent_directory + "stats_space_energy_" + m + ".root";
      file = TFile::Open(file_name.c_str());
      mass_point = atof(m.c_str());
      string ws_name = "w_" + m;
      w = (RooWorkspace*) file->Get(ws_name.c_str());  
      //w->Print(); could be use for check that everyting is in place 
      nsig_0 = (RooRealVar*) w->obj("nsig_0");
      string mc_name = "ModelConfig_" + m;  
      r = SimpleHypoTestInv(w,mc_name.c_str(), "data");

       
      double norm = pow(1e-3,4)/nsig_0->getValV(); // scaling factor to get the respective SM coupling value 
      double upperLimit = TMath::Sqrt(r->UpperLimit()*norm);
      double median = TMath::Sqrt(r->GetExpectedUpperLimit(0)*norm);
      double neg1sig =  TMath::Sqrt(r->GetExpectedUpperLimit(-1)*norm);
      double posi1sig =  TMath::Sqrt(r->GetExpectedUpperLimit(1)*norm);
      double neg2sig =  TMath::Sqrt(r->GetExpectedUpperLimit(-2)*norm);
      double posi2sig =  TMath::Sqrt(r->GetExpectedUpperLimit(2)*norm);
      
      std::cout << "Upper limit: " << upperLimit << std::endl;
      std::cout << "Expected limit: " << median << std::endl;
      std::cout << "Expected limit(+1 sig)" << posi1sig << std::endl;
      std::cout << "Expected limit(-1 sig)" << neg1sig << std::endl;

      gobs->SetPoint(counter, mass_point,  upperLimit);
      g0->SetPoint(counter, mass_point,  median);
      g_a01->SetPoint(counter,mass_point,a01_col->at(counter));
      g1->SetPoint(counter,mass_point,  median);
      g1->SetPointEYlow(counter, median - neg1sig); // -1 sigma errorr
      g1->SetPointEYhigh(counter, posi1sig - median);//+1 sigma error
      g2->SetPoint(counter, mass_point, median);
      g2->SetPointEYlow(counter, median - neg2sig);   // -2 -- -1 sigma error
      g2->SetPointEYhigh(counter, posi2sig - median);
      file->Close(); 
      counter+=1;

    }
  gobs->SetLineWidth(2);
  gobs->SetMarkerStyle(20);
  
  TMultiGraph* graph = new TMultiGraph("Limit","Limit");
  // set the graphics options and add in multi graph
  // orderof adding is drawing order
  g2->SetFillColor(kYellow); graph->Add(g2,"3");
  g1->SetFillColor(kGreen); graph->Add(g1,"3");
  g0->SetLineStyle(2); g0->SetLineWidth(1);
  g_a01->SetLineStyle(2); g_a01->SetLineWidth(1); g_a01->SetLineColor(kBlue);
  graph->Add(g0,"L");
  //graph->Add(g_a01,"L");
  
  TCanvas *cc = new TCanvas("cc","");
  cc->SetTicks(1, 1);
  cc->Draw();
  cc->SetLogx();
  cc->SetLogy();
  

  graph->Draw("A");
  if (graph->GetHistogram()) graph->GetHistogram()->SetTitle( "; M_{A'} [GeV/c^{2}]; #varepsilon^{2}" );
  
  gobs->Draw("SAME PL");
  
  cc->Update();

  string output_plot = parent_directory+"PLL_Limit_DT.pdf";
  //cc->Print("PLL_Limit_DT.eps");
  graph->GetXaxis()->SetLimits(1e-2,1e-1);
  graph->GetYaxis()->SetLimits(1e-11,1e-5);
  graph->GetXaxis()->SetLabelSize(0.04);
  graph->GetXaxis()->SetTitleOffset(0.98);
  graph->GetXaxis()->SetTitleSize(0.045);
  graph->GetYaxis()->SetLabelSize(0.04);
  graph->GetYaxis()->SetTitleOffset(1.12);
  graph->GetYaxis()->SetTitleSize(0.045);
  cc->Print(output_plot.c_str());


  
  gPad->Update();
}





HypoTestInverterResult* SimpleHypoTestInv( RooWorkspace* w,
                     const char* modelConfigName,
                     const char* dataName)
{

  /////////////////////////////////////////////////////////////
  // First part is just to access the workspace file 
  ////////////////////////////////////////////////////////////

  std::cout << "Running HypoTestInverter on the workspace " << w->GetName() << std::endl;
  //w->Print();

  // get the modelConfig out of the file
  RooStats::ModelConfig* mc = (RooStats::ModelConfig*) w->obj(modelConfigName);

  // get the modelConfig out of the file
  RooAbsData* data = w->data(dataName);
  
  //check the data
  if (!data) {
    Error("StandardHypoTestDemo","Not existing data %s",dataName);
    return 0;
  }
  else
    std::cout << "Using data set " << dataName << std::endl;

  ModelConfig*  sbModel = (RooStats::ModelConfig*) w->obj(modelConfigName);
  ModelConfig*  bModel  = (RooStats::ModelConfig*) w->obj(modelConfigName);


  //check the model
  if(!sbModel){
    Error("PLLTestLimits","Not existing ModelConfig %s",modelConfigName);
  }
  if(!sbModel->GetPdf()){
    Error("PLLTestLimits","Model %s has no pdf", modelConfigName);
  }

  if (!sbModel->GetSnapshot() ) {
    Info("PLLTestLimits","Model %s has no snapshot  - make one using model poi",modelConfigName);
    sbModel->SetSnapshot( *sbModel->GetParametersOfInterest() );
  }


  RooRealVar* poi = (RooRealVar*) sbModel->GetParametersOfInterest()->first();
  double orPoi = poi->getVal(); 
  cout << "POI original value: " << poi->GetName() << " = " << poi->getVal() << endl;  


  if (!bModel || bModel == sbModel) {
    Info("PLLTestLimits","The background model %s does not exist",modelConfigName);
    Info("PLLTestLimits","Copy it from ModelConfig %s and set POI to zero",modelConfigName);
    bModel = (ModelConfig*) sbModel->Clone();
    bModel->SetName(TString(modelConfigName)+TString("_with_poi_0"));
    poi->setVal(0);
    bModel->SetSnapshot( *poi );
    poi->setVal(orPoi); 
  }
  else{

    if (!bModel->GetSnapshot() ) {
      Info("PLLTestLimits","B model has no snapshot  - make one using model poi and 0 values ");
      poi->setVal(0);
      bModel->SetSnapshot( *poi);
      poi->setVal(orPoi);
    }
        
  }


  // create HypoTest calculator (data, alt model , null model)

  //FrequentistCalculator  fc(*data, *bModel, *sbModel);
  //fc.SetToys(1000,1000);

  // asymptotic calculator
  AsymptoticCalculator  ac(*data, *bModel, *sbModel, true);


  ac.SetOneSided(true);  // for one-side tests (limits)
  //  ac->SetQTilde(true);
  AsymptoticCalculator::SetPrintLevel(0);

  cout << "Calculator def" << endl;
  
#define UseAsym
  // create hypotest inverter 
  // passing the desired calculator 
#ifdef UseAsym
  HypoTestInverter calc(ac);    // for asymptotic
#else
  HypoTestInverter calc(fc);  // for frequentist
#endif
  // set confidence level (e.g. 95% upper limits)
  calc.SetConfidenceLevel(0.9);

  cout << "CL set" <<endl;

  // for CLS
  bool useCLs = true;
  calc.UseCLs(useCLs);
  calc.SetVerbose(true);

//  calc.SetParameter("EnableDetailedOutput", true);

#ifndef UseAsym
  // configure ToyMC Samler (needed only for frequentit calculator)
  ToyMCSampler *toymcs = (ToyMCSampler*)calc.GetHypoTestCalculator()->GetTestStatSampler();
#endif

  cout <<"Test statistics " << endl;

//  // profile likelihood test statistics
  ProfileLikelihoodTestStat profll(*sbModel->GetPdf());
//  // for CLs (bounded intervals) use one-sided profile likelihood
  if (useCLs) profll.SetOneSided(true);

  //profll.EnableDetailedOutput(); //detailed outputs

#ifndef UseAsym
  // set the test statistic to use 
  toymcs->SetTestStatistic(&profll);

  // if the pdf is not extended (e.g. in the Poisson model) 
  // we need to set the number of events
  if (!sbModel->GetPdf()->canBeExtended())
     toymcs->SetNEventsPerToy(1);
#endif

  int npoints = 100;  // number of points to scan
  // min and max (better to choose smaller intervals)
  double poimin = poi->getMin();
  double poimax = 10; //poi->getMax();
  //poimin = 0; poimax=10;

  std::cout << "Doing a fixed scan  in interval : " << poimin << " , " << poimax << std::endl;
  calc.SetFixedScan(npoints,poimin,poimax);
  

  HypoTestInverterResult * r = calc.GetInterval();

  double upperLimit = r->UpperLimit();

  std::cout << "The computed upper limit is: " << upperLimit << std::endl;
  RooRealVar* nsig_0 = (RooRealVar*) w->obj("nsig_0");
  std::cout << "The computed upper limit divide by nominal value: "<< nsig_0->getValV()<<" is: " << upperLimit/nsig_0->getValV() << std::endl;

  // compute expected limit
  std::cout << "Expected upper limits, using the B (alternate) model : " << std::endl;
  std::cout << " expected limit (median) " << r->GetExpectedUpperLimit(0) << std::endl;
  std::cout << " expected limit (-1 sig) " << r->GetExpectedUpperLimit(-1) << std::endl;
  std::cout << " expected limit (+1 sig) " << r->GetExpectedUpperLimit(1) << std::endl;
  
//  profll.GetDetailedOutput();

  // plot now the result of the scan 

  HypoTestInverterPlot *plot = new HypoTestInverterPlot("HTI_Result_Plot","HypoTest Scan Result",r);

  // plot in a new canvas with style
  TCanvas * c1 = new TCanvas("HypoTestInverter Scan"); 
  c1->SetLogy(false);

  plot->Draw("CLb 2CL");  // plot also CLb and CLs+b
//  plot->Draw("OBS");  // plot only observed p-value


  // plot also in a new canvas the test statistics distributions 
  
  // plot test statistics distributions for the two hypothesis
  // when distribution is generated (case of FrequentistCalculators)
  const int n = r->ArraySize();
  if (n> 0 &&  r->GetResult(0)->GetNullDistribution() ) { 
     TCanvas * c2 = new TCanvas("Test Statistic Distributions","",2);
     if (n > 1) {
        int ny = TMath::CeilNint( sqrt(n) );
        int nx = TMath::CeilNint(double(n)/ny);
        c2->Divide( nx,ny);
     }
     for (int i=0; i<n; i++) {
        if (n > 1) c2->cd(i+1);
        SamplingDistPlot * pl = plot->MakeTestStatPlot(i);
        pl->SetLogYaxis(true);
        pl->Draw();
        
     }
  }

  return r;
}

