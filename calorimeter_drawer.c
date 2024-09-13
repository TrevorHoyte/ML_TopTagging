/* call macro:
root -l mymacros/histogram2d.C'("delphes_output.root")'
*/
#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"

#else

class ExRootTreeReader;
class ExRootResult;

#endif
//------------------------------------------------------------------------------
struct TestPlots
{
 
  TH2 *fJet;
  TH2 *Jetimage1;
  TH2 *Jetimage2;
  TH2 *Jetimage3;



};
//------------------------------------------------------------------------------
void BookHistograms(ExRootResult *result, TestPlots *plots)
{
  TLegend *legend;
  TPaveText *comment;
 
  plots->fJet = result->AddHist2D(
    " all jet calorimeter info", "calorimeter rapidity vs phi angle",
    "rapidity", "azimuthal angle phi",
    100, -6,6,
    100,-TMath::Pi(),TMath::Pi());

 plots->Jetimage1 = result->AddHist2D(
    "jet Image 1", "jet image 1",
    "rapidity", "azimuthal angle phi",
     50, -6,6,
    50,-4,4);
  plots->Jetimage2 = result->AddHist2D(
     "jet Image 2", "jet image 2",
    "rapidity", "azimuthal angle phi",
    50, -6,6,
    50,-4,4);
    
  plots->Jetimage3 = result->AddHist2D(
     "jet Image 3", "jet image 3",
    "rapidity", "azimuthal angle phi",
    50, -6,6,
    50,-4,4);


/*
  plots->Jetimage1 = result->AddHist2D(
    "jet Image 1", "jet image 1",
    "rapidity", "azimuthal angle phi",
    40, -0.8,0.8,
    40,-TMath::PiOver4(),TMath::PiOver4());

  plots->Jetimage2 = result->AddHist2D(
     "jet Image 2", "jet image 2",
    "rapidity", "azimuthal angle phi",
    40, -0.8,0.8,
    40,-TMath::PiOver4(),TMath::PiOver4());
    */
}
//-----------------------------------------------------------------------

void AnalyseEvents(ExRootTreeReader *treeReader, TestPlots *plots)
{
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
  TClonesArray *branchEFlowPhoton = treeReader->UseBranch("EFlowPhoton");
  TClonesArray *branchEFlowNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
    Long64_t allEntries = treeReader->GetEntries();


  cout << "** Chain contains " << allEntries << " events" << endl;
  GenParticle *particle;
  Electron *electron;
  Photon *photon;
  Muon *muon;
  Track *track;
  Tower *tower,*tower2;
  Jet *jet;
  TObject *object;
  TLorentzVector momentum;
  Float_t Eem, Ehad;
  Bool_t skip;
  Long64_t entry;
  Int_t i,n, j, pdgCode,hardest1,hardest2,hardest3;
  Float_t psuedo_rapidity, phi_angle, energy;
  Float_t rotation_angle, distance, phi_2, rap_2, rotation2;
  // Loop over all events
  //branch entry to plot the individual pictures for debugging
  int entrynumber=7359;

  for(entry = 0; entry < allEntries; ++entry)
  {
    // Load selected branches with data from specified event  good entries 9345,
    treeReader->ReadEntry(entry);
    
    if(branchJet->GetEntriesFast() >= 1)
    {
      jet = (Jet*) branchJet->At(0);
    
     Float_t e1=-1.0,e2=-1.0,e3=-1.0;
     for(n = 0; n < jet->Constituents.GetEntriesFast(); ++n)
        {
            object = jet->Constituents.At(n);
            if(object == 0) continue;
            if(object->IsA()== Tower::Class()) 
            {
               tower = (Tower*) object;
                energy=tower->E;
                if(energy>e1){
                    e1=energy;
                    hardest1=n;
                }
                else if(energy>e2){
                    e2=energy;
                    hardest2=n;
                }
                else if(energy>e3){
                e3=energy;
                hardest3=n;
                }
            }
            
        }
    
        // Loop over all jet's constituents and plot them
        for(j = 0; j < jet->Constituents.GetEntriesFast(); ++j)
        {
            object = jet->Constituents.At(j);
            
            if(object == 0) continue;
            
            
            if(object->IsA() == Tower::Class())
            {
                tower = (Tower*) object;
                //get each jet constiuents locations
                momentum += tower->P4();
                psuedo_rapidity = tower->Eta;
                phi_angle = tower->Phi;
                energy =tower->E;
                plots->fJet->Fill(psuedo_rapidity,phi_angle,energy);
                //plot single jet images:
                
                if(entry==entrynumber){plots->Jetimage1->Fill(psuedo_rapidity,phi_angle,energy);}
                
                //define location of second hardest particle for rotation
                tower=jet->Constituents.At(hardest2);
                phi_2 = tower->Phi;
                rap_2 =tower->Eta;
                
                //Translation: shift hardest to center 
                tower=jet->Constituents.At(hardest1);
                psuedo_rapidity=psuedo_rapidity-(tower->Eta);
                phi_angle=phi_angle-(tower->Phi);
                if(entry==entrynumber){plots->Jetimage2->Fill(psuedo_rapidity,phi_angle,energy);}
               
               //shift second hardest particle
                phi_2 = phi_2-(tower->Phi);
                rap_2 = rap_2-(tower->Eta);

                //define modulus and rotation angle
                rotation_angle=TMath::ATan2(phi_angle,psuedo_rapidity);
                distance=TMath::Sqrt(phi_angle*phi_angle + psuedo_rapidity*psuedo_rapidity);
                //rotate 2nd hardest particle to pi/2 or 12 oclock
                rotation2=TMath::PiOver2()-TMath::ATan2(phi_2,rap_2);
                
                //add the rotation for all particles
                rotation_angle = rotation2 + rotation_angle;
                psuedo_rapidity=distance*TMath::Cos(rotation_angle);
                phi_angle=distance*TMath::Sin(rotation_angle);
                
                //plot specific graph after rotation and shift
                if(entry==entrynumber){plots->Jetimage3->Fill(psuedo_rapidity,phi_angle,energy);}
                
            }
            
        }
        
     
    }
  }
}

//------------------------------------------------------------------------------
void PrintHistograms(ExRootResult *result, TestPlots *plots)
{
  //result->Print("png");
  TCanvas *c1 = new TCanvas();
  plots->fJet->Draw("colz");
  TCanvas *c2 = new TCanvas();
  plots->Jetimage1->Draw("colz");
  TCanvas *c3 = new TCanvas();
  plots->Jetimage2->Draw("colz");
  TCanvas *c4 = new TCanvas();
  plots->Jetimage3->Draw("colz");
}
//------------------------------------------------------------------------------
void calorimeter_drawer(const char *inputFile)
{
  gSystem->Load("libDelphes");
  TChain *chain = new TChain("Delphes");
  chain->Add(inputFile);
  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
  ExRootResult *result = new ExRootResult();
  TestPlots *plots = new TestPlots;
  BookHistograms(result, plots);
  AnalyseEvents(treeReader, plots);
  PrintHistograms(result, plots);
  result->Write("results.root");
  cout << "** Exiting..." << endl;
  /*
  delete plots;
  delete result;
  delete treeReader;
  delete chain;
  */
}
//------------------------------------------------------------------------------