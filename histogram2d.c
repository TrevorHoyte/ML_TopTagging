
/* call macro:
root -l mymacros/histogram2d.C'("delphes_output.root")'
*/
#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
//#include "TMath.h
#include <stdbool.h>
#else
class ExRootTreeReader;
class ExRootResult;
#endif
//------------------------------------------------------------------------------
struct TestPlots
{
 
  TH2 *fJet;
  TH1 *ftops;

};
//------------------------------------------------------------------------------
void BookHistograms(ExRootResult *result, TestPlots *plots)
{
  TLegend *legend;
  TPaveText *comment;
 
  plots->fJet = result->AddHist2D(
    "jet calorimeter info", "calorimeter rapidity vs phi angle",
    "rapidity", "azimuthal angle phi",
    1000, -100,100,
    1000,0,360);
	plots->jtops = result->AddHist1D("status of top","topstatus","stat",100,-20,20)
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
  GenParticle *particle, *particle2;
  Electron *electron;
  Photon *photon;
  Muon *muon;
  Track *track, *track2;
  Tower *tower, *tower2;
  Jet *jet;
  TObject *daughter_object, *daughter_object2;
  TLorentzVector momentum;
  TLorentzVector jet_momentum , top_momentum ,daughter_momentum;
  //TLorentzVector *;
  Float_t Eem, Ehad;
  Bool_t skip;
  Long64_t entry;
  Int_t i,w, n, j, pdgCode;
  Double_t  deltaR,deltaR1, deltaR2 ;
  bool Jetisgood;
  
  // Loop over all events
  for(entry = 0; entry < allEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    
    // load all jets in event
    for(i = 0; i < branchJet->GetEntriesFast(); ++i)
    {
      jet = (Jet*) branchJet->At(i);
      Jetisgood=false;
      
      
      //only allow 550-650 GEV jets
      if(jet->PT>=550 && TMath::Abs(jet->Eta)<2)
      {
        //cout<<jet->Eta<<endl;
        
        //find jet 4vector
        jet_momentum = jet->P4();
        cout << jet->Eta << " the jets eta, "<<jet->Phi<<" the jets phi, " << "event #"<< entry << "jet number "<<i<<endl; 
        //check for top quark + all specifications
        for(w = 0; w < branchParticle->GetEntriesFast(); ++w)
        {
            GenParticle *gen = (GenParticle*) branchParticle->At(w);
            
            if(gen->PID==6 || gen->PID==-6)
            {
                plots->ftops->Fill(gen-Status());
                //find top 4 vector
                top_momentum = gen->P4();
                deltaR = top_momentum.DeltaR(jet_momentum);
                //cout<<", St: "<<gen->Status<<", PID: "<<gen->PID<< " event #"<< entry <<"  computed deltaR: "<<deltaR<< endl;
                 
                //check daughters
                if(deltaR<=0.8)
                {
                    daughter_object = branchParticle->At(gen->D1);
                    
                    if(daughter_object->IsA()== GenParticle::Class()){
                       particle = (GenParticle*) daughter_object; 
                       deltaR1 = particle->P4().DeltaR(jet_momentum);
                    }
                    else if(daughter_object->IsA()== Track::Class()){
                       track = (Track*) daughter_object; 
                       deltaR1 = track->P4().DeltaR(jet_momentum);
                    }
                    else if(daughter_object->IsA()== Tower::Class()){
                       tower = (Tower*) daughter_object; 
                       deltaR1 = tower->P4().DeltaR(jet_momentum);
                    }
                    
                    if(deltaR1<=0.8)
                    {
                        daughter_object2 = branchParticle->At(gen->D2);
                        if(daughter_object2->IsA()== GenParticle::Class()){
                        particle2 = (GenParticle*) daughter_object2; 
                        deltaR2 = particle2->P4().DeltaR(jet_momentum);
                        }
                        else if(daughter_object2->IsA()== Track::Class()){
                        track2 = (Track*) daughter_object2; 
                        deltaR2 = track2->P4().DeltaR(jet_momentum);
                        }
                        else if(daughter_object2->IsA()== Tower::Class()){
                        tower2 = (Tower*) daughter_object2; 
                        deltaR2 = tower2->P4().DeltaR(jet_momentum);
                        }
                        if(deltaR2<=0.8){Jetisgood=true;}
                        //cout<<"deltaR top: "<<deltaR<<" deltaR 1st daughter: "<<deltaR1<<" deltaR 2nd daughter: "<<deltaR2<<endl; 
                    }
                    
                    
                }
            } 

        }

       

      }

      //play with good jets
      if(Jetisgood)
      {
        for(j = 0; j < jet->Constituents.GetEntriesFast(); ++j)
        {
            object = jet->Constituents.At(j);
            // Check if the constituent is accessible
            if(object == 0) continue;
           
             if(object->IsA() == Tower::Class())
            {
            tower = (Tower*) object;
            cout<<"tower Eta: "<<tower->Eta<<"  tower phi: "<<tower->Phi<<endl;
            }
        }

      }  










    }

      
  }
  
}

//------------------------------------------------------------------------------
void PrintHistograms(ExRootResult *result, TestPlots *plots)
{
  result->Print("png");
}
//------------------------------------------------------------------------------
void histogram2d(const char *inputFile)
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