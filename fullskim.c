/* call macro:
root -l mymacros/histogram2d.C'("delphes_output.root")'
*/
#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#include <string.h>
#else
#include <stdio.h>
#include <string.h>
class ExRootTreeReader;
class ExRootResult;

#endif
//------------------------------------------------------------------------------
struct TestPlots
{

    TH2 *Jetimage1;
    TH2 *Jetimage2;
};
//------------------------------------------------------------------------------
void BookHistograms(ExRootResult *result, TestPlots *plots)
{
    TLegend *legend;
    TPaveText *comment;

    plots->Jetimage1 = result->AddHist2D(
        "jet Image 1", "jet image 1",
        "rapidity", "azimuthal angle phi",
        40, -0.8, 0.8,
        40, -TMath::PiOver4(), TMath::PiOver4());

    plots->Jetimage2 = result->AddHist2D(
        "jet Image 2", "jet image 2",
        "rapidity", "azimuthal angle phi",
        40, -0.8, 0.8,
        40, -TMath::PiOver4(), TMath::PiOver4());
}
//-----------------------------------------------------------------------

int AnalyseEvents(ExRootTreeReader *treeReader, TestPlots *plots, int entrynumber)
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
    Track *track, *track2;
    Tower *tower, *tower2;
    Jet *jet;
    TObject *object, *daughter_object, *daughter_object2;
    TLorentzVector momentum, jet_momentum, top_momentum, *print_momentum;
    Float_t Eem, Ehad;
    Bool_t skip;

    Long64_t entry;
    entry = (Long64_t *)entrynumber;

    Int_t i, n, w, j, pdgCode, hardest1, hardest2, hardest3;
    Float_t psuedo_rapidity, phi_angle, energy;
    Float_t rotation_angle, distance, phi_2, rap_2, rotation2, phi_1, rap_1;
    Double_t deltaR, deltaR1, deltaR2;
    // Loop over all events
    //branch entry to plot the individual pictures for debugging
    int Drawplots = 0;

    // for(entry = 0; entry < allEntries; ++entry)
    //{

    treeReader->ReadEntry(entry);

    if (branchJet->GetEntriesFast() >= 1) //&& entry==entrynumber)
    {
        cout << "came atleast1once" << endl;
        jet = (Jet *)branchJet->At(0); //fix this
        cout << jet->Eta << " the jets eta, " << jet->Phi << " the jets phi, "
             << "event #" << entry << endl;
        if (jet->PT >= 550 && jet->PT <= 650 && TMath::Abs(jet->Eta) < 2) //PT550-650,ETA<2
        {

            jet_momentum = jet->P4();
            //check for top quark + within R=0.8
            for (w = 0; w < branchParticle->GetEntriesFast(); ++w)
            {
                GenParticle *gen = (GenParticle *)branchParticle->At(w);

                if (gen->PID == 6 || gen->PID == -6)
                {
                    if (gen->Status == 22)
                    {

                        top_momentum = gen->P4();
                        deltaR = top_momentum.DeltaR(jet_momentum);
                        cout << ", Status: " << gen->Status << ", PID: " << gen->PID << " event #" << entry << "  computed deltaR: " << deltaR << endl;

                        if (deltaR <= 0.8) //check daughters
                        {

                            daughter_object = branchParticle->At(gen->D1);

                            if (daughter_object->IsA() == GenParticle::Class())
                            {
                                particle = (GenParticle *)daughter_object;
                                deltaR1 = particle->P4().DeltaR(jet_momentum);
                            }
                            else if (daughter_object->IsA() == Track::Class())
                            {
                                track = (Track *)daughter_object;
                                deltaR1 = track->P4().DeltaR(jet_momentum);
                            }
                            else if (daughter_object->IsA() == Tower::Class())
                            {
                                tower = (Tower *)daughter_object;
                                deltaR1 = tower->P4().DeltaR(jet_momentum);
                            }

                            if (deltaR1 <= 0.8)
                            {

                                daughter_object2 = branchParticle->At(gen->D2);
                                if (daughter_object2->IsA() == GenParticle::Class())
                                {
                                    particle2 = (GenParticle *)daughter_object2;
                                    deltaR2 = particle2->P4().DeltaR(jet_momentum);
                                }
                                else if (daughter_object2->IsA() == Track::Class())
                                {
                                    track2 = (Track *)daughter_object2;
                                    deltaR2 = track2->P4().DeltaR(jet_momentum);
                                }
                                else if (daughter_object2->IsA() == Tower::Class())
                                {
                                    tower2 = (Tower *)daughter_object2;
                                    deltaR2 = tower2->P4().DeltaR(jet_momentum);
                                }
                                if (deltaR2 <= 0.8)
                                {
                                    Drawplots += 1;
                                }
                            }
                        }
                    }
                }
            }
        }

        if (Drawplots > 0)
        {
            Float_t e1 = -1.0, e2 = -1.0, e3 = -1.0;
            for (n = 0; n < jet->Constituents.GetEntriesFast(); ++n)
            {
                object = jet->Constituents.At(n);
                if (object == 0)
                    continue;

                if (object->IsA() == Tower::Class())
                {
                    tower = (Tower *)object;
                    energy = tower->E;
                }

                if (object->IsA() == Track::Class())
                {
                    track = (Track *)object;
                    energy = track->P4().Energy();
                    //nergy=track->E;
                }

                if (object->IsA() == GenParticle::Class())
                {
                    particle = (GenParticle *)object;
                    energy = particle->E;
                }
                if (energy > e1)
                {
                    e1 = energy;
                    hardest1 = n;
                }
                else if (energy > e2)
                {
                    e2 = energy;
                    hardest2 = n;
                }
                else if (energy > e3)
                {
                    e3 = energy;
                    hardest3 = n;
                }
            }
        }

        if (Drawplots > 0)
        {

            char directory3[100];
            snprintf(directory3, sizeof(directory3), "/data/data065/thoyte/out/ML/jet3/%dzzzz.txt", entrynumber);
            ofstream myfile;
            myfile.open(directory3);
            int counter;
            int z;

            for (j = 0; j < jet->Constituents.GetEntriesFast(); ++j)
            {
                if (j >= 200)
                { //get rid of non 200 constiuents
                    cout << "more then 200 constituents in jet " << endl;
                }

                else
                {

                    object = jet->Constituents.At(j);

                    if (object == 0)
                        continue;

                    //get each jet constiuents locations
                    if (object->IsA() == Tower::Class())
                    {
                        tower = (Tower *)object;
                        
                        print_momentum = (TLorentzVector *)tower->P4();
                        psuedo_rapidity = tower->Eta;
                        phi_angle = tower->Phi;
                        energy = tower->E;
                        counter += 1;
                    }

                    else if (object->IsA() == Track::Class())
                    {
                        track = (Track *)object;
                        print_momentum = (TLorentzVector *)track->P4();
                        psuedo_rapidity = track->Eta;
                        phi_angle = track->Phi;
                        energy = track->P4().Energy();
                        counter += 1;
                    std::myfile << " v Track particle v" << std::endl;

                    }

                
                    //put momentum in txt file
                    std::myfile << "( " << print_momentum->Energy() << " , " << print_momentum->Px() << " , " << print_momentum->Py() << " , " << print_momentum->Pz() << " )" << std::endl;

                    //define location of second hardest particle for rotation
                    object = (TObject *)jet->Constituents.At(hardest2);

                    if (object->IsA() == Tower::Class())
                    {
                        tower = (Tower *)object;
                        phi_2 = tower->Phi;
                        rap_2 = tower->Eta;
                    }
                    if (object->IsA() == GenParticle::Class())
                    {
                        particle = (GenParticle *)object;
                        phi_2 = particle->Phi;
                        rap_2 = particle->Eta;
                    }

                    //Translation: shift hardest to center
                    object = (TObject *)jet->Constituents.At(hardest1);
                    if (object->IsA() == Track::Class())
                    {
                        track = (Track *)object;
                        phi_1 = track->Phi;
                        rap_1 = track->Eta;
                    }
                    if (object->IsA() == Tower::Class())
                    {
                        tower = (Tower *)object;
                        phi_1 = tower->Phi;
                        rap_1 = tower->Eta;
                    }
                    if (object->IsA() == GenParticle::Class())
                    {
                        particle = (GenParticle *)object;
                        phi_1 = particle->Phi;
                        rap_1 = particle->Eta;
                    }

                    psuedo_rapidity = psuedo_rapidity - rap_1;
                    phi_angle = phi_angle - phi_1;
                    plots->Jetimage1->Fill(psuedo_rapidity, phi_angle, energy);

                    //shift second hardest particle
                    phi_2 = phi_2 - phi_1;
                    rap_2 = rap_2 - rap_1;

                    //define modulus and rotation angle
                    rotation_angle = TMath::ATan2(phi_angle, psuedo_rapidity);
                    distance = TMath::Sqrt(phi_angle * phi_angle + psuedo_rapidity * psuedo_rapidity);
                    //rotate 2nd hardest particle to pi/2 or 12 oclock
                    rotation2 = TMath::PiOver2() - TMath::ATan2(phi_2, rap_2);

                    //add the rotation for all particles
                    rotation_angle = rotation2 + rotation_angle;
                    psuedo_rapidity = distance * TMath::Cos(rotation_angle);
                    phi_angle = distance * TMath::Sin(rotation_angle);

                    //plot specific graph after rotation and shift
                    plots->Jetimage2->Fill(psuedo_rapidity, phi_angle, energy);
                }
            }

            if(counter<200) //add zeros 
            {
                for(z=0; z < (200-counter) ; z++)
                {
                   std::myfile << "( 0 ,  0 , 0 , 0 )" <<std::endl;
                }

            }

            myfile.close();
        }
        //closing text file
    }

    return Drawplots;
}

//------------------------------------------------------------------------------
void PrintHistograms(ExRootResult *result, TestPlots *plots, int entrynumber, int Drawplots)
{
    if (Drawplots == 1)
    {

        //result->Print("png");
        TCanvas *c1 = new TCanvas();
        plots->Jetimage1->Draw("colz");
        TCanvas *c2 = new TCanvas();
        plots->Jetimage2->Draw("colz");

        char directory1[100];
        char directory2[100];

        snprintf(directory1, sizeof(directory1), "/data/data065/thoyte/out/ML/jet1/%dzzzz.jpg", entrynumber);
        snprintf(directory2, sizeof(directory2), "/data/data065/thoyte/out/ML/jet2/%dzzzz.jpg", entrynumber);

        c1->SaveAs(directory1);
        c2->SaveAs(directory2);

        delete c1;
        delete c2;
    }
    else
    {
        cout << "nothing to save" << endl;
    }
}
//------------------------------------------------------------------------------
void fullskim(const char *inputFile)
{
    gSystem->Load("libDelphes");
    TChain *chain = new TChain("Delphes");
    chain->Add(inputFile);

    int entrynumber; // =9298;
    int maxentry = 50000;
    int Drawplots;
    for (entrynumber = 0; entrynumber < maxentry; entrynumber++)
    {
        ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
        ExRootResult *result = new ExRootResult();
        TestPlots *plots = new TestPlots;
        BookHistograms(result, plots);
        Drawplots = AnalyseEvents(treeReader, plots, entrynumber);
        PrintHistograms(result, plots, entrynumber, Drawplots);
        delete plots;
        delete result;
        delete treeReader;
        //try deleteing canvases every time
    }
    //result->Write("results.root");
    cout << "** Exiting..." << endl;
    /*
  
  delete chain;
  */
}
//------------------------------------------------------------------------------