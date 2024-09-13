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

//----
int findDaughters1(TClonesArray *branchParticle, Long64_t entry, Int_t particle, int daughternumber)
{

    if (particle == -1)
    {
        //cout << "particle is -1" << endl;
        return -1;
    }
    else
    {
        TObject *obj = (TObject *)branchParticle->At(particle);

        if (obj->IsA() == GenParticle::Class())
        {
            GenParticle *gen = (GenParticle *)obj;
            //cout << daughternumber << " :Daughter number, particle numb: "<<particle<<" Status: " << gen->Status << ", PID: " << gen->PID <<"  "<< gen->M1<< ","<<gen->M2<<" :mother M1, M2  ,entry number :" << entry << endl;
            if (gen->PID == 24 || gen->PID == -24)
            {

                //cout << "found W boson" << endl;
                return particle;
            }
            else if (gen->PID == 5 || gen->PID == -5)
            {
                //cout << "found bottom quark" << endl;
                return particle;
            }

            else
            {
                return findDaughters1(branchParticle, entry, gen->D1, daughternumber + 1);
            }
        }
    }
}
//----------------------------------------------------------------------------------------------------------------

int findDaughters2(TClonesArray *branchParticle, Long64_t entry, Int_t particle, int daughternumber)
{

    if (particle == -1)
    {
        cout << "particle is -1" << endl;
        return -1;
    }
    else
    {
        TObject *obj = (TObject *)branchParticle->At(particle);

        if (obj->IsA() == GenParticle::Class())
        {
            GenParticle *gen = (GenParticle *)obj;
            //cout << daughternumber << " :Daughter number, particle numb: "<<particle<<" Status: " << gen->Status << ", PID: " << gen->PID <<"  "<< gen->M1<< ","<<gen->M2<<" :mother M1, M2  ,entry number :" << entry << endl;
            if (gen->PID == 24 || gen->PID == -24)
            {
                //cout << "found W boson" << endl;
                return particle;
            }
            else if (gen->PID == 5 || gen->PID == -5)
            {
                int part = 69;
                //cout << "found bottom quark" << endl;

                return particle;
            }

            else
            {
                return findDaughters2(branchParticle, entry, gen->D2, daughternumber + 1);
            }
        }
    }
}

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

void AnalyseEvents(ExRootTreeReader *treeReader)
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
    GenParticle *particle, *daughter_wboson, *daughter_bottom;
    Electron *electron;
    Photon *photon;
    Muon *muon;
    Track *track, *track2;
    Tower *tower, *tower2;
    Jet *jet;
    TObject *object, *daughter_object, *daughter_object2;
    TLorentzVector momentum, jet_momentum, top_momentum, *print_momentum, wboson_momentum, bottom_momentum;
    Float_t Eem, Ehad;
    Bool_t skip;

    Long64_t entry;
    

    Int_t i, n, w, j, pdgCode, wboson, bottom;
    Float_t psuedo_rapidity, phi_angle, energy;
    Float_t rotation_angle, distance, phi_2, rap_2, rotation2, phi_1, rap_1;
    Double_t deltaR, deltaR_bottom, deltaR_wboson;

    // Loop over all events
    //branch entry to plot the individual pictures for debugging

    for (entry = 2397; entry < allEntries; ++entry)
    {
        TestPlots *plots = new TestPlots;
        ExRootResult *result = new ExRootResult();
        BookHistograms(result, plots);
        int Drawplots = 0;

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

                                wboson = (Int_t *)findDaughters1(branchParticle, entry, gen->D1, 1);
                                bottom = (Int_t *)findDaughters2(branchParticle, entry, gen->D1, 1);
                                if (wboson != -1 && bottom != -1)
                                {

                                    daughter_wboson = (GenParticle *)branchParticle->At(wboson);
                                    daughter_bottom = (GenParticle *)branchParticle->At(bottom);

                                    wboson_momentum = daughter_wboson->P4();
                                    bottom_momentum = daughter_bottom->P4();
                                    ///compute deltaR
                                    deltaR_bottom = bottom_momentum.DeltaR(jet_momentum);
                                    deltaR_wboson = wboson_momentum.DeltaR(jet_momentum);
                                    cout << "Wboson and bottoms R delta R " << deltaR_wboson << " , " << deltaR_bottom << endl;
                                    if (deltaR_bottom <= 0.8 && deltaR_wboson <= 0.8)
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
                Int_t hardest1 = -1, hardest2 = -1, hardest3 = -1;
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
                        e3 = e2;
                        hardest3 = hardest2;
                        e2 = e1;
                        hardest2 = hardest1;

                        e1 = energy;
                        hardest1 = n;
                    }
                    else if (energy > e2)
                    {
                        e3 = e2;
                        hardest3 = hardest2;

                        e2 = energy;
                        hardest2 = n;
                    }
                    else if (energy > e3)
                    {
                        e3 = energy;
                        hardest3 = n;
                    }
                }
                cout << "hardest entry number:energy  h1: " << hardest1 << " : " << e1 << "  , h2:  " << hardest2 << " : " << e2 << " ,h3  " << hardest3 << " : " << e3 << endl;
            }

            if (Drawplots > 0)
            {
                cout<<" 1 debug "<<endl;
                char directory3[100];
                cout<<" 2 debug "<<endl;
                snprintf(directory3, sizeof(directory3), "/data/data065/thoyte/out/ML/test2/ttbar%dxxxxx.txt", entry);
                ofstream myfile;
                cout<<" 3 debug "<<endl;
                myfile.open(directory3);
                int counter;
                int foundobj;
                int z;
                //TLorentzVector momentum_sum;
                //momentum_sum.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);

                for (j = 0; j < jet->Constituents.GetEntriesFast(); ++j)
                {
                    if (j >= 200)
                    { //get rid of non 200 constiuents
                        cout << "more then 200 constituents in jet " << endl;
                    }

                    else
                    {
                        foundobj = 0;

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
                            foundobj += 1;
                            //momentum_sum+=tower->P4();
                            std::myfile << "( " << print_momentum->Energy() << " , " << print_momentum->Px() << " , " << print_momentum->Py() << " , " << print_momentum->Pz() << " )" << std::endl;
                        }

                        else if (object->IsA() == Track::Class())
                        {
                            track = (Track *)object;
                            print_momentum = (TLorentzVector *)track->P4();
                            psuedo_rapidity = track->Eta;
                            phi_angle = track->Phi;
                            energy = track->P4().Energy();
                            counter += 1;
                            foundobj += 1;
                            //momentum_sum+=track->P4();
                            //std::myfile << " v Track particle v" << std::endl;
                            std::myfile << "( " << print_momentum->Energy() << " , " << print_momentum->Px() << " , " << print_momentum->Py() << " , " << print_momentum->Pz() << " )" << std::endl;
                        }

                        
                        if (foundobj > 0) // plot calorimeter objects
                        {

                            //define location of second hardest particle for rotation
                            object = (TObject *)jet->Constituents.At(hardest2);

                            if (object->IsA() == Tower::Class())
                            {
                                tower = (Tower *)object;
                                phi_2 = tower->Phi;
                                rap_2 = tower->Eta;
                            }
                            if (object->IsA() == Track::Class())
                            {
                                track = (Track *)object;
                                phi_2 = track->Phi;
                                rap_2 = track->Eta;
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
                }
                PrintHistograms(result, plots, entry);

                if (counter < 200) //add zeros
                {
                    for (z = 0; z < (200 - counter); z++)
                    {
                        std::myfile << "( 0 ,  0 , 0 , 0 )" << std::endl;
                    }
                }
                //std::myfile << "( " << momentum_sum.Energy() << " , " << momentum_sum.Px() << " , " << momentum_sum.Py() << " , " << momentum_sum.Pz() << " )" << std::endl;
                //std::myfile<<"jet energy: "<<jet->P4()->Energy()<<"  jet PT: "<<jet->PT<<std::endl;
                myfile.close();
            }
            //closing text file
        }

        delete plots;
        delete result;
        delete gROOT->FindObject("Jetimage1");
        delete gROOT->FindObject("Jetimage2");
    }
}

//------------------------------------------------------------------------------
void PrintHistograms(ExRootResult *result, TestPlots *plots, Int_t entrynumber)
{

    //result->Print("png");
    TCanvas *c1 = new TCanvas();
    plots->Jetimage1->Draw("colz");
    TCanvas *c2 = new TCanvas();
    plots->Jetimage2->Draw("colz");

    char directory1[100];
    char directory2[100];

    //snprintf(directory1, sizeof(directory1), "/data/data065/thoyte/out/ML/jet2/%dzzzz.jpg", entrynumber);
    snprintf(directory2, sizeof(directory2), "/data/data065/thoyte/out/ML/test2/ttbar%dxxxxx.jpg", entrynumber);

    //c1->SaveAs(directory1);
    c2->SaveAs(directory2);

    //TCanvas *c1 = (TCanvas *)gROOT->GetListOfCanvases()->FindObject("c1");
    //TCanvas *c2 = (TCanvas *)gROOT->GetListOfCanvases()->FindObject("c2");
    delete c1;
    delete c2;
}
//------------------------------------------------------------------------------
void ttbarfix(const char *inputFile)
{
    gSystem->Load("libDelphes");
    TChain *chain = new TChain("Delphes");
    chain->Add(inputFile);
    ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
    AnalyseEvents(treeReader);
    delete treeReader;
    delete chain;
}
//------------------------------------------------------------------------------