//
//  Expected_B8_neutrino_events.C
//  
//
//  Created by Patrick Hatch on 2015-07-01.
//
//

#include <iostream>
#include <stdio.h>
#include <cmath>
#include <vector>
#include <sstream>
#include <fstream>
#include <string>
#include "TGraph.h"

// The purpose of this code is simply to calculate the number of B8 (boron 8) neutrino we expect to see in a water
// detector (specifically SNO+). It requires ROOT to run.

using namespace std;

vector<double> q_vec, P_q;
double q_arr[150], lambda_q[150], x[1500], y[1500];

const double sin2Wein = 0.23116; // The squared sine of the Weinberg angle
const double gL = 0.5 + sin2Wein; // See Bahcall (1989) p 217 for details on these constants
const double gR = sin2Wein;
const double gL2 = gL*gL;
const double gR2 = gR*gR;
const double eXSec = 88.083e-46; // in cm^2
const double mec2 = 0.511; // MeV (The rest energy of an electron)
const double B8Flux = 5.8*1e6; // cm^-2 s^-1 (The flux of B8 neutrinos)
const double fullWaterWeight = 0.905; //kton (The amount of water in the fully filled SNO+ detector)
const double partialFilledWeight = 0.185*fullWaterWeight; //kton (The fraction of water in our partial filled run)
const double lengthOfRun = (7.0/365.0); //years (this value is roughtly a week)

double TMax(double TMaxSet, double TMinSet, double q)
{
    double TMax = (2.0*q*q)/(mec2*(1.0 + (2.0*q)/(mec2))); // See Bahcall (1989) p 217 eq. 8.32, converted from units of mec2
    if(TMax > TMaxSet) // TMax can't go above our upper energy limit
        TMax = TMaxSet;
    if(TMax < TMinSet) // TMax can't be below TMin, thus this would set the eXSecTotQ to zero
        TMax = TMinSet;
    return TMax;
}

double eXSecTotQ(double TMaxSet, double TMinSet, double q)
{
    // See Bahcall (1989) p 217 eq. 8.33, again converted from units of mec2
    double eXSecTotQ = (eXSec/mec2)*(((gL2 + gR2)*(TMax(TMaxSet,TMinSet,q) - TMinSet))
                                    - (((gR2/q) + (gL*gR*mec2)/(2.0*q*q))*(TMax(TMaxSet,TMinSet,q)*TMax(TMaxSet,TMinSet,q) - TMinSet*TMinSet))
                                    + ((gR2/(3.0*q*q))*(TMax(TMaxSet,TMinSet,q)*TMax(TMaxSet,TMinSet,q)*TMax(TMaxSet,TMinSet,q) - TMinSet*TMinSet*TMinSet)));
    return eXSecTotQ;
}

void Expected_B8_neutrino_events()
{
    vector<double> q_vec, P_q;
    double q_arr[150], lambda_q[150], x[1500], y[1500];
    
    // Here we store the numerical table given on p 156 of Bahcall (1989), which is the normalized B8 neutrino spectrum, as 2 vectors
    string probFunc = "/Users/patrickhatch/Desktop/Summer_2015/B8_spec_v2.txt"; // This file is a text file of the above mentioned table
    ifstream mfile;
    mfile.open(probFunc.data());
    
    if (!mfile.is_open())
    {
        cerr << "File cannot be opened!" << endl;
        return;
    }
    
    string dummy;
    double mdummy;
    
    while(getline(mfile, dummy))
    {
        istringstream is(dummy);
        int colcount = 0;
        
        while(is >> mdummy)
        {
            if(colcount == 0)
                q_vec.push_back(mdummy);
            
            if(colcount == 1)
                P_q.push_back(mdummy);
            
            colcount++;
            
            if(colcount == 2)
                colcount = 0;
        }
    }
    
    for(unsigned int i = 0; i < q_vec.size(); i++)
    {
        lambda_q[i] = P_q[i];
        q_arr[i] = q_vec[i];
    }
    
    
    // Here we just set up a simple numerical integral to perform the integration indicated in Bahcall (1989) p 218 eq. 8.34
    const double interval = 0.00001;
    const double interval_length = 15.0;
    const double stepNum = interval_length/interval;
    const double deltaQ = interval_length/stepNum;
    
    double sum = 0;
    double x = 0;
    
    // Here we just fit a ninth order polynominal to the normalized B8 neutrino spectrum
    TGraph *graph = new TGraph(q_vec.size(),q_arr,lambda_q);
    graph->Fit("pol9");
    
    const double OurTMin = 4.0; //MeV
    const double OurTMax = 10.0; //MeV
    
    for(int k = 0; k < static_cast<int>(stepNum); k++)
    {
        x = k*interval;
        if(x != 0)
            sum += (graph->GetFunction("pol9")->Eval(x))*eXSecTotQ(OurTMax,OurTMin,x)*deltaQ;
    }
    
    // See Bahcall (1989) p 382 eq. 13.2 for details
    double numOfEvents = 610*(sum/(1e-44))*(B8Flux/(5.8*1e6))*partialFilledWeight*lengthOfRun;
    
    cout << "The number of expected events is " << numOfEvents << endl;
    
}
// References: Bahcall, J. (1989). Neutrino astrophysics. Cambridge: Cambridge University Press.
