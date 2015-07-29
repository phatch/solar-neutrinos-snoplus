/* Pi Monte Carlo Program 2015 */

#include <cmath>

using namespace std;


int pi_monte_carlo(){
    TRandom1 rn(643278);
    const int tot = 300000;
    int count = 0;
    TCanvas *c = new TCanvas();
    double pi[tot],trial[tot];
    
    for(int i = 1; i <= tot; i++){
        double rn1 = rn.Rndm();
        double rn2 = rn.Rndm();
        if (sqrt(pow(rn1, 2.0) + pow(rn2, 2.0)) < 1)
            count = ++count;
        pi[i-1] = 4*(static_cast<double>(count)/(i));
        trial[i-1] = i;
    }
    
    TGraph *g = new TGraph(tot,trial,pi);
    g->Draw("Ap");
    g->SetTitle("$\\pi$ Monte\\ Carlo");
    g->GetXaxis()->SetTitle("Number of trials");
    g->GetYaxis()->SetTitle("$\\pi$ approximation");
    g->GetXaxis()->CenterTitle();
    g->GetYaxis()->CenterTitle();
    cout << "Pi is approximately " << 4*(static_cast<double>(count)/static_cast<double>(tot));
    return 0;
    
}
