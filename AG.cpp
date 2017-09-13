// Gillespie Algorithm, solving reation-diffusion equation.
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <vector>
#include <random>
using namespace std;



// Structure for input parameters
typedef struct initParameters {
    
    int nsample;
    int nstep;
    int nmonomer;
    float length;
    int compartment;
    int primnuc;
    float kp;
    float km;
    
    
    initParameters( int nsampleIn, int nmonomerIn, float lengthIn, int compartmentIn, int primnucIn, float kpIn, float kmIn, int nstepIn ) {
        this->nsample = nsampleIn;
        this->nstep = nstepIn;
        this->nmonomer = nmonomerIn;
        this->length = lengthIn;
        this->compartment = compartmentIn;
        this->primnuc = primnucIn;
        this->kp = kpIn;
        this->km = kmIn;
        
    }
    
    
}initParameters;


double diffconst(int i, int primnuc, float length, int compartment);


vector<double> propensity( vector<int> a0, vector<int> a, int i, float kp, float km, int primnuc, float length, int compartment, double (*d)(int , int , float, int) );



int main(){
    
    //================== Read initial parameters from parameters.ini =======================
    ifstream paraini("parameters.ini");
    string name;
    double para_ini[8+1];
    
    int ipara = 1;
    while( paraini >> name >> para_ini[ipara] ){
        ipara++;
    }
    
    const int nsample = para_ini[1];
    const int nstep = para_ini[2];
    const int nmonomer = para_ini[3];
    const float length = para_ini[4];
    const int compartment = para_ini[5];
    const int primnuc = para_ini[6];
    const float kp = para_ini[7];
    const float km = para_ini[8];
    
    
    //=======================================================================================
    
    
    cout << diffconst(1, primnuc, length, compartment) << "\n";
    
    
    vector<double> a(3,0);
    
    cout << a[0] << "\n";
    
    vector<double> c(2,0);
    //c = propensity(&a,1,kp,km);
    
    //cout << c[0] << "\t" << c[1] << "\n";
    
    
    
    //=============== initialize vectors for monomer and nucleation seed ===============
    vector<int> adatarow(compartment, 0);
    vector< vector<int> > adata(2, adatarow);
    
    double totalpropensity = 0;
    
    cout << adata.size() << endl;
    cout << adata.at(1).at(5) << endl;
    
    for(int i = 20; i < 30; i++){
        adata.at(0).at(i) = nmonomer/10;
    }
    
    adata.at(1).at(24) = 10;
    adata.at(1).at(25) = 30;
    //================================================================================
    
    //vector<double> j(4,0);
    //j = propensity(adata.at(0), adata.at(1), 1, kp, km, primnuc, length, compartment, diffconst);
    
    
    for(int istep=0; istep < nstep +1; istep++ ){
        
        unsigned seeds = std::chrono::system_clock::now().time_since_epoch().count();
        default_random_engine generator (seeds);
        uniform_real_distribution<float> random (0.000001,1);
        
        float r1, r2;
        
        r1 = random(generator);
        r2 = random(generator);
        
        cout << r1 << "\t" << r2 << endl;
        
        
        
        
        
        
        
    }
    
    
    
    
    
    
    
    

    
    
}


double diffconst(int i, int primnuc, float length, int compartment){

    int l;
    float h = length/compartment;
    double d;
    if(i == 0)
        l = 1;
    else
        l = i + primnuc - 1;
    
    double R;
    R = l * 10 * 0.001/ (2*log10(2*l*10/8));
    
    d = 0.483 / (h*h*R);
    
    return d;
}

vector<double> propensity( vector<int> a0, vector<int> a, int i, float kp, float km, int primnuc, float length, int compartment, double (*d)(int, int, float, int) ){
    
    double diffconst;
    diffconst = (*d)(i, primnuc, length, compartment);
    
    // ip = ( ipright, ipleft, ipdissolve, ipplusone )
    vector<double> ip(4,0);
    
    int a0sum = 0;
    int asum = 0;
    int a0asum = 0;
    
    for(int i = 0; i<a0.size();i++){
        a0sum = a0sum + a0.at(i);
        asum = asum + a.at(i);
        a0asum = a0asum + a0.at(i) * a.at(i);
    }
    
    
    ip.at(0) = diffconst * ( asum - a.back() );
    ip.at(1) = diffconst * ( asum - a.front() );
    
    if(i == 0 or i == 1 ){
        ip.at(2) = 0;
        ip.at(3) = kp * a0asum * i;
    }
    else if( 1 < i < primnuc+1 ){
        ip.at(2) = 2. * km * asum * 1;
        ip.at(3) = kp * a0asum;
    }
    else{ // i - (nc+1) + 1 + 2
        ip.at(2) = km * asum * (i - primnuc + 2 );
        ip.at(3) = kp * a0asum;
    }
    
    
    return ip;
}







