#include <RcppArmadillo.h>
#include <string.h>
#include <iostream> 
#include <sstream> 
#include <utility>
#include <unordered_map>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;
// [[Rcpp::plugins("cpp11")]]

void pairsort1(arma::imat & A, arma::ivec & b){
    arma::ivec  a = A.col(0);
    int         n = a.n_elem;
    int         m = b.n_elem;
    if(n != m){
        exit(EXIT_SUCCESS);
    }
    pair<int, int> pairt[n];
    for (int i = 0; i < n; i++){ 
        pairt[i].first  = a[i]; 
        pairt[i].second = b[i]; 
    } 
    
    sort(pairt, pairt + n); 
    
    // Modifying original arrays 
    for (int i = 0; i < n; i++){ 
        a[i] = pairt[i].first; 
        b[i] = pairt[i].second; 
    } 
}

void pairsort2(arma::uvec & a, arma::uvec & b){
    //arma::uvec  a = x; 
    int         n = a.n_elem;
    int         m = b.n_elem;
    if(n != m){
        exit(EXIT_SUCCESS);
    }
    pair<int, int> pairt[n];
    for (int i = 0; i < n; i++){ 
        pairt[i].first  = a[i]; 
        pairt[i].second = b[i]; 
    } 
    
    sort(pairt, pairt + n); 
    
    // Modifying original arrays 
    for (int i = 0; i < n; i++){ 
        //a[i] = pairt[i].first; 
        b[i] = pairt[i].second; 
    } 
}

// [[Rcpp::export]]
arma::ivec OrderPed(const StringVector & Anim, const StringVector & Sire, const StringVector & Dam){
    unordered_map <string, int> links; // Use unordered map for fast lookup
    int                         nanim   = Anim.size();
    int                         maxiter = 1000;
    int                         count   = 0;
    int                         i       = 1;
    string                      ks;
    string                      kd;
    int                         js;
    int                         jd;
    int                         gen;
    int                         delta;
    string                      aid;
    arma::ivec                  Old     = ones<ivec>(nanim, 1);
    arma::ivec                  New     = Old;
    wall_clock timer;
    /*
     * Initialize unordered map
     */
    timer.tic();
    for(int ia = 0; ia < nanim; ia++){
        aid = Anim[ia];
        links.insert({aid, ia});
    }
    /*
     * Main loop
     */
    while(i > 0){
        for(int j = 0; j < nanim; j++){
            ks  = Sire[j];
            kd  = Dam[j];
            gen = New[j] + 1;
            if(ks != "0"){
                js = links[ks];
                if(gen > New[js]){
                    New[js] = gen;
                }
            }
            if(kd != "0"){
                jd = links[kd];
                if(gen > New[jd]){
                    New[jd] = gen;
                }
            }
        }
        delta   = accu(New - Old);
        Old     = New;
        i       = delta;
        count   ++;
        if(count > maxiter){
            i = 0;
        }
    }
    double tt = timer.toc();
    cout << "Pedigree sorted in " << count << " iterations (" << tt << " seconds)" << endl;
    return(New);
}

// [[Rcpp::export]]
arma::umat RecodePed(const StringVector & ID, const StringVector & Sire, const StringVector & Dam){
    int                         nanim   = Sire.size();
    int                         m       = Dam.size();
    int                         j       = 0;
    int                         k       = 0;
    string                      aid;
    string                      sid;
    string                      did;
    unordered_map <string, int> links;
    arma::umat                  RecPed  = zeros<umat>(nanim, 3);
    wall_clock                  timer;
    if(nanim != m){
        exit(EXIT_SUCCESS);
    }
    /*
     * Initialize unordered_map 
     */
    timer.tic();
    for(int i = 0; i < nanim; i++){
        aid = ID[i];
        links.insert({aid, i + 1});
    }
    /*
     * Filling Recoded pedigree
     */
    for(int i = 0; i < nanim; i++){
        //aid             = ID[i];
        //links.insert({aid, i + 1});
        sid             = Sire[i];
        did             = Dam[i];
        j               = links[sid];
        k               = links[did];
        RecPed.at(i, 0) = i + 1;
        RecPed.at(i, 1) = j;
        RecPed.at(i, 2) = k;
    }
    double tt = timer.toc();
    cout << "Pedigree recoded. " << endl;
    cout << "Pedigree size:    " << nanim << " individuals" << endl;
    cout << "Time elapsed:     " << tt << " seconds." << endl;
    return(RecPed);
}

arma::vec InbreedingMC2(int n, int m, arma::imat Ped){
    int         i;
    int         j;
    int         k;
    int         rN;
    int         rS;
    int         S;
    int         D;
    int         MIP;
    
    arma::imat  rPed    = zeros<imat>(m + 1, 2);
    arma::ivec  SId     = zeros<ivec>(n + 1, 1);
    arma::ivec  Link    = zeros<ivec>(n + 1, 1);
    arma::ivec  MaxIdP  = zeros<ivec>(m + 1, 1);
    arma::vec   F       = zeros(n + 1);
    arma::vec   B       = zeros(m + 1);
    arma::vec   x       = zeros(m + 1);
    
    F[0]    = -1.0;
    x[0]    = 0.0;
    Link[0] = 0;
    
    for(rN = i = 1; i <= n; i++){
        SId[i]      = i;
        Link[i]     = 0;
        if(i <= m){
            x[i]    = 0.0;
        }
        S = Ped(i, 0);
        D = Ped(i, 1);
        if(S && !Link[S]){
            MaxIdP[rN]      = Link[S] = rN;
            rPed(rN, 0)     = Link[Ped(S, 0)];
            rPed(rN++, 1)   = Link[Ped(S, 1)];
        }
        if(D && !Link[D]){
            Link[D]         = rN;
            rPed(rN, 0)     = Link[Ped(D, 0)];
            rPed(rN++, 1)   = Link[Ped(D, 1)];
        }
        if(MaxIdP[Link[S]] < Link[D]){
            MaxIdP[Link[S]] = Link[D];
        }
    }
    // Sort animals according to ID of their sires into SId
    pairsort1(Ped, SId);
    for(k = i = 1; i <= n;){
        if(!Ped(SId[i], 0)){
            F[SId[i++]] = 0.0;
        } else {
            S       = Ped(SId[i], 0);
            rS      = Link[S];
            MIP     = MaxIdP[rS];
            x[rS]   = 1.0;
            for(; k <= S; k++){
                if(Link[k]){
                    B[Link[k]] = 0.5 - 0.25*(F[Ped(k, 0)] + F[Ped(k, 1)]);
                }
            }
            for(j = rS; j; --j){
                if(x[j]){
                    if(Ped(j, 0)){
                        x[rPed(j, 0)] += x[j]*0.5;
                    }
                    if(Ped(j, 1)){
                        x[rPed(j, 1)] += x[j]*0.5;
                    }
                    x[j] *= B[j];
                }
            }
            for(j = 1; j <= MIP; j++){
                x[j] += (x[rPed(j, 0)] + x[rPed(j, 1)])*0.5;
            }
            for(; i <= n; i++){
                if(S != Ped(SId[i], 0)){
                    break;
                } else {
                    F[SId[i]] = x[Link[Ped(SId[i], 1)]]*0.5;
                }
            }
            for(j = 1; j <= MIP; j++){
                x[j] = 0.0;
            }
        }
    }
    return(F);
}

// [[Rcpp::export]]
arma::vec InbreedingMC(arma::uvec siretmp, arma::uvec damtmp){
    int         n   = siretmp.n_rows;
    int         m   = 0; 
    int         i   = 0;
    int         j   = 0;
    int         k   = 0;
    int         rn  = 0;
    int         rs  = 0;
    int         S   = 0;
    int         D   = 0;
    int         MIP = 0;
    
    arma::uvec  tmp1    = unique(siretmp);
    arma::uvec  tmp2    = unique(damtmp);
    m                   = tmp1.n_rows + tmp2.n_rows - 1;
    arma::uvec  Sire    = zeros<uvec>(n + 1, 1);
    arma::uvec  Dam     = zeros<uvec>(n + 1, 1);
    arma::uvec  Rsire   = zeros<uvec>(m + 1, 1);
    arma::uvec  Rdam    = zeros<uvec>(m + 1, 1);
    arma::uvec  SID     = zeros<uvec>(n + 1, 1);
    arma::uvec  Link    = zeros<uvec>(n + 1, 1);
    arma::uvec  MaxIDP  = zeros<uvec>(m + 1, 1);
    arma::vec   F       = zeros(n + 1);
    arma::vec   B       = zeros(m + 1);
    arma::vec   x       = zeros(m + 1);
    /*
     * Copying tmp vectors of sires and dams into Sire and Dam
     */
    Sire.tail(n) = siretmp;
    Dam.tail(n)  = damtmp;
    
    F[0] = -1;
    for(rn = i = 1; i <= n; i++){
        SID[i]      = i;
        Link[i]     = 0;
        if(i <= m){
            x[i]    = 0.0;
        }
        S = Sire[i];
        D = Dam[i];
        if(S && !Link[S]){
            MaxIDP[rn]      = Link[S] = rn;
            Rsire[rn]       = Link[Sire[S]];
            Rdam[rn]        = Link[Dam[S]];
            rn++;
        }
        if(D && !Link[D]){
            Link[D]         = rn;
            Rsire[rn]       = Link[Sire[D]];
            Rdam[rn]        = Link[Dam[D]];
            rn++;
        }
        if(MaxIDP[Link[S]] < Link[D]){
            MaxIDP[Link[S]] = Link[D];
        }
    }
    // Sort animals according to ID of their sires into SId
    pairsort2(Sire, SID);
    for(k = i = 1; i <= n;){
        if(!Sire[SID[i]]){
            F[SID[i]] = 0.0;
            i++;
        } else {
            S       = Sire[SID[i]];
            rs      = Link[S];
            MIP     = MaxIDP[rs];
            x[rs]   = 1.0;
            for(; k <= S; k++){
                if(Link[k]){
                    B[Link[k]] = 0.5 - 0.25*(F[Sire[k]] + F[Dam[k]]);
                }
            }
            for(j = rs; j; --j){
                if(x[j]){
                    if(Sire[j]){
                        x[Rsire[j]] += x[j]*0.5;
                    }
                    if(Dam[j]){
                        x[Rdam[j]] += x[j]*0.5;
                    }
                    x[j] *= B[j];
                }
            }
            for(j = 1; j <= MIP; j++){
                x[j] += 0.5*(x[Rsire[j]] + x[Rdam[j]]);
            }
            for(; i <= n; i++){
                if(S != Sire[SID[i]]){
                    break;
                } else {
                    F[SID[i]] = x[Link[Dam[SID[i]]]]*0.5;
                }
            }
            for(j = 1; j <= MIP; j++){
                x[j] = 0.0;
            }
        }
    }
    return(F.tail(n));
}