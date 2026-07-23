#include "cbs.hpp"
#include "stats.hpp"
#include "correction.hpp"
#include <algorithm>
#include <random>

static void rotate_seg(const std::vector<double>& x, int s, int e, int r, std::vector<double>& out){ int n=e-s+1; out.resize(n); for(int i=0;i<n;++i){ int src=s+(i+r)%n; out[i]=x[src]; } }
static double max1_temp(const std::vector<double>& tmp){ std::vector<double> c; build_csum(tmp, c); auto ms=max_cusum_1cp(tmp,c,0,(int)tmp.size()-1); return ms.stat; }
static double max2_temp(const std::vector<double>& tmp, int min_seg){ std::vector<double> c; build_csum(tmp, c); auto ms=max_cusum_2cp(tmp,c,0,(int)tmp.size()-1,min_seg); return ms.stat; }

void segment_chromosome(const Series& S, const Args& a, std::vector<Segment>& out, std::mt19937_64& rng){
    const auto &x0 = S.x; std::vector<double> x = x0; // local copy (for correction)

    // Optional corrections
    if(a.center){ // per-chrom median center
        std::vector<double> tmp=x; std::nth_element(tmp.begin(), tmp.begin()+tmp.size()/2, tmp.end()); double med=tmp[tmp.size()/2]; if(tmp.size()%2==0){ std::nth_element(tmp.begin(), tmp.begin()+tmp.size()/2-1, tmp.end()); med=0.5*(tmp[tmp.size()/2]+tmp[tmp.size()/2-1]); }
        for(double &v:x) v-=med;
    }
    if(a.correction=="rolling") rolling_median_detrend(x, a.roll_window);
    else if(a.correction=="linear"){
        std::unordered_map<std::string,std::vector<double>> covars;
        if(!a.gc_col.empty()){
            auto it=S.covars.find(a.gc_col); if(it!=S.covars.end()) covars["GC"]=it->second; }
        if(!a.wave_col.empty()){
            auto it=S.covars.find(a.wave_col); if(it!=S.covars.end()) covars["WAVE"]=it->second; }
        std::vector<double> beta; if(!covars.empty()) linear_covariate_regress(x, covars, beta);
    }

    int N=(int)x.size(); if(N==0){ return; }
    std::vector<double> csum; build_csum(x, csum);

    struct Node{ int s,e,level; }; std::vector<Node> stack; stack.push_back({0,N-1,0});
    std::uniform_int_distribution<int> dist;

    while(!stack.empty()){
        Node nd=stack.back(); stack.pop_back(); int s=nd.s,e=nd.e,level=nd.level; int n=e-s+1;
        if(n<=2*a.min_seg_len || level>=a.max_depth){ double mean=seg_mean(csum,s,e); out.push_back({s,e,mean,level}); continue; }

        if(a.method=="2cp"){
            auto ms = max_cusum_2cp(x, csum, s, e, a.min_seg_len);
            if(ms.k1<0 || ms.k2<0){ double mean=seg_mean(csum,s,e); out.push_back({s,e,mean,level}); continue; }
            // Permutation test using circular rotations, test max 2-CP statistic
            double obs = ms.stat; int extreme=1; int total=a.perms+1; std::vector<double> tmp; tmp.reserve(n);
            dist = std::uniform_int_distribution<int>(0, n-1);
            for(int p=0;p<a.perms;++p){ int r = dist(rng); rotate_seg(x,s,e,r,tmp); double t=max2_temp(tmp, a.min_seg_len); if(t>=obs) ++extreme; }
            double pval=(double)extreme/total;
            if(pval<a.alpha){ stack.push_back({ms.k2+1,e,level+1}); stack.push_back({ms.k1+1,ms.k2,level+1}); stack.push_back({s,ms.k1,level+1}); }
            else { double mean=seg_mean(csum,s,e); out.push_back({s,e,mean,level}); }
        } else {
            auto ms = max_cusum_1cp(x, csum, s, e);
            if(ms.k<0){ double mean=seg_mean(csum,s,e); out.push_back({s,e,mean,level}); continue; }
            int nL=ms.k-s+1; int nR=e-(ms.k+1)+1; if(nL<a.min_seg_len || nR<a.min_seg_len){ double mean=seg_mean(csum,s,e); out.push_back({s,e,mean,level}); continue; }
            double obs=ms.stat; int extreme=1; int total=a.perms+1; std::vector<double> tmp; tmp.reserve(n);
            dist = std::uniform_int_distribution<int>(0, n-1);
            for(int p=0;p<a.perms;++p){ int r=dist(rng); rotate_seg(x,s,e,r,tmp); double t=max1_temp(tmp); if(t>=obs) ++extreme; }
            double pval=(double)extreme/total;
            if(pval<a.alpha){ stack.push_back({ms.k+1,e,level+1}); stack.push_back({s,ms.k,level+1}); }
            else { double mean=seg_mean(csum,s,e); out.push_back({s,e,mean,level}); }
        }
    }
}
