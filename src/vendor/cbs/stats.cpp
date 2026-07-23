#include "stats.hpp"
#include <cmath>
#include <algorithm>

void build_csum(const std::vector<double>& x, std::vector<double>& csum){
    csum.assign(x.size()+1, 0.0);
    for(size_t i=0;i<x.size();++i) csum[i+1]=csum[i]+x[i];
}

// Standardized 1-CP max CUSUM
MaxStat max_cusum_1cp(const std::vector<double>& x, const std::vector<double>& csum, int s, int e){
    int n = e-s+1; if (n<2) return {0.0,-1};
    double mean = seg_mean(csum,s,e);
    double best=-1.0; int bestk=-1; double pref=0.0;
    for(int i=s;i<=e-1;++i){ pref += x[i]-mean; int L=i-s+1; int R=n-L; double denom = std::sqrt((double)L*R/n); double t=std::abs(pref)/denom; if(t>best){best=t; bestk=i;} }
    return {best,bestk};
}

// Approximate 2-CP scan: search k1<k2 maximizing between-segment variability vs pooled
MaxStat2 max_cusum_2cp(const std::vector<double>& x, const std::vector<double>& csum, int s, int e, int min_seg_len){
    int n=e-s+1; if(n<3*min_seg_len) return {0.0,-1,-1};
    double best=-1.0; int bk1=-1,bk2=-1; double totMean=seg_mean(csum,s,e);
    // Coarse grid to keep O(n^2) manageable; scan every step=1 for n<=5000 else step>1
    int step = (n<=5000)?1: (n<=20000?2:5);
    for(int k1=s+min_seg_len-1; k1<=e-2*min_seg_len+1; k1+=step){
        for(int k2=k1+min_seg_len; k2<=e-min_seg_len+1; k2+=step){
            // means
            double m1 = seg_mean(csum,s,k1);
            double m2 = seg_mean(csum,k1+1,k2);
            double m3 = seg_mean(csum,k2+1,e);
            int n1 = k1-s+1; int n2 = k2-(k1+1)+1; int n3 = e-(k2+1)+1;
            double between = n1*(m1-totMean)*(m1-totMean) + n2*(m2-totMean)*(m2-totMean) + n3*(m3-totMean)*(m3-totMean);
            // Use sqrt(between) as monotone proxy; more stable under permutations
            double t = std::sqrt(between);
            if(t>best){ best=t; bk1=k1; bk2=k2; }
        }
    }
    return {best,bk1,bk2};
}
