#include "correction.hpp"
#include <algorithm>
#include <cmath>

static double median_vec(std::vector<double> v){ if(v.empty()) return 0.0; size_t n=v.size(); std::nth_element(v.begin(), v.begin()+n/2, v.end()); double m=v[n/2]; if(n%2==0){ std::nth_element(v.begin(), v.begin()+n/2-1, v.end()); m=0.5*(v[n/2]+v[n/2-1]); } return m; }

void rolling_median_detrend(std::vector<double>& x, int window){
    if(window<3) return; if(window%2==0) window+=1; int n=(int)x.size(); std::vector<double> y(n);
    int h=window/2; std::vector<double> buf; buf.reserve(window);
    for(int i=0;i<n;++i){ int s=std::max(0,i-h); int e=std::min(n-1,i+h); buf.clear(); for(int j=s;j<=e;++j) buf.push_back(x[j]); double med=median_vec(buf); y[i]=x[i]-med; }
    x.swap(y);
}

// Ordinary Least Squares for small #covariates (<=4 typical): X=[1, c1, c2, ...]
static void ols_solve(const std::vector<std::vector<double>>& X, const std::vector<double>& y, std::vector<double>& beta){
    // Solve normal equations X^T X beta = X^T y with tiny dimension — use naive Gauss-Jordan
    int n = (int)y.size(); int p = (int)X.size(); // X stored row-major by column vectors here: X[p][n]
    // Build XtX and Xty
    std::vector<std::vector<double>> A(p, std::vector<double>(p,0.0));
    std::vector<double> b(p,0.0);
    for(int i=0;i<p;++i){ for(int j=0;j<p;++j){ double s=0; for(int k=0;k<n;++k) s += X[i][k]*X[j][k]; A[i][j]=s; } }
    for(int i=0;i<p;++i){ double s=0; for(int k=0;k<n;++k) s += X[i][k]*y[k]; b[i]=s; }
    // Augment
    for(int i=0;i<p;++i){ A[i].push_back(b[i]); }
    // Gauss-Jordan
    for(int i=0;i<p;++i){
        // pivot
        int piv=i; for(int r=i+1;r<p;++r) if (std::abs(A[r][i])>std::abs(A[piv][i])) piv=r; std::swap(A[i],A[piv]);
        double diag=A[i][i]; if(std::abs(diag)<1e-12) continue; for(size_t c=0;c<A[i].size();++c) A[i][c]/=diag;
        // eliminate
        for(int r=0;r<p;++r){ if(r==i) continue; double f=A[r][i]; if(std::abs(f)<1e-18) continue; for(size_t c=0;c<A[r].size();++c) A[r][c]-=f*A[i][c]; }
    }
    beta.assign(p,0.0); for(int i=0;i<p;++i) beta[i]=A[i][p];
}

void linear_covariate_regress(std::vector<double>& x, const std::unordered_map<std::string,std::vector<double>>& covars,
                              std::vector<double>& out_beta){
    if(covars.empty()) return; int n=(int)x.size();
    // Construct design with intercept + each covariate normalized to mean 0, var 1
    std::vector<std::vector<double>> X; X.reserve(covars.size()+1);
    std::vector<double> ones(n,1.0); X.push_back(std::move(ones));
    for(const auto& kv: covars){ const auto& v = kv.second; if((int)v.size()!=n) continue; double mu=0, s2=0; for(double z:v) mu+=z; mu/=n; for(double z:v) s2+=(z-mu)*(z-mu); s2 = (s2<=0?1.0:std::sqrt(s2/n)); std::vector<double> z(n); for(int i=0;i<n;++i) z[i]=(v[i]-mu)/s2; X.push_back(std::move(z)); }
    std::vector<double> beta; ols_solve(X,x,beta);
    // subtract fitted
    std::vector<double> fit(n,0.0);
    for(size_t j=0;j<X.size();++j){ for(int i=0;i<n;++i) fit[i]+=beta[j]*X[j][i]; }
    for(int i=0;i<n;++i) x[i]-=fit[i];
    out_beta = std::move(beta);
}
