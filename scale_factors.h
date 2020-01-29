#include <TH1.h>
#include <string>

class scale_factors {
public:
    float get_scale_factor(double eta,double pt,bool isGap=false);
    scale_factors(bool isMC,std::string era,std::string particle_type);
    ~electron_scale_factors();
private:
    void select_scale_factor_files(const char*,const char*);
    TH2F* scale_factor_not_gap;
    TH2F* scale_factor_is_gap;
protected:
    bool isMC_;
};

class scale_factors_and_efficiencies : public scale_factors {
public:
    scale_factors_and_efficiencies(bool isMC,std::string era,std::string particle_type);
    float get_efficiency(double eta,double pt);
private:
    void select_efficiency_files(const char* low_et,const char* high_et);
    TH2F* efficiency_low_et;
    TH2F* efficiency_high_et;
};
