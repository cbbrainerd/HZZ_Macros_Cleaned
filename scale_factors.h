#include <TH2.h>
#include <string>

class scale_factors {
public:
    float get_scale_factor(double eta,double pt,bool isGap=false);
    scale_factors(bool isMC,std::string era,std::string particle_type);
    ~scale_factors();
private:
    void select_scale_factor_files(const char* gap,const char* no_gap=nullptr);
    float overflow_pt;
    TH2F* scale_factor_no_gap;
    TH2F* scale_factor_is_gap;
    const char* SFHistName;
protected:
    bool isMC_;
};

class scale_factors_and_efficiencies : public scale_factors {
public:
    scale_factors_and_efficiencies(bool isMC,std::string era,std::string particle_type);
    float get_efficiency(double eta,double pt);
    ~scale_factors_and_efficiencies();
private:
    void select_efficiency_files(const char* low_pt,const char* high_pt);
    TH2F* efficiency_low_pt;
    TH2F* efficiency_high_pt;
    const char* EffHistName;
};
