class pileup_corrector {
public:
    double get_pileup_weight(int num_pu_vertices);
    pileup_corrector(bool isMC,std::string era);
    ~pileup_corrector();
private:
    void fail(int line);
    TH1D* pileup_ratio_;
    bool isMC_;
    int year_;
};
