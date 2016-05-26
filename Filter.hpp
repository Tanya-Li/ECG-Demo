//Tianyang Li, V00814119

#ifndef Filter__hpp
#define Filter__hpp

#include <vector>
#include <cstdlib>

const double PI = 3.14159265358979323846;

//Filter class design
class Filter {
public:
    Filter(int, double, double);
    ~Filter();
    
    void ecg_filter(double*, size_t);
private:
    int _filter_length; //filter coef length
    double _filter_para_low; //bandpass low para
    double _filter_para_high; //bandpass high para
    std::vector<double> _filter_coef; //bandpass filter coefficiences
};


#endif
