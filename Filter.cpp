//Tianyang Li, V00814119

#include <cmath>
#include "Filter.hpp"
#include <cstdlib>

//Digital filter design part- bandpass fileter
Filter::Filter(int length, double para_low, double para_high) {
    _filter_length = length;
    //parameters for bandpass filter design
    _filter_para_low = para_low * PI;
    _filter_para_high = para_high * PI;
    //filter coefficients
    for (int i = 0; i < _filter_length / 2; ++i) {
        _filter_coef.push_back(((sin((i - _filter_length / 2) * _filter_para_high) -
                                 sin((i - _filter_length / 2) * _filter_para_low)) /
                                ((i - _filter_length / 2) * PI) ) *
                               (0.54 - 0.46 * cos(PI * i / (_filter_length / 2))));
    }
    
    _filter_coef.push_back((_filter_para_high - _filter_para_low) / PI);
    
    _filter_coef.insert(_filter_coef.end(), _filter_coef.rbegin() + 1, _filter_coef.rend());
}

Filter::~Filter() {
    
}
//ECG data goes through bandpass filter -- coefficients convolution
void Filter::ecg_filter(double* ecg_data, size_t length) {
    double tmp;
    for (size_t i = length - 1; i >= _filter_length - 1; --i) {
        tmp = 0;
        for (int j = 0; j < _filter_length; ++j) {
            tmp += ecg_data[i - j] * _filter_coef[j];
        }
        ecg_data[i] = tmp;
    }
    double sum = 0;
    for (size_t i = 0; i < length; ++i) {
        sum += (ecg_data[i]/(double)length);
    }
    for (size_t i = 0; i < length; ++i) {
        ecg_data[i] -= sum;
    }
}
