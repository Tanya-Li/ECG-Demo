//Tianyang Li, V00814119

#include <iostream>
#include <fstream>
#include <cmath>
#include "FatalDetect.hpp"
#include "Filter.hpp"

//FatalDetect class constuctor
FatalDetect::FatalDetect(std::string path, int fs):
_ecg_path(path), _fs(fs), _ecg_data(0),_ecg_data_length(0), _fatalResult(false) {
    
}
//Destructor
FatalDetect::~FatalDetect() {
    delete [] _ecg_data;
}

//function for read data in
void FatalDetect::load_ECG_Data() {
    std::string line;
    std::ifstream myfile (_ecg_path);
    
    if(myfile.is_open()) {
        while (getline(myfile, line)){
            double temp = atof(line.c_str());
            fataldata.push_back(temp);
        }
        myfile.close();
    }else {
        std::cerr << "File can't be open successfully";
    }
    _ecg_data_length = fataldata.size();
    _ecg_data = new double [_ecg_data_length];
    for (size_t i = 0; i < _ecg_data_length; ++i) {
        _ecg_data[i] = fataldata[i];
    }
}
//obtian data in a vector
std::vector<double> FatalDetect::get_buffer_plot() const {
    return fataldata;
}
//get data length
size_t FatalDetect::get_ecg_data_length() const {
    return _ecg_data_length;
}
//set data sampling rate
void FatalDetect::setFS(int fs) {
    _fs = fs;
}
//get data sampling rate
int FatalDetect::getFS() const {
    return _fs;
}
//get the result of detection process
bool FatalDetect::fatal() {
    return _fatalResult;
}

void FatalDetect::run() {
    load_ECG_Data();
    fatal_detect();
}
//key detection code -- idea from relevant papers
void FatalDetect::fatal_detect() {
    
    //preprocessing - moving average and so on
    double ecg_mean = 0;
    if (_ecg_data_length > 0) {
        for (int idx = 0; idx < _ecg_data_length; ++idx) {
            ecg_mean = ecg_mean * idx / (idx + 1) + _ecg_data[idx] * 1.0 / (idx + 1);
        }
    }
    
    //base_line wander elimination
    for (int i = 0; i < _ecg_data_length; ++i) {
        _ecg_data[i] -= ecg_mean;
    }
    
    int Nw = 3;
    //moving window averaging filter
    std::vector<double> ecg_s;
    
    ecg_s.push_back(_ecg_data[0]);
    
    for (int idx = 1; idx < _ecg_data_length-1; ++idx) {
        double ecg_tmp = (1.0/Nw)*(_ecg_data[idx-1]+_ecg_data[idx]+_ecg_data[idx+1]);
        ecg_s.push_back(ecg_tmp);
    }
    ecg_s.push_back(_ecg_data[_ecg_data_length-1]);
    
    for (int i = 1; i < _ecg_data_length-1; ++i) {
        _ecg_data[i] = ecg_s[i-1];
    }
    
    //process with designed filters
    Filter filter(21, 0.008, 0.24);
    filter.ecg_filter(_ecg_data, _ecg_data_length);
    
    std::vector<double> ecgn;
    std::vector<double> mav_a;
    std::vector<double> ecgs;
    std::vector<double> mav;
    std::vector<double> ecg_section;
    
    //normalized ECG signal
    int fs = getFS();
    int l = (int)(_ecg_data_length/fs);
    int num_section = (int)(l/8);
    double sum;
    double mav_sum;
    
    for (int idx =1; idx <= num_section; ++idx) {
        ecgs.clear();
        
        for (int i = 0; i < 8*fs; ++i) {
            ecgs.push_back(_ecg_data[(idx-1)*8*fs+i]);
        }
        mav.clear();
        
        for (int k = 1; k < 8; ++k) {
            ecg_section.clear();
            double max_2s = 0;
            for (int kdx = 0; kdx < 2*fs; kdx++) {
                ecg_section.push_back(std::abs(ecgs[(k-1)*fs+kdx]));
                max_2s = std::max(max_2s, ecg_section[kdx]);
            }
            for (int g = 0; g < ecg_section.size(); ++g) {
                ecg_section[g] = ecg_section[g]/max_2s;
            }
            sum = 0;
            for (int h=0; h < ecg_section.size(); h++) {
                sum += ecg_section[h];
            }
            sum = sum/(2*fs);
            mav.push_back(sum);
        }
        mav_sum = 0;
        for (int jdx = 0; jdx < 7; ++jdx) {
            mav_sum += mav[jdx];
        }
        mav_a.push_back((mav_sum/7.0));
    }
    
    //obtain the time period of VF signal
    
    int VFcount = 0;
    for (int tdx = 0; tdx < mav_a.size(); ++tdx) {
        if (mav_a[tdx] > 0.27) {
            ++VFcount;
        }
    }
    //if VF signal found, this is a fatal signal
    if (VFcount){
        _fatalResult = true;
    }
    
}

