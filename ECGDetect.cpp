//Tianyang Li, V00814119

#include <iostream>
#include <fstream>
#include <cmath>
#include "ECGDetect.hpp"
#include "Filter.hpp"

//constructor
ECGDetect::ECGDetect(std::string path, int fs):
_ecg_path(path), _fs(fs), _ecg_data(0) /*nullptr*/ , _ecg_data_length(0), _HR(0),
_QRS_duration(0), _QT_interval(0), _PR_interval(0), _QTc(0){
    
}
//destructor
ECGDetect::~ECGDetect() {
    if (_ecg_data)
        delete [] _ecg_data;
}
//read ECG data in, bin file
void ECGDetect::load_ECG_Data() {
    size_t file_len = -1;
    std::ifstream infile;
    char* buffer = 0; // buffer
    // ecg_path
    infile.open(_ecg_path.c_str(), std::ios::in | std::ios::binary | std::ios::ate);
    
    if(infile.is_open()) {
        
        infile.seekg(0, std::ios::end);
        file_len = infile.tellg();
        infile.seekg(0, std::ios::beg);
        
        buffer = new char[file_len];
        
        if (!buffer) {
            std::cerr<<"Memory allocation error!\n";
            infile.close();
            return;
        }
        
        infile.read(buffer, file_len);
        infile.close();
        
        _ecg_data = new double[file_len];
        for (size_t i = 0; i < file_len; ++i) {
            _ecg_data[i] = buffer[i] & 0x7F + (buffer[i] & 0x80 ? 128 : 0);
        }
        
        _ecg_data_length = file_len;
        
        delete [] buffer;
    }
    else {
        std::cerr<<"Unable to open ecg file\n";
    }
}
//obtain data length
size_t ECGDetect::get_ecg_data_length() const {
    return _ecg_data_length;
}
//set sampling frequency
void ECGDetect::setFS(int fs) {
    _fs = fs;
}
//get sampling frequency
int ECGDetect::getFS() const {
    return _fs;
}
//obtain R point values
std::vector<int> ECGDetect::get_R_i() const {
    return _R_i;
}
//obtain S point values
std::vector<int> ECGDetect::get_S_i() const {
    return _S_i;
}
//obtain Q point values
std::vector<int> ECGDetect::get_Q_i() const {
    return _Q_i;
}
//obtain T point values
std::vector<int> ECGDetect::get_T_i() const {
    return _T_i;
}
//obtain S end point values
std::vector<int> ECGDetect::get_S_end() const {
    return _S_end;
}
//obtain T end point values
std::vector<int> ECGDetect::get_T_end() const {
    return _T_end;
}
//obtain value for plot
std::vector<double> ECGDetect::get_buffer_plot() const {
    return _buffer_plot;
}
//obtain heart rate value
double ECGDetect::get_HR() const {
    return _HR;
}
//obtain QRS values
double ECGDetect::get_QRS_duration() const {
    return _QRS_duration;
}
//obtain QT interval value
double ECGDetect::get_QT_interval() const {
    return _QT_interval;
}
//obtain PR interval value
double ECGDetect::get_PR_interval() const {
    return _PR_interval;
}
//obtain QTc value
double ECGDetect::get_QTc() const {
    return _QTc;
}

void ECGDetect::run() {
    load_ECG_Data();
    _detect();
}
//important data processing part 
void ECGDetect::_detect() {
    //set parameters initial value
    double time_scale = _ecg_data_length / _fs;
    int window = _fs / 25;
    int state = 0;
    double weight = 1.8;
    int dum = 0;
    int co = 0;
    int S_on = 0;
    int c = 0;
    int T_on = 0;
    int T_on1 = 0;
    int sleep = 0;
    double slope_q = 0;
    
    double buffer_T;
    double buffer_mean;
    std::vector<double> buffer_base;
    double buffer_long;
    double mean_online;
    double current_max = 0;
    int ind = 0;
    double thres2;
    int min_QRST_length;
    int ecg_length;
    int dist_ST;
    int I = 0;
    double min_slope;
    int qi_index = 0;
    double slope;
    int I0;
    double slope_p;
    int po_index = 0;
    int uni_dim;
    
    std::vector<double> R_amp;
    std::vector<double> S_amp;
    std::vector<double> T_amp;
    std::vector<double> Q_amp;
    std::vector<double> S_end_amp;
    std::vector<double> T_end_amp;
    std::vector<int> P_peak;
    std::vector<double> P_peak_amp;
    std::vector<double> thres2_p;
    std::vector<int> thres2_p_i;
    std::vector<int> ST_test;
    std::vector<double> grad_one;
    std::vector<int> Q_onset;
    std::vector<double> Q_onset_amp;
    std::vector<int> indexh;
    std::vector<int> index_order;
    std::vector<double> height;
    std::vector<int> dist;
    std::vector<double> max_slope;
    std::vector<int> index_max_slope;
    std::vector<double> height_diff;
    std::vector<int> dist_x;
    std::vector<int> P_onset;
    std::vector<double> P_onset_amp;
    std::vector<int> PR_mat;
    std::vector<int> PR_inter;
    std::vector<int> QT_mat;
    std::vector<int> QRS_mat;
    std::vector<int> QT_inter;
    std::vector<int> QRS_inter;
    
    int idx;
    int jdx;
    int count;
    
    // ----------------------------------------------------------------------------//
    //date goes through designed filters
    Filter filter(21, 0.0006, 0.09);
    filter.ecg_filter(_ecg_data, _ecg_data_length);
    buffer_T = 0;
    for (idx = 0, count = 0; idx < 2 * _fs; ++idx, ++count) {
        buffer_T = buffer_T * count / (count + 1) + _ecg_data[idx] / (count + 1);
    }
    buffer_mean = 0;
    for (idx = 0, count = 0; idx < 2 * _fs; ++idx, ++count) {
        buffer_mean = buffer_mean * count / (count + 1) + fabs(_ecg_data[idx] - buffer_T) / (count + 1);
    }
    
    buffer_long = 0;
    
    // start online inference (Assuming the signal is being acquired online)
    for (idx = 0, count = 1; idx < _ecg_data_length; ++idx, ++count) {
        
        buffer_long += _ecg_data[idx];
        
        // Renew the mean and adapt it to the signal after 1 second of processing
        if (count == 2 * _fs) {
            count = 0;
            buffer_T = 0;
            for (jdx = 0; jdx < 2 * _fs; ++jdx) {
                buffer_T = buffer_T * jdx / (jdx + 1) + _ecg_data[idx - jdx] / (jdx + 1);
            }
            buffer_mean = 0;
            for (jdx = 0; jdx < 2 * _fs; ++jdx) {
                buffer_mean = buffer_mean * jdx / (jdx + 1) + fabs(_ecg_data[idx - jdx] - buffer_T) / (jdx + 1);
            }
        }
        
        // smooth the signal by taking the average of 15 samples and add the new upcoming samples
        if (idx >= window - 1) {
            mean_online = buffer_long / window;
            _buffer_plot.push_back(mean_online);
            
            // Enter state 1 (putative R wave) as soon as that the mean exceeds the double time of threshold
            if (state == 0) {
                if (_buffer_plot.size() >= 3) {
                    if ((mean_online > buffer_mean * weight) && (_buffer_plot[idx - 1 - window] > _buffer_plot[idx - window])) {
                        state = 1;
                        current_max = _buffer_plot[idx - 1 - window];
                        ind = idx - 1 - window;
                        
                    }
                    else {
                        state = 0;
                    }
                }
            }
            
            // Locate the R wave location by finding the highest local maximum
            if (state == 1) {
                if (current_max > _buffer_plot[idx - window]) {
                    dum++;
                    if (dum > 4) {
                        _R_i.push_back(ind);
                        R_amp.push_back(_buffer_plot[ind]);
                        
                        // Locate Q wave
                        if (ind > (int)(0.04 * _fs + 0.5)) {
                            int Q_ti = ind;
                            double Q_tamp = _buffer_plot[ind];
                            for (jdx = 1; jdx <= (int)(0.04 * _fs + 0.5); ++jdx) {
                                if (Q_tamp > _buffer_plot[ind - jdx]) {
                                    Q_tamp = _buffer_plot[ind - jdx];
                                    Q_ti = ind - jdx;
                                }
                            }
                            _Q_i.push_back(Q_ti);
                            Q_amp.push_back(Q_tamp);
                        }
                        
                        if (_R_i.size() > 8) {
                            weight = 0;
                            for (jdx = 0; jdx < 8; ++jdx) {
                                weight += 0.3 * R_amp[_R_i.size() - jdx - 1];
                            }
                            weight = weight / buffer_mean / 8;
                        }
                        
                        state = 2;
                        dum = 0;
                    }
                }
                else {
                    dum = 0;
                    state = 0;
                }
            }
            
            // Check whether the signal drops below the threshold to look for S wave
            if (state == 2) {
                if (mean_online <= buffer_mean) {
                    state = 3;
                }
            }
            
            // Enter S wave detection state 3 (S detection)
            if (state == 3) {
                co++;
                if (co < (int)(0.2 * _fs + 0.5)) {
                    if (_buffer_plot[idx - window - 1] <= _buffer_plot[idx - window]) {
                        S_on++;
                        if (S_on >= (int)(0.012 * _fs + 0.5)) {
                            _S_i.push_back(idx - window - 4);
                            S_amp.push_back(_buffer_plot[idx - window - 4]);
                            state = 4;
                            S_on = 0;
                            co = 0;
                        }
                    }
                }
                else {
                    state = 4;
                    co = 0;
                }
            }
            
            // enter state 4 possible T wave detection
            if (state == 4) {
                if (mean_online < buffer_mean) {
                    state = 6;
                }
            }
            
            // Enter State 6 which is T wave possible detection
            if (state == 6) {
                c++;
                if (c <= (int)(0.7 * _fs + 0.5)) {
                    thres2 = fabs(fabs(buffer_T) - fabs(S_amp[_S_i.size() - 1])) * 0.75 + S_amp[_S_i.size() - 1];
                    thres2_p.push_back(thres2);
                    thres2_p_i.push_back(ind);
                    if (mean_online > thres2) {
                        T_on++;
                        if (T_on >= (int)(0.012 * _fs + 0.5)) {
                            if (_buffer_plot[idx - window - 1] >= _buffer_plot[idx - window]) {
                                T_on1++;
                                if (T_on1 > (int)(0.032 * _fs + 0.5)) {
                                    _T_i.push_back(idx - window - 11);
                                    T_amp.push_back(_buffer_plot[idx - window - 11]);
                                    state = 5;
                                    T_on = 0;
                                    T_on1 = 0;
                                }
                            }
                        }
                    }
                }
                else {
                    state = 5;
                }
            }
            
            // This state is for avoiding the detection of a highly variate noise or another peak
            // this avoids detection of two peaks R waves less than half a second
            if (state == 5) {
                sleep = sleep + c + 1;
                c = 0;
                if (sleep * 1.0 / _fs >= 0.4) {
                    state = 0;
                    sleep = 0;
                }
            }
            
            buffer_long -= _ecg_data[idx - window + 1];
        }
    }
    
    // Self-design Error control part 1
    if (_R_i[0] < _Q_i[0]) {
        for (idx = 0; idx < _R_i.size() - 1; ++idx)
            _R_i[idx] = _R_i[idx + 1];
        for (idx = 0; idx < _S_i.size() - 1; ++idx)
            _S_i[idx] = _S_i[idx + 1];
        for (idx = 0; idx < _T_i.size() - 1; ++idx)
            _T_i[idx] = _T_i[idx + 1];
    }
    // Self-design Error control part 2
    min_QRST_length = _R_i.size();
    if (min_QRST_length < _S_i.size())
        min_QRST_length = _S_i.size();
    if (min_QRST_length < _T_i.size())
        min_QRST_length = _T_i.size();
    if (min_QRST_length < _Q_i.size())
        min_QRST_length = _Q_i.size();
    min_QRST_length--;
    ecg_length = _T_i[min_QRST_length - 1] + 120;
    
    // detect the end of s wave
    for (idx = 0; idx < min_QRST_length; ++idx) {
        dist_ST = (int)ceil((_T_i[idx] - _S_i[idx]) / 2.5);
        ST_test.push_back(dist_ST + _S_i[idx]);
        // compute each gradient during ST_test
        for (jdx = 3; jdx <= dist_ST - 2; ++jdx) {
            if (_S_i[idx] + jdx + 2 <= ecg_length) {
                if (jdx - 3 >= grad_one.size()) {
                    grad_one.push_back(fabs(-2 * _ecg_data[_S_i[idx] + jdx - 2] - _ecg_data[_S_i[idx] + jdx - 1]
                                            + _ecg_data[_S_i[idx] + jdx + 1] + 2 * _ecg_data[_S_i[idx] + jdx + 2]));
                }
            }
        }
        if (grad_one.size() > 0) {
            I = 0;
            for (jdx = 1; jdx < grad_one.size(); ++jdx) {
                if (grad_one[jdx] < grad_one[I]) {
                    I = jdx;
                }
            }
        }
        if (_S_i[idx] + I + 2 <= ecg_length) {
            _S_end.push_back(_S_i[idx]);
            S_end_amp.push_back(_S_end[_S_end.size() - 1]);
        }
    }
    
    // Q onset detection --window 40ms
    for (idx = 0; idx < min_QRST_length; ++idx) {
        min_slope = 10;
        slope_q = 999999.99;
        for (jdx = 3; jdx <= 30 - 2; ++jdx) {
            if (_Q_i[idx] - jdx - 2 >= 0) {
                slope_q = fabs(-2 * _ecg_data[_Q_i[idx] - jdx - 2] - _ecg_data[_Q_i[idx] - jdx - 1]
                               + _ecg_data[_Q_i[idx] - jdx + 1] + 2 * _ecg_data[_Q_i[idx] - jdx + 2]);
            }
            if (slope_q < min_slope) {
                min_slope = slope_q;
                qi_index = jdx;
            }
        }
        Q_onset.push_back(_Q_i[idx] - qi_index);
        Q_onset_amp.push_back(_ecg_data[Q_onset[idx]]);
    }
    
    // preparation work
    indexh.push_back(thres2_p_i[0]);
    index_order.push_back(1);
    for (idx = 0; idx < thres2_p.size() - 1; ++idx) {
        if (thres2_p_i[idx] != thres2_p_i[idx + 1]) {
            indexh.push_back(thres2_p_i[idx + 1]);
            index_order.push_back(idx + 1);
        }
    }
    
    for (idx = 0; idx < indexh.size(); ++idx) {
        height.push_back(thres2_p[index_order[idx]]);
    }
    
    // T wave end detection
    for (idx = 0; idx < min_QRST_length; ++idx) {
        dist.push_back((int)ceil(0.65 * (_T_i[idx] - _S_i[idx])));
        max_slope.clear();
        for (jdx = 0; jdx < _T_i.size(); ++jdx) {
            max_slope.push_back(0.0);
        }
        
        index_max_slope.push_back(0);
        for (jdx = 0; jdx < dist[idx]; ++jdx) {
            if (_T_i[idx] + jdx + 1 < _ecg_data_length) {
                slope = fabs(_ecg_data[_T_i[idx] + jdx] - _ecg_data[_T_i[idx] + jdx + 1]);
                if (slope > max_slope[idx]) {
                    max_slope[idx] = slope;
                    index_max_slope[idx] = jdx;
                }
            }
        }
        height_diff.push_back(_ecg_data[_T_i[idx] + index_max_slope[idx]] - height[idx]);
        dist_x.push_back((int)ceil(height_diff[idx] / max_slope[idx]));
        if (_T_i[idx] + index_max_slope[idx] + dist_x[idx] <= ecg_length) {
            _T_end.push_back(_T_i[idx] + index_max_slope[idx] + dist_x[idx]);
            T_end_amp.push_back(_ecg_data[_T_end[_T_end.size() - 1]]);
        }
    }
    
    // P peak detection
    for (idx = 0; idx < min_QRST_length; ++idx) {
        if (_R_i[idx] >= 250) {
            I0 = _R_i[idx];
            for (jdx = 0; jdx < 250; ++jdx) {
                if (_ecg_data[_R_i[idx] - jdx - 1] > _ecg_data[I0]) {
                    I0 = _R_i[idx] - jdx - 1;
                }
            }
        }
        else {
            continue;
        }
        P_peak.push_back(I0);
        P_peak_amp.push_back(_ecg_data[I0]);
    }
    
    // P onset detection based on time range threshold
    for (idx = 0; idx < min_QRST_length; ++idx) {
        min_slope = 10;
        for (jdx = 1; jdx <= 35; ++jdx) {
            slope_p = fabs(_ecg_data[P_peak[idx] - 80 + jdx] - _ecg_data[P_peak[idx] - 80 + jdx + 1]);
            if (slope_p < min_slope) {
                min_slope = slope_p;
                po_index = jdx;
            }
        }
        P_onset.push_back(P_peak[idx] - 80 + po_index);
        P_onset_amp.push_back(_ecg_data[P_onset[idx]]);
    }
    //-----------------------------------------------------------------------------------------//
    //Error control on parameters
    uni_dim = _S_end.size() < _T_end.size() ? _S_end.size() : _T_end.size();
    
    if (Q_onset[0] < P_onset[0]) {
        for (idx = 0; idx < uni_dim - 1; ++idx) {
            PR_mat.push_back(Q_onset[idx + 1] - P_onset[idx]);
        }
        for (idx = 0; idx < uni_dim - 1; ++idx) {
            if (60 < PR_mat[idx] && PR_mat[idx] < 250) {
                PR_inter.push_back(PR_mat[idx]);
            }
        }
    }
    else {
        for (idx = 0; idx < uni_dim; ++idx) {
            PR_mat.push_back(Q_onset[idx] - P_onset[idx]);
        }
        for (idx = 0; idx < uni_dim; ++idx) {
            if (60 < PR_mat[idx] && PR_mat[idx] < 250) {
                PR_inter.push_back(PR_mat[idx]);
            }
        }
    }
    
    for (idx = 0; idx < uni_dim; ++idx) {
        QT_mat.push_back(_T_end[idx] - Q_onset[idx]);
        QRS_mat.push_back(_S_end[idx] - Q_onset[idx]);
    }
    
    for (idx = 0; idx < uni_dim; ++idx) {
        if (QT_mat[idx] < 420) {
            QT_inter.push_back(QT_mat[idx]);
        }
        if (30 < QRS_mat[idx] && QRS_mat[idx] < 180) {
            QRS_inter.push_back(QRS_mat[idx]);
        }
    }
    //-----------------------------------------------------------------------------------------//
    // compute heart rate
    _HR = _R_i.size() / time_scale * 60;
    // compute PR interval
    _PR_interval = 0.0;
    if (PR_inter.size() > 0) {
        for (idx = 0; idx < PR_inter.size(); ++idx) {
            _PR_interval = _PR_interval * idx / (idx + 1) + PR_inter[idx] * 1.0 / (idx + 1);
        }
    }
    // compute QRS duration
    _QRS_duration = 0.0;
    if (QRS_inter.size() > 0) {
        for (idx = 0; idx < QRS_inter.size(); ++idx) {
            _QRS_duration = _QRS_duration * idx / (idx + 1) + QRS_inter[idx] * 1.0 / (idx + 1);
        }
    }
    // compute QT interval
    _QT_interval = 0.0;
    if (QT_inter.size() > 0) {
        for (idx = 0; idx < QT_inter.size(); ++idx) {
            _QT_interval = _QT_interval * idx / (idx + 1) + QT_inter[idx] * 1.0 / (idx + 1);
        }
    }
    // compute QTc based on Bazett
    _QTc = _QT_interval / (sqrt(60 / _HR));
}
