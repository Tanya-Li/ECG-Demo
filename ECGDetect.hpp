//Tianyang Li, V00814119

#ifndef ECGDetect__hpp
#define ECGDetect__hpp

#include <string>
#include <vector>
//ECGDetect class design -- for parameters extraction
class ECGDetect {
public:
    ECGDetect(std::string, int);
    ~ECGDetect();
    
    void load_ECG_Data();
    size_t get_ecg_data_length() const;
    void setFS(int);
    int getFS() const;
    
    std::vector<int> get_R_i() const;
    std::vector<int> get_S_i() const;
    std::vector<int> get_T_i() const;
    std::vector<int> get_Q_i() const;
    std::vector<int> get_S_end() const;
    std::vector<int> get_T_end() const;
    
    std::vector<double> get_buffer_plot() const;
    
    double get_HR() const;
    double get_QRS_duration() const;
    double get_QT_interval() const;
    double get_PR_interval() const;
    double get_QTc() const;
    
    void run();
    
private:
    void _detect();
    
private:
    //difine parameters
    std::string _ecg_path;
    int _fs;
    
    double* _ecg_data;
    size_t _ecg_data_length;
    
    double _HR;
    double _QRS_duration;
    double _QT_interval;
    double _PR_interval;
    double _QTc;
    
    std::vector<int> _R_i;
    std::vector<int> _S_i;
    std::vector<int> _T_i;
    std::vector<int> _Q_i;
    std::vector<int> _S_end;
    std::vector<int> _T_end;
    //save plot information in a vector
    std::vector<double> _buffer_plot;
};

#endif 
