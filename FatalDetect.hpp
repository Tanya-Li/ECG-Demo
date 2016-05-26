
//Tianyang Li, V00814119

#ifndef FatalDetect__hpp
#define FatalDetect__hpp

#include <string>
#include <vector>

class FatalDetect {
public:
    FatalDetect(std::string, int);
    ~FatalDetect();
    
    void load_ECG_Data();
    size_t get_ecg_data_length() const;
    void setFS(int);
    int getFS() const;
    bool fatal();
    std::vector<double> get_buffer_plot() const;
    void run();
    
private:
    void fatal_detect();
    
private:
    std::string _ecg_path;
    int _fs;
    size_t _ecg_data_length;
    
    double* _ecg_data;
    std::vector<double> fataldata;
    
    bool _fatalResult; //if false - not fatal signal
    
};


#endif


