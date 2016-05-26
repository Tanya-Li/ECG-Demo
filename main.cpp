
//Tianyang Li, V00814119

#include <iostream>
#include <GL/glut.h>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>
#include <algorithm>
#include "ECGDetect.hpp"
#include "FatalDetect.hpp"

//define global variables for output
double HR;
double QRS;
double PR;
double QTc;
double dataLen;
int samplingRate;
std::vector<double> ecgDisplay;
std::string dataName("ecg1.bin");
int keyHit = 0;
bool fatalSignal = false;

//For text output-style 1
void drawText(const char* text, int length, GLfloat x, GLfloat y){
    glMatrixMode(GL_PROJECTION);
    double *matrix = new double[16];
    glGetDoublev(GL_PROJECTION_MATRIX, matrix);
    glLoadIdentity();
    glOrtho(0, 800, 0, 600, -5, -5);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glPushMatrix();
    glLoadIdentity();
    glRasterPos2f(x, y);
    
    for(int i=0; i<length; i++){
        glutBitmapCharacter(GLUT_BITMAP_8_BY_13, text[i]);
    }
    glPopMatrix();
    glMatrixMode(GL_PROJECTION);
    glLoadMatrixd(matrix);
    glMatrixMode(GL_MODELVIEW);
}

//For text output-style 2
void drawText2(const char* text, int length, GLfloat x, GLfloat y){
    glMatrixMode(GL_PROJECTION);
    double *matrix = new double[16];
    glGetDoublev(GL_PROJECTION_MATRIX, matrix);
    glLoadIdentity();
    glOrtho(0, 800, 0, 600, -5, -5);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glPushMatrix();
    glLoadIdentity();
    glRasterPos2f(x, y);
    
    for(int i=0; i<length; i++){
        glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, text[i]);
    }
    glPopMatrix();
    glMatrixMode(GL_PROJECTION);
    glLoadMatrixd(matrix);
    glMatrixMode(GL_MODELVIEW);
}

//Display different notes according to different cases
std::string displayNote(){
    std::string note;
    int state = 0;
    if (fatalSignal){
        state = 2;
    }
    if (HR >= 90){
        state = 1;
    }
    
    switch (state) {
        case 0:
            note = "Congratulations! Your heart is healthy.Keep it up!";
            break;
        case 1:
            note = "Oops! Heart rate is higher than the standard.";
            break;
        case 2:
            note = "Fatal ECG signal! Sudden death may occur.";
            break;
    }
    return note;
}

//Plot ECG waveform
void ecgPlot() {
    int ecgLen;
    //obtian ECG data length due to different sampling rate
    if (fatalSignal == false){
	ecgLen = (int)(dataLen * 1000);
    } else {
	ecgLen = (int)(dataLen * 250);
    }
    int state;
    if (keyHit+4000 > ecgLen) {
        state = 1;
    }else if (keyHit < 0){
        state = 0;
    }else {
        state = 2;
    }
    
    
    if (state == 1){
       keyHit = ecgLen - 4000;
    }else if (state == 0) {
        keyHit = 0;
    }
    
    std::vector<double> dataPlot;
     for (int i = keyHit; i < (4000+keyHit); ++i) {
         dataPlot.push_back(ecgDisplay[i]);
     }
     double max_data = *std::max_element(dataPlot.begin(), dataPlot.end());
     double min_data = *std::min_element(dataPlot.begin(), dataPlot.end());
     for (int i = 0; i < 4000; ++i) {
         dataPlot[i] -= min_data;
     }
     double range = max_data - min_data;
     //draw ecg plot
     int index = 0;
     for(double i=0.04;i<0.96-0.00046;i=i+0.00023)
     {
         double y=0.55;
         double y_=0.55;
         y += (dataPlot[index]/range)*0.40;
         y_ +=(dataPlot[index+1]/range)*0.40;
         glBegin(GL_LINES);
         glVertex2f(i,y);
         glVertex2f(i+0.00023,y_);
         glEnd();
         ++index;
     }
    
    
}


//important display function
void display() {
    glClearColor(1, 0.627451, 0.478431, 0.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    // Draw a quad.
    glBegin(GL_QUADS);
    glColor3f(0.25098, 0.878431, 0.815686);
    glVertex2f(0.08, 0.05);
    glVertex2f(0.23, 0.05);
    glVertex2f(0.23, 0.12);
    glVertex2f(0.08, 0.12);
    glEnd();
    
    // Draw a quad.
    glBegin(GL_QUADS);
    glColor3f(0.25098, 0.878431, 0.815686);
    glVertex2f(0.31, 0.05);
    glVertex2f(0.46, 0.05);
    glVertex2f(0.46, 0.12);
    glVertex2f(0.31, 0.12);
    glEnd();
    
    // Draw a quad.
    glBegin(GL_QUADS);
    glColor3f(0.25098, 0.878431, 0.815686);
    glVertex2f(0.54, 0.05);
    glVertex2f(0.69, 0.05);
    glVertex2f(0.69, 0.12);
    glVertex2f(0.54, 0.12);
    glEnd();
    
    // Draw a quad.
    glBegin(GL_QUADS);
    glColor3f(0.25098, 0.878431, 0.815686);
    glVertex2f(0.77, 0.05);
    glVertex2f(0.92, 0.05);
    glVertex2f(0.92, 0.12);
    glVertex2f(0.77, 0.12);
    glEnd();
    
    // Draw a quad.
    glBegin(GL_QUADS);
    glColor3f(0.564706, 0.933333, 0.564706);
    glVertex2f(0.07, 0.20);
    glVertex2f(0.31, 0.20);
    glVertex2f(0.31, 0.28);
    glVertex2f(0.07, 0.28);
    glEnd();

    // Draw a quad.
    glBegin(GL_QUADS);
    glColor3f(0.564706, 0.933333, 0.564706);
    glVertex2f(0.38, 0.20);
    glVertex2f(0.62, 0.20);
    glVertex2f(0.62, 0.28);
    glVertex2f(0.38, 0.28);
    glEnd();
    
    // Draw a quad.
    glBegin(GL_QUADS);
    glColor3f(0.564706, 0.933333, 0.564706);
    glVertex2f(0.69, 0.20);
    glVertex2f(0.93, 0.20);
    glVertex2f(0.93, 0.28);
    glVertex2f(0.69, 0.28);
    glEnd();
    
    //draw a note space
    if (!fatalSignal){
        glBegin(GL_QUADS);
        glColor3f(0.901961, 0.901961, 0.980392);
        glVertex2f(0.10, 0.34);
        glVertex2f(0.90, 0.34);
        glVertex2f(0.90, 0.52);
        glVertex2f(0.10, 0.52);
        glEnd();
    } else {
        glBegin(GL_QUADS);
        glColor3f(1, 0.0784314, 0.576471);
        glVertex2f(0.10, 0.34);
        glVertex2f(0.90, 0.34);
        glVertex2f(0.90, 0.52);
        glVertex2f(0.10, 0.52);
        glEnd();
    }
    
    //draw note decoration
    glBegin(GL_LINES);
    glColor3f(0.517647, 0.439216, 1);
    glVertex2f(0.095, 0.335);
    glVertex2f(0.905, 0.335);
    glEnd();
    
    glBegin(GL_LINES);
    glColor3f(0.517647, 0.439216, 1);
    glVertex2f(0.905, 0.335);
    glVertex2f(0.905, 0.525);
    glEnd();
    
    glBegin(GL_LINES);
    glColor3f(0.517647, 0.439216, 1);
    glVertex2f(0.095, 0.525);
    glVertex2f(0.905, 0.525);
    glEnd();
    
    glBegin(GL_LINES);
    glColor3f(0.517647, 0.439216, 1);
    glVertex2f(0.095, 0.335);
    glVertex2f(0.095, 0.525);
    glEnd();
    
    //draw ECG plot window
    glBegin(GL_QUADS);
    glColor3f(1.0, 1.0, 1.0);
    glVertex2f(0.04, 0.54);
    glVertex2f(0.96, 0.54);
    glVertex2f(0.96, 0.96);
    glVertex2f(0.04, 0.96);
    glEnd();
    //draw plot window decoration
    glPushMatrix();
    glLineWidth(2);
    glBegin(GL_LINES);
    glColor3f(1, 0.411765, 0.705882);
    glVertex2f(0.035, 0.535);
    glVertex2f(0.965, 0.535);
    glEnd();
    
    glBegin(GL_LINES);
    glColor3f(1, 0.411765, 0.705882);
    glVertex2f(0.965, 0.535);
    glVertex2f(0.965, 0.965);
    glEnd();
    
    glBegin(GL_LINES);
    glColor3f(1, 0.411765, 0.705882);
    glVertex2f(0.035, 0.965);
    glVertex2f(0.965, 0.965);
    glEnd();
    
    glBegin(GL_LINES);
    glColor3f(1, 0.411765, 0.705882);
    glVertex2f(0.035, 0.965);
    glVertex2f(0.035, 0.535);
    glEnd();
    glPopMatrix();
    
    //write text to screen
    std::string text;
    text = "HR(bpm)";
    glColor3f(0.0, 0.0, 0.0);
    drawText(text.data(), text.size(), -0.8, -0.73);
    text = "QRS(ms)";
    drawText(text.data(), text.size(), -0.35, -0.73);
    text = "QTc(ms)";
    drawText(text.data(), text.size(), 0.1, -0.73);
    text = "PR(ms)";
    drawText(text.data(), text.size(), 0.55, -0.73);
    
    text = "Data Name";
    drawText(text.data(), text.size(), -0.75, -0.4);
    text = "Data Length(s)";
    drawText(text.data(), text.size(), -0.25, -0.4);
    text = "Sampling Rate(Hz)";
    drawText(text.data(), text.size(), 0.35, -0.4);
   
    
    //output the parameters value
    glColor3f(0.0, 0.0, 0.0);
    text = std::to_string((int)HR);
    drawText2(text.data(), text.size(), -0.75, -0.86);
    text = std::to_string(QRS);
    drawText2(text.data(), 5, -0.33, -0.86);
    text = std::to_string(QTc);
    drawText2(text.data(), 5, 0.14, -0.86);
    text = std::to_string(PR);
    drawText2(text.data(), 5, 0.60, -0.86);
    
    //display note
    glPushMatrix();
    glPointSize(2);
    text = displayNote();
    drawText(text.data(), text.size(), -0.75, -0.1);
    glPopMatrix();
    
    //output the information of ECG data
    text = std::to_string(dataLen);
    drawText2(text.data(), 5, -0.1, -0.56);
    text = std::to_string(samplingRate);
    drawText2(text.data(), text.size(), 0.52, -0.56);
    //draw data name according to pop-up menu
    drawText2(dataName.data(), dataName.size(), -0.82, -0.56);
    
    glPointSize(1);
    glColor3f(0.0,0.0,0.0);
    
   //ecg waveform plot
    glEnable (GL_LINE_SMOOTH);
    ecgPlot();
    glDisable(GL_LINE_SMOOTH);
   
    glutSwapBuffers();
}


//reshape function
void reshape(GLint width, GLint height)
{
    // Compute the aspect ratio, avoiding the possibility of division by zero.
    GLfloat aspectRatio = static_cast<GLfloat>(width) /
    ((height) ? height : 1.0);
    
    // Set the viewport to the entire window.
    glViewport(0, 0, width, height);
    
    // Initialize the projection matrix.
    // This is done in such a way to maintain the aspect ratio in case
    // the window shape is not square.
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    if (width >= height) {
        gluOrtho2D(0.0, 1.0 * aspectRatio, 0.0, 1.0);
    } else {
        gluOrtho2D(0.0, 1.0, 0.0, 1.0 / aspectRatio);
    }
    // Initialize the modelview matrix.
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

//key board function
void keyboard(unsigned char key, int x, int y){
    switch (key) {
    case 'q':
            exit(0);
            break;
    }
}
//special keys function  
void processSpecialKeys(int key, int x, int y){
    switch (key){
        case GLUT_KEY_LEFT:
            keyHit -= 2000;
            glutPostRedisplay();
            break;
        case GLUT_KEY_RIGHT:
            keyHit += 2000;
            glutPostRedisplay();
            break;
    }
    
}

//calls for ECGDetect class
void computePara(){
    ECGDetect ecg(dataName, 1000);
    ecg.run();
    HR = ecg.get_HR();
    QRS = ecg.get_QRS_duration();
    PR = ecg.get_PR_interval();
    QTc = ecg.get_QTc();
    samplingRate = ecg.getFS();
    dataLen = ((double)ecg.get_ecg_data_length())/((double)samplingRate);
    ecgDisplay = ecg.get_buffer_plot();
    keyHit = 0;
    fatalSignal = false;
    
}
//Calls for FatalDetect class
void fatalPara() {
    FatalDetect fatal(dataName, 250);
    fatal.run();
    fatalSignal = fatal.fatal();
    HR = NULL;
    QRS = NULL;
    PR = NULL;
    QTc = NULL;
    samplingRate = fatal.getFS();
    dataLen = ((double)fatal.get_ecg_data_length())/(double)samplingRate;
    ecgDisplay.clear();
    ecgDisplay = fatal.get_buffer_plot();
    keyHit = 0;

}
//pop-up menu design
void dataNameMenu (int id){
   
    switch (id){
        case 1:
            dataName = "ecg1.bin"; //1-6 normal ECG data
            computePara();
            glutPostRedisplay();
            break;
        case 2:
            dataName = "ecg2.bin";
            computePara();
            glutPostRedisplay();
            break;
        case 3:
            dataName = "ecg3.bin";
            computePara();
            glutPostRedisplay();
            break;
        case 4:
            dataName = "ecg4.bin";
            computePara();
            glutPostRedisplay();
            break;
        case 5:
            dataName = "ecg5.bin";
            computePara();
            glutPostRedisplay();
            break;
        case 6:
            dataName = "ecg6.bin";
            computePara();
            glutPostRedisplay();
            break;
	 case 7:
            dataName = "ecg7.txt"; //fatal signal data
            fatalPara();
            glutPostRedisplay();
            break;
        case 8:
            dataName = "ecg8.txt"; //fatal signal data
            fatalPara();
            glutPostRedisplay();
            break; 
	case 9:
            dataName = "ecg_n1.bin"; //data with noise
            computePara();
            glutPostRedisplay();
            break;
	
	case 10:
            dataName = "ecg_n2.bin"; //data with noise
            computePara();
            glutPostRedisplay();
            break;
        
    }
    
}

int main (int argc, char** argv) {
    
    computePara();
    const int winWidth = 512; // The nominal window width.
    const int winHeight = 512; // The nominal window height.
    int data_name_menu;

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(winWidth, winHeight);
    glutCreateWindow("ECG Detection");
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);
    glutSpecialFunc(processSpecialKeys);
    
    data_name_menu = glutCreateMenu(dataNameMenu);
        glutAddMenuEntry("ecg1.bin", 1);
        glutAddMenuEntry("ecg2.bin", 2);
        glutAddMenuEntry("ecg3.bin", 3);
        glutAddMenuEntry("ecg4.bin", 4);
        glutAddMenuEntry("ecg5.bin", 5);
        glutAddMenuEntry("ecg6.bin", 6);
	glutAddMenuEntry("ecg7.txt", 7);
        glutAddMenuEntry("ecg8.txt", 8);
	glutAddMenuEntry("ecg_n1.bin", 9);
	glutAddMenuEntry("ecg_n2.bin", 10);
	
    glutAttachMenu(GLUT_RIGHT_BUTTON);
    
    glutMainLoop();
    
    return 0;
}
