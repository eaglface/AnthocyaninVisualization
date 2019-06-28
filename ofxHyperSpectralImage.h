#pragma once

#include "ofxOpenCv.h"

using namespace cv;

enum DType
{
    dt_byte8 = 1,
    dt_int16 = 2,
    dt_int32 = 3,
    dt_float32 = 4,
    dt_float64 = 5,
    dt_complex_2x32 = 6,
    dt_complex_2x64 = 9,
    dt_uint16 = 12,
    dt_uint32 = 13,
    dt_uint64 = 14,
    dt_uint64_long = 15  // what's the difference between this ant dt_unint64??
};

enum InterleaveType
{
    it_BSQ,
    it_BIP,
    it_BIL
};

class hyperSpectralImage{

public:
    hyperSpectralImage();
    ~hyperSpectralImage();

    int load(string path);
    int headerLoad(string path);
	int imageLoad(string path);

    int loadColorMatchingFunctionTable(string path);
    int loadColorMapTable(string path);

    int buildColorMatchTable();

    void convertLinearXYZ();
    void correlation1D(vector <float> *corrCoef);
    void makeColoredImagefromCorrelation();

	int thresholdAndRegionCount(float thres, cv::Mat& inputImage, cv::Mat& thresholdImage, bool countHigh = 1);


    void saveReferenceData();
    void loadReferenceData();

public:
    vector<Mat> images;
	cv::Mat rgbImage;
    int width;
    int height;
    int bands;
    int headerOffset;

    vector <double> wavelengths;
    vector <float> colorMatchTable[3];

    vector <ofVec4f> colorMatchingFunction;
    vector <ofVec4f> colorMapTable;

    std::string description;
    int samples;
    int lines;
  //  int bands;
    int header_offset;
    std::string file_type;
    DType data_type;
    InterleaveType interleave;
    bool big_endian;
    int x_start;
    int y_start;

    // Unhandled information, but still read in
    std::string band_names;
    std::string wavelength;
    std::string misc;  // everything else that is unrecognized

    Mat displayImage;
    Mat tempImage;
    Mat corrImage;
    double minCorrValue, maxCorrValue;


    // anthocyanin analysis
    vector <double> referenceData;
    Mat absorbanceAnthoCyanin;
    Mat colormapLegendImage;

    float calcAbsorbanceAnthocyanin();
    void makeColoredImagefromAnthocyanin(double minValue, double maxValue);


};

