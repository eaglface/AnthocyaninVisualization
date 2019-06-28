#include "ofxHyperSpectralImage.h"

#include <iostream>
#include <fstream>
#include <opencv2\opencv.hpp>

static bool split(const string& str, char split_on, string& left, string& right)
{
    int pos = str.find(split_on, 0);
    if (pos == string::npos)
        return false;

    left = str.substr(0, pos);
    right = str.substr(pos + 1);
    return true;
}

static void trim(string& str, const string& chars_to_trim)
{
    while (chars_to_trim.find(str[0]) != string::npos)
        str = str.substr(1);
    while (chars_to_trim.find(str[str.size() - 1]) != string::npos)
        str.resize(str.size() - 1);
}


hyperSpectralImage::hyperSpectralImage()
{

}
hyperSpectralImage::~hyperSpectralImage()
{

}

int hyperSpectralImage::loadColorMatchingFunctionTable(string path)
{
    colorMatchingFunction.clear();

    bool success = false;
    ifstream file;
    // Try opening the file
    file.open(path);
    if (file.is_open())
    {
        cout << "Loading color matching function table..."; cout.flush();
        // Read the entire file, one line at a time
        while (file)
        {
            string line;
            getline(file, line);

            vector <string> cmfLine;
            cmfLine = ofSplitString(line, ",", true, true);
            if (cmfLine.size() == 4) {
                ofVec4f temp;
                temp.x = atof(cmfLine[0].c_str());
                temp.y = atof(cmfLine[1].c_str());
                temp.z = atof(cmfLine[2].c_str());
                temp.w = atof(cmfLine[3].c_str());

                colorMatchingFunction.push_back(temp);
            }
        }
        cout << "done. size = " << colorMatchingFunction.size() << endl;
        file.close();
        return 1;

    }
    else
        return 0;
    
}

int hyperSpectralImage::loadColorMapTable(string path)
{
    colorMapTable.clear();

    bool success = false;
    ifstream file;
    // Try opening the file
    file.open(path);
    if (file.is_open())
    {
        cout << "Loading color map function table..."; cout.flush();
        // Read the entire file, one line at a time
        while (file)
        {
            string line;
            getline(file, line);

            vector <string> cmfLine;
            cmfLine = ofSplitString(line, ",", true, true);
            if (cmfLine.size() == 4) {
                ofVec4f temp;
                temp.x = atof(cmfLine[0].c_str());
                temp.y = atof(cmfLine[1].c_str());
                temp.z = atof(cmfLine[2].c_str());
                temp.w = atof(cmfLine[3].c_str());

                colorMapTable.push_back(temp);
            }
        }
        cout << "done. size = " << colorMapTable.size() << endl;
        file.close();
        return 1;

    }
    else
        return 0;

}


int hyperSpectralImage::buildColorMatchTable()
{
    //colorMatchTable[3];
    colorMatchTable[0].resize(wavelengths.size());
    colorMatchTable[1].resize(wavelengths.size());
    colorMatchTable[2].resize(wavelengths.size());

    for (int ii = 0; ii < wavelengths.size(); ii++){
        double currentWavelength = wavelengths[ii];
        double minError = 5;
        int selectedIndex = -1;
        for (int jj = 0; jj < colorMatchingFunction.size(); jj++){
            double error = fabs(colorMatchingFunction[jj].x - currentWavelength);
            if (minError>error) {
                minError = error;
                selectedIndex = jj;
            }
        }
        if (selectedIndex != -1) {
            colorMatchTable[0][ii] = colorMatchingFunction[selectedIndex].y;
            colorMatchTable[1][ii] = colorMatchingFunction[selectedIndex].z;
            colorMatchTable[2][ii] = colorMatchingFunction[selectedIndex].w;
        }
        else {
            colorMatchTable[0][ii] = 0;
            colorMatchTable[1][ii] = 0;
            colorMatchTable[2][ii] = 0;
        }
    }

    return 1;

}


int hyperSpectralImage::headerLoad(string path)
{
    bool success = false;
    ifstream file;
    // Try opening the file
    file.open(path);
    if (file.is_open())
    {
        // Read the entire file, one line at a time
        while (file)
        {
            string line, keyword, value;

            // Read a line
            getline(file, line);

            // Split on '=', ignore if there wasn't a '='
            if (split(line, '=', keyword, value))
            {
                trim(keyword, " \t\n");
                trim(value, " \t\n");

                if (keyword == "description")
                    description = value;
                else if (keyword == "samples")
                    samples = atoi(value.c_str());
                else if (keyword == "lines")
                    lines = atoi(value.c_str());
                else if (keyword == "bands")
                    bands = atoi(value.c_str());
                else if (keyword == "header offset")
                    header_offset = atoi(value.c_str());
                else if (keyword == "file type")
                    file_type = value;
                else if (keyword == "data type")
                    data_type = (DType)atoi(value.c_str());
                else if (keyword == "interleave")
                {
                    if (value == "bsq")
                        interleave = it_BSQ;
                    else if (value == "bip")
                        interleave = it_BIP;
                    else if (value == "bil")
                        interleave = it_BIL;
                }
                else if (keyword == "byte order")
                    big_endian = atoi(value.c_str()) == 1;
                else if (keyword == "x-start")
                    x_start = atoi(value.c_str());
                else if (keyword == "y-start")
                    y_start = atoi(value.c_str());
                else if (keyword == "map info")
                    misc += value;
                else if (keyword == "projection info")
                    misc = value;
                else if (keyword == "band names")
                    band_names = value;
                else if (keyword == "wavelength units")
                    wavelength = value.c_str();
                else if (keyword == "wavelength") {
                    bool inParenthesis = true;
                    string line1;
                    wavelengths.clear();
                    while (inParenthesis) {
                        getline(file, line1);
                        // check if the line1 has }
                        int index = line1.find("}");
                        if (line1.npos == index) {  // no ending parenthesis
                            
                        }
                        else {// this is end of the parenthesis
                            line1 = line1.substr(0, index - 1); 
                            inParenthesis = false;
                        }

                        // tokenizing
                        vector <string> wavelengthString;
                        wavelengthString = ofSplitString(line1, ",", true, true);
                        for (int ii = 0; ii < wavelengthString.size(); ii++){
                            wavelengths.push_back(atof(wavelengthString[ii].c_str()));
                        }

                    }

                }
            }
        }

		//for (int ii = 0; ii < wavelengths.size(); ii++)
		//	cout << ii << " -- " << wavelengths[ii] << endl;
  //      cout << endl;

        success = true;
    }

    return success;
}

int hyperSpectralImage::imageLoad(string path)
{
	rgbImage = cv::imread(path);
	return 0;
}

int hyperSpectralImage::load(string path)
{
    height = lines;// 696;
    width = samples;// 520;
   // bands = 128;
    headerOffset = header_offset;//  165888;// 32768; //

    images.resize(bands);
    for (int ii = 0; ii < bands; ii++){
        images[ii].create(width, height, CV_16UC1);
    }

    FILE *fp = fopen(path.c_str(), "rb");
    if (fp == NULL) {
        cout << "File " << path << " read error. quit()" << endl;
        return 0;
    }


    fseek(fp, 0, SEEK_END);
    unsigned int fileLen = ftell(fp);
    fseek(fp, 0, SEEK_SET);

    cout << "file size = " << fileLen << endl;

    cout << "opening cube file" << endl;

    char *buf = new char[fileLen + 1];
    uint16_t *buffer;

    fread(buf, fileLen, 1, fp);

    fclose(fp);

    cout << "Buffer read success.." << fileLen << endl;

    // bil read
    buffer = (uint16_t*)&buf[headerOffset];
    int blocksize = width * 2;
    int rowsize = width * 2 * bands;
    for (int yy = 0; yy < height; yy++){
        for (int channel = 0; channel < bands; channel++){
            for (int xx = 0; xx < width; xx++) {
                images[channel].at<uint16_t>(width-xx-1, height-1-yy) =
                    // (uint16_t)(buf[headerOffset + yy*rowsize + blocksize*channel + xx*2]);
                    buffer[yy*width*bands + width*channel + xx];
            }
            //void *pt = images[channel].ptr<uint16_t>(yy);
            //memcpy(pt, &buf[headerOffset + yy*rowsize + blocksize*channel], blocksize);
        }
    }

    cout << "load done.." << fileLen << endl;

    delete buf;
}


void hyperSpectralImage::convertLinearXYZ()
{

    double minx = 10000, maxx = 0;
    double miny = 10000, maxy = 0;
    double minz = 10000, maxz = 0;

    tempImage.create(width, height, CV_32FC3);
    displayImage.create(width, height, CV_8UC3);
    for (int yy = 0; yy < height; yy++)
    for (int xx = 0; xx < width; xx++) {
        double x = 0, y = 0, z = 0;
        for (int ch = 0; ch < bands; ch++) {
            x += images[ch].at<uint16_t>(xx, yy)*colorMatchTable[0][ch];
            y += images[ch].at<uint16_t>(xx, yy)*colorMatchTable[1][ch];
            z += images[ch].at<uint16_t>(xx, yy)*colorMatchTable[2][ch];
        }

        minx = min(minx, x); maxx = max(maxx, x);
        miny = min(miny, y); maxy = max(maxy, y);
        minz = min(minz, z); maxz = max(maxz, z);

        tempImage.at<Vec3f>(xx, yy)[0] = x / 100000.0;
        tempImage.at<Vec3f>(xx, yy)[1] = y / 100000.0;
        tempImage.at<Vec3f>(xx, yy)[2] = z / 100000.0;
    }

    //cout << "minx = " << minx << " ,maxx = " << maxx << endl;
    //cout << "miny = " << miny << " ,maxy = " << maxy << endl;
    //cout << "minz = " << minz << " ,maxz = " << maxz << endl;

    // think about the scaling

    cvtColor(tempImage, displayImage, COLOR_XYZ2BGR);
    imshow("rgb", displayImage);
}

int hyperSpectralImage::thresholdAndRegionCount(float thres, cv::Mat& inputImage, cv::Mat& thresholdImage, bool countHigh)
{
	//thresholdImage.create(inputImage.size(), CV_8UC1);
	cv::Mat tempImage;
	//tempImage.create(inputImage.size(), CV_8UC1);

	cvtColor(inputImage, tempImage, CV_BGR2GRAY);
	//cv::Mat temp1;
	//tempImage.convertTo(temp1, CV_8UC1, 255.0);
	cv::threshold(tempImage, thresholdImage, thres, 255, THRESH_BINARY);
	//cv::adaptiveThreshold(temp1, thresholdImage, 255, ADAPTIVE_THRESH_MEAN_C, THRESH_BINARY, 51,10);//)
	//floodFill(thresholdImage, Point(tempImage.cols / 2, tempImage.rows / 2), 0.5);
	imshow("gray", tempImage);
	imshow("thres", thresholdImage);

	// count numbers of the white
	int count = 0;
	float maxI = 0, minI = 255;
	for (int yy = 0; yy < thresholdImage.rows; yy++)
	for (int xx = 0; xx < thresholdImage.cols; xx++) {
		//printf("%d ", tempImage.at<unsigned char>(yy, xx));;// << " ";
		maxI = maxI<tempImage.at<float>(yy, xx) ? tempImage.at<float>(yy, xx) : maxI;
		minI = minI>tempImage.at<float>(yy, xx) ? tempImage.at<float>(yy, xx) : minI;

		//if (thresholdImage.at<float>(yy, xx) == 0.5)
		//	count++;
		//else
		//	thresholdImage.at<float>(yy, xx) = 0;

		if (countHigh) {
		//	if (thresholdImage.at<unsigned char>(yy, xx)>128)
			if (tempImage.at<unsigned char>(yy, xx)>thres)
				count++;
		}
		else  {
			if (tempImage.at<unsigned char>(yy, xx) < thres)
				count++;
		}
	}

	cout << "min Intensity = " << minI << ", max Intensity = " << maxI << endl; 
	return count;
}


void hyperSpectralImage::correlation1D(vector <float> *corrCoef)
{
    double minx = 10000, maxx = 0;
    cout << "correlation1D" << endl;
    corrImage.create(width, height, CV_32FC1);

    for (int yy = 0; yy < height; yy++)
    for (int xx = 0; xx < width; xx++) {
        double x = 0;
        double mindata = 10000, maxdata = 0;
        for (int ch = 0; ch < bands; ch++) {
            double dd = images[ch].at<uint16_t>(xx, yy);
            mindata = min(dd, mindata); maxdata = max(dd, maxdata);
        }
        for (int ch = 0; ch < bands; ch++) {
            x += images[ch].at<uint16_t>(xx, yy)/(maxdata-mindata)*corrCoef->at(ch);
        }

        minx = min(minx, x); maxx = max(maxx, x);
        corrImage.at<float>(xx, yy) = x;
    }

    minCorrValue = minx;
    maxCorrValue = maxx;
    cout << "min = " << minx << ", max = " << maxx << endl;
}

void hyperSpectralImage::makeColoredImagefromCorrelation()
{
    cout << "makeColoredImagefromCorrelation" << endl;

    vector <int> votes;
    votes.resize(colorMapTable.size());
    for (int ii = 0; ii < votes.size(); ii++)
        votes[ii] = 0;

    // normalize the data into [0 1]
    // think about the normalization strategy
    displayImage.create(width, height, CV_8UC3);

    if (minCorrValue < 0) {
        for (int yy = 0; yy < height; yy++)
        for (int xx = 0; xx < width; xx++) {
            corrImage.at<float>(xx, yy) = corrImage.at<float>(xx, yy) + minCorrValue;
        }
        maxCorrValue += minCorrValue;
    }

    for (int yy = 0; yy < height; yy++)
    for (int xx = 0; xx < width; xx++) {
        float normalizedValue;
        normalizedValue = corrImage.at<float>(xx, yy) / maxCorrValue;

        // find index
        int index = normalizedValue*(colorMapTable.size()-1);
        displayImage.at<Vec3b>(xx, yy)[0] = colorMapTable[index].w*255.0;
        displayImage.at<Vec3b>(xx, yy)[1] = colorMapTable[index].z*255.0;
        displayImage.at<Vec3b>(xx, yy)[2] = colorMapTable[index].y*255.0;
        votes[index]++;
    }


    imshow("responseMap", displayImage);

    // generate colormap legend
    colormapLegendImage.create(colorMapTable.size(), 130, CV_8UC3);
    for (int yy = 0; yy < 100; yy++)
    for (int xx = 0; xx < colorMapTable.size(); xx++) {
        colormapLegendImage.at<Vec3b>(xx, yy)[0] = colorMapTable[xx].w*255.0;
        colormapLegendImage.at<Vec3b>(xx, yy)[1] = colorMapTable[xx].z*255.0;
        colormapLegendImage.at<Vec3b>(xx, yy)[2] = colorMapTable[xx].y*255.0;
    }

    imshow("Color Map Legend", colormapLegendImage);

}

void hyperSpectralImage::saveReferenceData()
{
    if (referenceData.size() == 0) {
        cout << " No reference data. quit." << endl;
        return;
    }
    else cout << "referenceData size = " << referenceData.size() << endl;
    FILE *fp = fopen("referenceData.txt", "wt");
    if (fp != 0) {
        for (int ii = 0; ii < referenceData.size(); ii++)
            fprintf(fp, "%ld\n", (uint16_t)referenceData[ii]);

        fclose(fp);
        cout << "Reference Data saved in \"referenceData.txt \" " << endl;
    }
}

void hyperSpectralImage::loadReferenceData()
{
    uint16_t data;
    FILE *fp = fopen("referenceData.txt", "rt");
    
    if (fp != 0) {
        referenceData.clear();
        while (fscanf(fp, "%ld", &data) != EOF)
            referenceData.push_back((double)data);
        fclose(fp);
        cout << "Reference Data loaded from \"referenceData.txt \" : data number = " << referenceData.size() << endl;
    }
    else {
        cout << "No reference data file" << endl;
    }
}

float hyperSpectralImage::calcAbsorbanceAnthocyanin()
{
    double minx = 10000, maxx = 0;
    cout << "calcAbsorbanceAnthocyanin" << endl;

    // memory allocation
    absorbanceAnthoCyanin.create(width, height, CV_32FC1);

    //int index = 29; // 519.4 nm
	int index = 35; // 550.14 nm

	double referenceValue = 2773; // referenceData[index];
    cout << " referenceValue = " << referenceValue << endl;

    for (int yy = 0; yy < height; yy++)
    for (int xx = 0; xx < width; xx++) {

        double dd = images[index].at<uint16_t>(xx, yy);
        double x = -log((double)dd / (double)referenceValue);

        minx = min(minx, x); maxx = max(maxx, x);
        absorbanceAnthoCyanin.at<float>(xx, yy) = x;
    }

    cout << "min = " << minx << ", max = " << maxx << endl;
	return maxx;
}


void hyperSpectralImage::makeColoredImagefromAnthocyanin(double minValue, double maxValue)
{
    cout << "makeColoredImagefromAnthocyanin" << endl;

    // normalize the data into [0 1]
    // think about the normalization strategy
    displayImage.create(width, height, CV_8UC3);

    for (int yy = 0; yy < height; yy++)
    for (int xx = 0; xx < width; xx++) {
        float normalizedValue;
        //normalizedValue = (absorbanceAnthoCyanin.at<float>(xx, yy)) / (maxValue);
		normalizedValue = (absorbanceAnthoCyanin.at<float>(xx, yy) - minValue) / (maxValue - minValue);
        if (normalizedValue < 0) normalizedValue = 0;

        // find index
        int index = normalizedValue*(colorMapTable.size() - 1);
        displayImage.at<Vec3b>(xx, yy)[0] = colorMapTable[index].w*255.0;
        displayImage.at<Vec3b>(xx, yy)[1] = colorMapTable[index].z*255.0;
        displayImage.at<Vec3b>(xx, yy)[2] = colorMapTable[index].y*255.0;
    }


    imshow("Anthocyanin Response", displayImage);

    // generate colormap legend
    colormapLegendImage.create(70, colorMapTable.size(), CV_8UC3);
    colormapLegendImage.setTo(255);
    for (int yy = 0; yy < 50; yy++)
    for (int xx = 0; xx < colorMapTable.size(); xx++) {
        colormapLegendImage.at<Vec3b>(yy, xx)[0] = colorMapTable[xx].w*255.0;
        colormapLegendImage.at<Vec3b>(yy, xx)[1] = colorMapTable[xx].z*255.0;
        colormapLegendImage.at<Vec3b>(yy, xx)[2] = colorMapTable[xx].y*255.0;
    }
    cv::Point textOrg(0, 65);
    cv::Point textOrg1(colorMapTable.size()-25, 65);
    cv::Point textOrg2(colorMapTable.size() / 2-7, 65);

    int fontFace = FONT_HERSHEY_PLAIN;

	char midValueString[100], maxValueString[100];
	sprintf(midValueString, "%1.1f", maxValue / 2.0);
	sprintf(maxValueString, "%1.1f", maxValue);

    cv::putText(colormapLegendImage, "0", textOrg, fontFace, 1, Scalar::all(0),1);
	cv::putText(colormapLegendImage, midValueString, textOrg2, fontFace, 1, Scalar::all(0), 1);
	cv::putText(colormapLegendImage, maxValueString, textOrg1, fontFace, 1, Scalar::all(0), 1);

    imshow("Color Map Legend", colormapLegendImage);
}
