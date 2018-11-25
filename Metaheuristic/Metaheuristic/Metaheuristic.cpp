#include "pch.h"
#include "Point.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <windows.h>
#include <vector>

const double PI = 3.141592;
const double RRR = 6378.388;

int citiesNumber;
int** edgesMatrix;

static int nint(double);
static double convert_to_geo(double);
int TWOD_geo_distance(const Point a, const Point b);
int** readAsymetricTSPData(std::string);
int** readSymetricTSPData(std::string);

int main()
{
	readAsymetricTSPData("ATSP_Data/br17.atsp");

	readSymetricTSPData("TSP_Data/burma14.tsp");
}

int** readAsymetricTSPData(std::string fileName)
{
	if (fileName.find(".atsp") == std::string::npos) {
		std::cout << "Uncorrect file format." << '\n';
		return nullptr;
	}

	std::fstream sourceDataFile;
	sourceDataFile.open(fileName, std::ios::in);

	int citiesNumber = 0;
	int** edgesMatrix = nullptr;

	if (sourceDataFile.good())
	{
		std::cout << "Succesfully opened a file, begin reading data process..." << std::endl;

		int citiesNumber;

		std::string readLine;
		do
		{
			sourceDataFile >> readLine;
		} while (readLine.find("DIMENSION:") == std::string::npos);

		std::cout << readLine;
		sourceDataFile >> citiesNumber;
		std::cout << "Number of cities: " << citiesNumber << std::endl;

		for (int i = 1; i <= 5; i++) sourceDataFile >> readLine;

		edgesMatrix = new int*[citiesNumber];
		for (int i = 0; i < citiesNumber; i++) edgesMatrix[i] = new int[citiesNumber];
		for (int i = 0; i < citiesNumber; i++)
		{
			for (int j = 0; j < citiesNumber; j++)
			{
				edgesMatrix[i][j] = -1;
			}
		}

		int edgeCost = 0;
		for (int i = 0; i < citiesNumber; i++)
		{
			for (int j = 0; j < citiesNumber; j++)
			{
				sourceDataFile >> edgeCost;
				edgesMatrix[i][j] = edgeCost;
			}
		}
	}

	return edgesMatrix;
}

int** readSymetricTSPData(std::string fileName)
{
	if (fileName.find(".tsp") == std::string::npos) {
		std::cout << "Uncorrect file format." << '\n';
		return nullptr;
	}

	std::fstream sourceDataFile;
	sourceDataFile.open(fileName, std::ios::in);

	int citiesNumber = 0;
	int** edgesMatrix = nullptr;

	if (sourceDataFile.good())
	{

		std::cout << "Succesfully opened a file, begin reading data process..." << std::endl;

		std::string readLine;
		do
		{
			sourceDataFile >> readLine;
		} while (readLine.find("DIMENSION:") == std::string::npos);

		std::cout << readLine;
		sourceDataFile >> citiesNumber;
		std::cout << "Number of cities: " << citiesNumber << std::endl;

		edgesMatrix = new int*[citiesNumber];
		for (int i = 0; i < citiesNumber; i++) edgesMatrix[i] = new int[citiesNumber];
		for (int i = 0; i < citiesNumber; i++)
		{
			for (int j = 0; j < citiesNumber; j++) edgesMatrix[i][j] = -1;
		}

		int currentCityIndex;
		double latitude;
		double longitude;

		for (int i = 1; i <= 7; i++)  sourceDataFile >> readLine;

		std::vector<Point> pointsVector;

		for (int i = 0; i < citiesNumber; i++)
		{
			sourceDataFile >> currentCityIndex >> latitude >> longitude;
			pointsVector.push_back(Point(longitude, latitude));
		}

		for (int i = 0; i < citiesNumber; i++)
		{
			Point currentCity = pointsVector.at(i);

			for (int j = 0; j < citiesNumber; j++)
			{
				Point pointDupa = pointsVector.at(j);
				if (i == j) edgesMatrix[i][j] = 0;
				else edgesMatrix[i][j] = TWOD_geo_distance(currentCity, pointDupa);
			}
		}

		std::cout << "\nRead edges Matrix: " << std::endl;
		for (int i = 0; i < citiesNumber; i++)
		{
			for (int j = 0; j < citiesNumber; j++)
			{
				std::cout << std::setw(5) << edgesMatrix[i][j];
			}
			std::cout << std::endl;
		}
	}

	return edgesMatrix;
}

int nint(double d)
{
	return std::floor(d + 0.5);
}

double convert_to_geo(double x) {

	int deg = nint(x);
	return PI * (deg + 5.0 * (x - deg) / 3.0) / 180.0;
}

int TWOD_geo_distance(const Point a, const Point b) {

	double dupa = a.x;

	Point a_geo( convert_to_geo(a.x), convert_to_geo(a.y));
	Point b_geo( convert_to_geo(b.x), convert_to_geo(b.y));

	double q1 = std::cos(a_geo.y - b_geo.y);
	double q2 = std::cos(a_geo.x - b_geo.x);
	double q3 = std::cos(a_geo.x + b_geo.x);
	return (int)RRR * std::acos(0.5 * ((1.0 + q1)*q2 - (1.0 - q1) * q3)) + 1.0;
}