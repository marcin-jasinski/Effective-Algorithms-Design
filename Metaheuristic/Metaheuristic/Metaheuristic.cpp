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

#include "SimulatedAnnealing.h"
#include "TabuSearch.h"

const double PI = 3.141592;
const double RRR = 6378.388;

int citiesNumber;

static int nint(double);
static double convert_to_geo(double);
int TWOD_geo_distance(const Point a, const Point b);

int** readAsymetricTSPData(std::string);
int** readSymetricTSPData(std::string);

int main()
{
	std::srand(unsigned(std::time(0)));

	int** asymetricEdgesMatrix = readAsymetricTSPData("ATSP_Data/ftv64.atsp");

	// int** symetricEdgesMatrix = readSymetricTSPData("TSP_Data/burma14.tsp");

	/*for (double maxTemperature = 250; maxTemperature >= 10; maxTemperature -= 10)
	{
		std::cout << "\nStarting temperature: " << maxTemperature << std::endl;
		double tempChange = 0.9999;
		std::vector<int> bestSolution = simulatedAnnealing(maxTemperature, tempChange, citiesNumber, symetricEdgesMatrix);
		std::cout << "Best solution: " << std::endl;
		for (int i = 0; i < citiesNumber; i++)
		{
			std::cout << std::setw(4) << bestSolution.at(i);
		}
		std::cout << std::setw(4) << bestSolution.at(0);
		int bestCost = getRouteCost(bestSolution, symetricEdgesMatrix, citiesNumber);
		std::cout << "\nBest cost: " << bestCost << std::endl;
	}*/

	for (int i = 1; i < 2; i++)
	{
		std::vector<int> bestSolution = tabuSearch(citiesNumber, asymetricEdgesMatrix);
		std::cout << "Best solution: " << std::endl;
		for (int i = 0; i < citiesNumber; i++)
		{
			std::cout << std::setw(4) << bestSolution.at(i);
		}
		std::cout << std::setw(4) << bestSolution.at(0);
		int bestCost = getRouteCost(bestSolution, asymetricEdgesMatrix, citiesNumber);
		std::cout << "\nBest cost: " << bestCost << std::endl;
	}
}

int** readAsymetricTSPData(std::string fileName)
{
	if (fileName.find(".atsp") == std::string::npos) {
		std::cout << "Uncorrect file format." << '\n';
		return nullptr;
	}

	std::fstream sourceDataFile;
	sourceDataFile.open(fileName, std::ios::in);

	int **edgesMatrix = nullptr;

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

	int **edgesMatrix = nullptr;

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
		double x, y;
		
		do
		{
			sourceDataFile >> readLine;
		} while (readLine.find("NODE_COORD_SECTION") == std::string::npos);

		std::vector<Point> pointsVector;

		for (int i = 0; i < citiesNumber; i++)
		{
			sourceDataFile >> currentCityIndex >> x >> y;
			pointsVector.push_back(Point(x, y));
		}

		for (int i = 0; i < citiesNumber; i++)
		{
			Point currentCity = pointsVector.at(i);

			for (int j = 0; j < citiesNumber; j++)
			{
				Point nextPoint = pointsVector.at(j);
				if (i == j) edgesMatrix[i][j] = 0;
				else edgesMatrix[i][j] = TWOD_geo_distance(currentCity, nextPoint);
			}
		}
	}

	return edgesMatrix;
}

int nint(double d)
{
	return std::floor(d + 0.5);
}

double convert_to_geo(double x) 
{
	int deg = nint(x);
	return PI * (deg + 5.0 * (x - deg) / 3.0) / 180.0;
}

int TWOD_geo_distance(const Point a, const Point b)
{
	Point a_geo( convert_to_geo(a.x), convert_to_geo(a.y));
	Point b_geo( convert_to_geo(b.x), convert_to_geo(b.y));

	double q1 = std::cos(a_geo.y - b_geo.y);
	double q2 = std::cos(a_geo.x - b_geo.x);
	double q3 = std::cos(a_geo.x + b_geo.x);

	return (int)RRR * std::acos(0.5 * ((1.0 + q1)*q2 - (1.0 - q1) * q3)) + 1.0;
}