#include "pch.h"
#include "Point.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <string>
#include <random>

int** edgesMatrix;
int citiesNumber;
int correctSolution;

const double PI = 3.141592;
const double RRR = 6378.388;

std::vector<int> solution;

std::vector<int> getNextSpecimen(int);
std::vector<std::vector<int>> getStartingPopulation(int, int);

int main()
{
	srand(time(NULL));
    std::cout << "Hello World!\n"; 

	std::getchar();
	return 0;
}

std::vector<int> getNextSpecimen(int problemSize)
{
	int genome;
	std::vector<int> specimen;

	for (int i = 1; i <= problemSize; i++)
	{
		do
		{
			genome = std::rand() % problemSize + 1;
		} while (std::find(specimen.begin(), specimen.end(), genome) != specimen.end());

		specimen.push_back(genome);
	}
	
	return specimen;
}

std::vector<std::vector<int>> getStartingPopulation(int populationSize, int problemSize)
{
	std::vector<std::vector<int>> population;

	for (int i = 1; i <= populationSize; i++)
	{
		std::vector<int> newSpecimen = getNextSpecimen(problemSize);
		population.push_back(newSpecimen);
	}

	return population;
}

void printPath(std::vector<int> solution)
{
	for (int i = 0; i < solution.size(); i++)
	{
		std::cout << solution.at(i) << " -> ";
	}

	std::cout << solution.at(0) << std::endl;
}

int** readAsymetricTSPData(std::string fileName)
{
	if (fileName.find(".atsp") == std::string::npos) {
		std::cout << "Uncorrect file format." << '\n';
		return nullptr;
	}

	std::fstream sourceDataFile;
	sourceDataFile.open(fileName, std::ios::in);

	edgesMatrix = nullptr;

	if (sourceDataFile.good())
	{
		std::string readLine;
		do
		{
			sourceDataFile >> readLine;
		} while (readLine.find("DIMENSION:") == std::string::npos);

		sourceDataFile >> citiesNumber;

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

	edgesMatrix = nullptr;

	if (sourceDataFile.good())
	{
		std::string readLine;
		do
		{
			sourceDataFile >> readLine;
		} while (readLine.find("DIMENSION:") == std::string::npos);

		sourceDataFile >> citiesNumber;

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
	Point a_geo(convert_to_geo(a.x), convert_to_geo(a.y));
	Point b_geo(convert_to_geo(b.x), convert_to_geo(b.y));

	double q1 = std::cos(a_geo.y - b_geo.y);
	double q2 = std::cos(a_geo.x - b_geo.x);
	double q3 = std::cos(a_geo.x + b_geo.x);

	return (int)RRR * std::acos(0.5 * ((1.0 + q1)*q2 - (1.0 - q1) * q3)) + 1.0;
}
