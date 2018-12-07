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
#include "Stopwatch.h"

int correctSolution;
int** edgesMatrix;

double coolingFactor[5] = { 0.99, 0.999, 0.9999, 0.99999, 0.999999 };

int maxIterationsTab[] = { 100, 300, 500, 1000, 2000, 5000, 7500, 10000, 25000, 50000 };
int tabuMultiplierTab[] = { 3, 5, 7, 9, 11 };

std::vector<int> solution;

Stopwatch timer = Stopwatch();

const double PI = 3.141592;
const double RRR = 6378.388;

int citiesNumber;

static int nint(double);
static double convert_to_geo(double);
int TWOD_geo_distance(const Point a, const Point b);

int** readAsymetricTSPData(std::string);
int** readSymetricTSPData(std::string);

void runSimulatedAnnealingBenchmarks();

int main()
{
	std::srand(time(NULL));

	runSimulatedAnnealingBenchmarks();
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

	int **edgesMatrix = nullptr;

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
	Point a_geo( convert_to_geo(a.x), convert_to_geo(a.y));
	Point b_geo( convert_to_geo(b.x), convert_to_geo(b.y));

	double q1 = std::cos(a_geo.y - b_geo.y);
	double q2 = std::cos(a_geo.x - b_geo.x);
	double q3 = std::cos(a_geo.x + b_geo.x);

	return (int)RRR * std::acos(0.5 * ((1.0 + q1)*q2 - (1.0 - q1) * q3)) + 1.0;
}

void benchmarkTS(std::fstream &file)
{
	file << "Number of cities: " << citiesNumber << std::endl;
	for (int i = 0; i < 5; i++)
	{
		double tabuMulti = tabuMultiplierTab[i];
		for (int j = 0; j < 10; j++)
		{
			int maxIterations = maxIterationsTab[i];
			double averageTime = 0;
			system("cls");
			std::cout << "\nNumber of cities: " << citiesNumber << std::endl;
			file << "\nMax iterations," << maxIterations << ",tabu multiplier," << tabuMulti << std::endl;
			std::cout << "Max iterations: " << maxIterations << ", tabu multiplier: " << tabuMulti << std::endl;
			std::cout << std::endl;
			file << "Elapsed time,solution,error rate" << std::endl;
			for (int i = 1; i <= 10; i++)
			{
				timer.StartCounter();
				solution = tabuSearch(maxIterations, tabuMulti, citiesNumber, edgesMatrix);
				double endTime = timer.GetCounter();
				averageTime += endTime;
				double pathCost = getRouteCost(solution, edgesMatrix, citiesNumber);
				double errorRate = ((pathCost - correctSolution) / correctSolution) * 100;
				std::cout << "Elapsed time: " << endTime << " seconds." << std::endl;
				std::cout << "Solution:     " << pathCost << std::endl;
				std::cout << "Error:        " << errorRate << "%" << std::endl;

				file << endTime << "," << pathCost << "," << errorRate << "%" << std::endl;
			}
			
			file << "Average time," << averageTime / 10 << std::endl;
		}
	}
}

void benchmarkSA(std::fstream &file)
{
	file << "Number of cities: " << citiesNumber << std::endl;
	for (int i = 0; i < 5; i++)
	{
		double clf = coolingFactor[i];
		for (int startingTemperature = 250; startingTemperature >= 25; startingTemperature -= 25)
		{
			double averageTime = 0;
			system("cls");
			std::cout << "\nNumber of cities: " << citiesNumber << std::endl;
			file << "\nTemperature," << startingTemperature << ",cooling factor," << clf << std::endl;
			std::cout << "Temperature: " << startingTemperature << ", cooling factor: " << clf << std::endl;
			std::cout << std::endl;
			file << "Elapsed time,solution,error rate" << std::endl;
			for (int i = 1; i <= 10; i++)
			{
				timer.StartCounter();
				solution = simulatedAnnealing(startingTemperature, clf, citiesNumber, edgesMatrix);
				double endTime = timer.GetCounter();
				averageTime += endTime;
				double pathCost = getRouteCost(solution, edgesMatrix, citiesNumber);
				double errorRate = ((pathCost - correctSolution) / correctSolution) * 100;
				std::cout << "Elapsed time: " << endTime << " seconds." << std::endl;
				std::cout << "Solution:     " << pathCost << std::endl;
				std::cout << "Error:        " << errorRate << "%" << std::endl;

				file << endTime << "," << pathCost << "," << errorRate << "%" << std::endl;
			}

			file << "Average time," << averageTime / 10 << std::endl;
		}
	}
}

void runSimulatedAnnealingBenchmarks()
{
	std::fstream annealingSA_ATSP;
	annealingSA_ATSP.open("annealingBenchmark_ATSP.txt", std::ios::out | std::ios::trunc);
	if (!annealingSA_ATSP.good())
	{
		std::cout << "Error opening file..." << std::endl;
		return;
	}
	else
	{
		annealingSA_ATSP << "=== ATSP files benchmark (SA) ===" << std::endl;

		edgesMatrix = readAsymetricTSPData("ATSP_Data/br17.atsp");
		correctSolution = 39;
		benchmarkSA(annealingSA_ATSP);

		edgesMatrix = readAsymetricTSPData("ATSP_Data/ftv35.atsp");
		correctSolution = 1473;
		benchmarkSA(annealingSA_ATSP);

		edgesMatrix = readAsymetricTSPData("ATSP_Data/ftv64.atsp");
		correctSolution = 1839;
		benchmarkSA(annealingSA_ATSP);

		edgesMatrix = readAsymetricTSPData("ATSP_Data/kro124p.atsp");
		correctSolution = 36230;
		benchmarkSA(annealingSA_ATSP);

		edgesMatrix = readAsymetricTSPData("ATSP_Data/ftv170.atsp");
		correctSolution = 2755;
		benchmarkSA(annealingSA_ATSP);

		annealingSA_ATSP.close();
	}
	
	std::fstream annealingSA_TSP;
	annealingSA_TSP.open("annealingBenchmark_TSP.txt", std::ios::out | std::ios::trunc);
	if (!annealingSA_TSP.good())
	{
		std::cout << "Error opening file..." << std::endl;
		return;
	}
	else
	{
		annealingSA_TSP << "=== TSP files benchmark (SA) ===" << std::endl;

		edgesMatrix = readSymetricTSPData("TSP_Data/burma14.tsp");
		correctSolution = 3323;
		benchmarkSA(annealingSA_TSP);

		edgesMatrix = readSymetricTSPData("TSP_Data/ulysses22.tsp");
		correctSolution = 7013;
		benchmarkSA(annealingSA_TSP);

		edgesMatrix = readSymetricTSPData("TSP_Data/gr96.tsp");
		correctSolution = 55209;
		benchmarkSA(annealingSA_TSP);

		edgesMatrix = readSymetricTSPData("TSP_Data/gr137.tsp");
		correctSolution = 69853;
		benchmarkSA(annealingSA_TSP);

		edgesMatrix = readSymetricTSPData("TSP_Data/gr202.tsp");
		correctSolution = 40160;
		benchmarkSA(annealingSA_TSP);

		annealingSA_TSP.close();
	}
}

void runTabuSearchBenchmarks()
{
	std::fstream tabuSearch_ATSP;
	tabuSearch_ATSP.open("annealingBenchmark_ATSP.txt", std::ios::out | std::ios::trunc);
	if (!tabuSearch_ATSP.good())
	{
		std::cout << "Error opening file..." << std::endl;
		return;
	}
	else
	{
		tabuSearch_ATSP << "=== ATSP files benchmark (SA) ===" << std::endl;

		edgesMatrix = readAsymetricTSPData("ATSP_Data/br17.atsp");
		correctSolution = 39;
		benchmarkTS(tabuSearch_ATSP);

		edgesMatrix = readAsymetricTSPData("ATSP_Data/ftv35.atsp");
		correctSolution = 1473;
		benchmarkTS(tabuSearch_ATSP);

		edgesMatrix = readAsymetricTSPData("ATSP_Data/ftv64.atsp");
		correctSolution = 1839;
		benchmarkTS(tabuSearch_ATSP);

		edgesMatrix = readAsymetricTSPData("ATSP_Data/kro124p.atsp");
		correctSolution = 36230;
		benchmarkTS(tabuSearch_ATSP);

		edgesMatrix = readAsymetricTSPData("ATSP_Data/ftv170.atsp");
		correctSolution = 2755;
		benchmarkTS(tabuSearch_ATSP);

		tabuSearch_ATSP.close();
	}

	std::fstream tabuSearch_TSP;
	tabuSearch_TSP.open("annealingBenchmark_TSP.txt", std::ios::out | std::ios::trunc);
	if (!tabuSearch_TSP.good())
	{
		std::cout << "Error opening file..." << std::endl;
		return;
	}
	else
	{
		tabuSearch_TSP << "=== TSP files benchmark (SA) ===" << std::endl;

		edgesMatrix = readSymetricTSPData("TSP_Data/burma14.tsp");
		correctSolution = 3323;
		benchmarkTS(tabuSearch_TSP);

		edgesMatrix = readSymetricTSPData("TSP_Data/ulysses22.tsp");
		correctSolution = 7013;
		benchmarkTS(tabuSearch_TSP);

		edgesMatrix = readSymetricTSPData("TSP_Data/gr96.tsp");
		correctSolution = 55209;
		benchmarkTS(tabuSearch_TSP);

		edgesMatrix = readSymetricTSPData("TSP_Data/gr137.tsp");
		correctSolution = 69853;
		benchmarkTS(tabuSearch_TSP);

		edgesMatrix = readSymetricTSPData("TSP_Data/gr202.tsp");
		correctSolution = 40160;
		benchmarkTS(tabuSearch_TSP);

		tabuSearch_TSP.close();
	}
}
