#include "pch.h"
#include "Point.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <string>
#include <random>
#include <chrono>

void printSpecimen(std::vector<int>);
void printPath(std::vector<int>);
int** readAsymetricTSPData(std::string);
int** readSymetricTSPData(std::string);
int nint(double d);
double convert_to_geo(double);
int TWOD_geo_distance(const Point, const Point);

int** edgesMatrix;
int citiesNumber;
int correctSolution;

const double PI = 3.141592;
const double RRR = 6378.388;

std::vector<int> solution;

std::vector<int> queryModelTable(std::vector<int>);
std::vector<int> getNextSpecimen(int);
std::vector<std::vector<int>> getStartingPopulation(int, int);
std::vector<std::vector<int>> selectSpecimenForCrossing(std::vector<std::vector<int>>, int, int);
std::vector<std::vector<int>> crossoverPopulation(std::vector<std::vector<int>>, double, int);
std::vector<std::vector<int>> mutatePopulation(std::vector<std::vector<int>>, double, int);

void printPopulation(std::vector<std::vector<int>>);
void printSpecimen(std::vector<int>);

int getRouteCost(std::vector<int>, int**, int);

int main()
{
	srand(time(NULL));

	std::cout << "Starting genetic algorithm... " << std::endl;

	int populationSize = 100;
	double crossoverRatio = 0.2;
	double mutationRatio = 0.05;

	edgesMatrix = readSymetricTSPData("TSP_Data/burma14.tsp");
	correctSolution = 3323;

	std::vector<std::vector<int>> population;
	std::vector<std::vector<int>> newParentsSet;

	population = getStartingPopulation(populationSize, citiesNumber);

	for (int i = 1; i <= 1000; i++)
	{
		std::cout << "Population: " << i << ", ";

		newParentsSet = selectSpecimenForCrossing(population, populationSize, citiesNumber);
		population = crossoverPopulation(newParentsSet, crossoverRatio, populationSize);
		population = mutatePopulation(population, mutationRatio, populationSize);

		int bestCost = INT_MAX;
		int bestCandidate = 0;

		for (int i = 0; i < population.size(); i++)
		{
			std::vector<int> candidate = population[i];
			int candidateCost = getRouteCost(candidate, edgesMatrix, citiesNumber);

			if (candidateCost < bestCost)
			{
				bestCost = candidateCost;
				bestCandidate = i;
			}
		}

		std::cout << "best cost: " << bestCost << std::endl;
	}
	
	int bestCost = INT_MAX;
	int bestCandidate = 0;

	for (int i = 0; i < population.size(); i++)
	{
		std::vector<int> candidate = population[i];
		int candidateCost = getRouteCost(candidate, edgesMatrix, citiesNumber);

		if (candidateCost < bestCost)
		{
			bestCost = candidateCost;
			bestCandidate = i;
		}
	}

	std::cout << "Best specimen: " << std::endl;
	printSpecimen(population[bestCandidate]);
	
	std::cout << "Correct solution:" << correctSolution << std::endl;
	std::cout << "Error: " << ((bestCost - correctSolution) * 100.0 / correctSolution) << "%" << std::endl;
	
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
			genome = (std::rand() % problemSize) + 1;
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

std::vector<std::vector<int>> selectSpecimenForCrossing(std::vector<std::vector<int>> population, int populationSize, int problemSize)
{
	std::vector<std::vector<int>> newParentsSet;

	int targetParentsCount = populationSize;
	int currentPopulationSize = 0;
	
	std::vector<int> bufferedSpecimen;
	std::vector<std::vector<int>> populationBuffer;

	while (currentPopulationSize < targetParentsCount)
	{
		while (populationBuffer.size() <= population.size() / 4)
		{
			int specimenCandidateIndex = rand() % populationSize;

			if (std::find(bufferedSpecimen.begin(), bufferedSpecimen.end(), specimenCandidateIndex) == bufferedSpecimen.end())
			{
				bufferedSpecimen.push_back(specimenCandidateIndex);
				populationBuffer.push_back(population.at(specimenCandidateIndex));
			}
		}

		int bestCost = INT_MAX;
		int bestCandidate = 0;

		for (int i = 0; i < bufferedSpecimen.size(); i++)
		{
			std::vector<int> candidate = populationBuffer[i];
			int candidateCost = getRouteCost(candidate, edgesMatrix, citiesNumber);

			if (candidateCost < bestCost)
			{
				bestCost = candidateCost;
				bestCandidate = i;
			}
		}

		newParentsSet.push_back(population[bestCandidate]);
		currentPopulationSize++;

		bufferedSpecimen.clear();
		populationBuffer.clear();
	}

	return newParentsSet;
}

std::vector<std::vector<int>> crossoverPopulation(std::vector<std::vector<int>> population, double crossoverRatio, int populationSize)
{
	int parentPopulationSize = populationSize / 2;
	int targetPopulationSize = populationSize;
	int currentPopulationSize = 0;

	std::vector<std::vector<int>> newPopulation;

	// PMX
	while (currentPopulationSize < targetPopulationSize)
	{
		std::vector<int> parent1, parent2;

		do
		{
			int parent1_index = (rand() % parentPopulationSize);
			int parent2_index = (rand() % parentPopulationSize);

			parent1 = population.at(parent1_index);
			parent2 = population.at(parent2_index);

		} while (parent1 == parent2);
		
		std::vector<int> child1(parent1.begin(), parent1.end());
		std::vector<int> child2(parent2.begin(), parent2.end());

		std::mt19937_64 rng;
		uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
		std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
		rng.seed(ss);
		std::uniform_real_distribution<double> unif(0, 1);

		double crossoverFactor = unif(rng);

		if (crossoverFactor < crossoverRatio)
		{
			
			int crossStartPoint = rand() % citiesNumber + 1;
			int crossEndPoint = rand() % citiesNumber + 1;

			if (crossStartPoint > crossEndPoint) std::swap(crossStartPoint, crossEndPoint);

			std::vector<std::pair<int, int>> modelTable;

			for (int i = crossStartPoint; i < crossEndPoint; i++)
			{
				auto modelPair = std::make_pair(0, 0);
				modelPair.first = parent1.at(i);
				modelPair.second = parent2.at(i);
				modelTable.push_back(modelPair);
			}

			std::fill(child1.begin(), child1.end(), 0);
			std::fill(child2.begin(), child2.end(), 0);

			// initial crossing
			for (int i = crossStartPoint; i < crossEndPoint; i++)
			{
				child1[i] = parent2[i];
				child2[i] = parent1[i];
			}

			// inserting "no conflict" values from original parents
			for (int i = 0; i < citiesNumber; i++)
			{
				if (child1.at(i) == 0)
				{
					int valueInParent1 = parent1.at(i);
					if (std::find(child1.begin(), child1.end(), valueInParent1) == child1.end())
					{
						child1[i] = valueInParent1;
					}
				}

				if (child2.at(i) == 0)
				{
					int valueInParent2 = parent2.at(i);
					if (std::find(child2.begin(), child2.end(), valueInParent2) == child2.end())
					{
						child2[i] = valueInParent2;
					}
				}
			}

			// inserting values from model table
			for (int i = 0; i < citiesNumber; i++)
			{
				if (child1.at(i) == 0)
				{
					int valueInParent1 = parent1.at(i);

					while (true)
					{
						int indexInParent2 = std::distance(parent2.begin(), std::find(parent2.begin(), parent2.end(), valueInParent1));

						int mappedTo = parent1.at(indexInParent2);

						if (std::find(child1.begin(), child1.end(), mappedTo) == child1.end())
						{
							child1[i] = mappedTo;
							break;
						}
						else
						{
							valueInParent1 = mappedTo;
						}
					}
				}

				if (child2.at(i) == 0)
				{
					int valueInParent2 = parent2.at(i);

					while (true)
					{
						int indexInParent1 = std::distance(parent1.begin(), std::find(parent1.begin(), parent1.end(), valueInParent2));

						int mappedTo = parent2.at(indexInParent1);

						if (std::find(child2.begin(), child2.end(), mappedTo) == child2.end())
						{
							child2[i] = mappedTo;
							break;
						}
						else
						{
							valueInParent2 = mappedTo;
						}
					}
				}
			}

			modelTable.clear();
		}

		newPopulation.push_back(child1);
		newPopulation.push_back(child2);
		currentPopulationSize += 2;
	}

	return newPopulation;
}

std::vector<std::vector<int>> mutatePopulation(std::vector<std::vector<int>> population, double mutationRatio, int populationSize)
{
	std::mt19937_64 rng;
	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
	rng.seed(ss);
	std::uniform_real_distribution<double> unif(0, 1);

	// Inversion
	for (int i = 0; i < populationSize; i++)
	{
		double mutationFactor = unif(rng);

		for (int j = 0; j < citiesNumber; j++)
		{
			if (mutationFactor < mutationRatio)
			{
				int mutationStartPoint = j;
				int mutationEndPoint = rand() % citiesNumber + 1;

				if (mutationStartPoint > mutationEndPoint) std::swap(mutationStartPoint, mutationEndPoint);

				for (int k = mutationStartPoint; k < mutationEndPoint / 2; k++)
				{
					std::swap(population[i][k], population[i][mutationEndPoint - k - 1]);
				}

				break;
			}
		}
	}

	return population;
}

int getRouteCost(std::vector<int> specimen, int** edgesMatrix, int citiesNumber)
{
	int totalCost = 0;
	int startingCityIndex = specimen.at(0) - 1;
	int lastCityIndex = specimen.at(specimen.size() - 1) - 1;

	for (int i = 0; i < specimen.size() - 1; i++)
	{
		int currentCity = specimen.at(i) - 1;
		int nextOnRoute = specimen.at(i + 1) - 1;

		totalCost += edgesMatrix[currentCity][nextOnRoute];
	}
	totalCost += edgesMatrix[lastCityIndex][startingCityIndex];

	return totalCost;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

void printSpecimen(std::vector<int> specimen)
{
	for (int i = 0; i < specimen.size(); i++)
	{
		std::cout << std::setw(3) << specimen[i] << " ";
	}
	std::cout << "   cost: " << getRouteCost(specimen, edgesMatrix, citiesNumber) << std::endl;
}

void printPopulation(std::vector<std::vector<int>> population)
{
	for (int i = 0; i < population.size(); i++)
	{
		printSpecimen(population[i]);
	}
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
