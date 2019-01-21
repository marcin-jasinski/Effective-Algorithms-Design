#include "pch.h"
#include "Point.h"
#include "Stopwatch.h"
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

std::vector<std::vector<int>> runGeneticAlgorithm(int, int, int, double, int, double);
void runGAbenchmarks();

std::vector<int> queryModelTable(std::vector<int>);
std::vector<int> getNextSpecimen(int);
std::vector<std::vector<int>> getStartingPopulation(int, int);
std::vector<std::vector<int>> selectSpecimenForCrossing(std::vector<std::vector<int>>, int, int);
std::vector<std::vector<int>> crossoverPopulation(std::vector<std::vector<int>>, double, int, int);
std::vector<std::vector<int>> mutatePopulation(std::vector<std::vector<int>>, double, int, int);

void printPopulation(std::vector<std::vector<int>>);
void printSpecimen(std::vector<int>);

int getRouteCost(std::vector<int>, int**, int);
std::pair<std::vector<int>, int> getBestSpecimen(std::vector<std::vector<int>>);
double averageInPopulation(std::vector<std::vector<int>>);

/////////////////////////////////////////////////////////////////////////////////////////////////////
// main

int main()
{
	srand(time(NULL));

	edgesMatrix = readSymetricTSPData("TSP_Data/burma14.tsp");
	correctSolution = 3323;

	// edgesMatrix = readSymetricTSPData("TSP_Data/gr96.tsp");
	// correctSolution = 55209;

	int populationSize = 50;
	int maxPopulations = 1000;
	int crossoverType = 1;
	double crossRatio = 0.6;
	int mutationType = 1;
	double mutationRatio = 0.1;

	std::cout << "Starting genetic algorithm..." << std::endl; 

	std::vector<std::vector<int>> population = runGeneticAlgorithm(populationSize, maxPopulations, crossoverType, crossRatio, mutationType, mutationRatio);

	std::pair<std::vector<int>, int> bestSpecimen = getBestSpecimen(population);
	int bestCost = bestSpecimen.second;
	double error = ((bestSpecimen.second - correctSolution) * 100.0 / correctSolution);

	std::cout << "Done!" << std::endl;
	std::cout << "Best cost: " << bestCost << std::endl;
	std::cout << "Error:     " << error << "%" << std::endl;
	std::cout << "Best specimen: ";
	printSpecimen(bestSpecimen.first);
	return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
// genetic algorithm

std::vector<std::vector<int>> runGeneticAlgorithm(int populationSize, int maxPopulations, int crossoverType, double crossoverRatio, int mutationType, double mutationRatio)
{
	std::vector<std::vector<int>> population;
	std::vector<std::vector<int>> newParentsSet;

	population = getStartingPopulation(populationSize, citiesNumber);

	for (int i = 1; i <= maxPopulations; i++)
	{
		newParentsSet = selectSpecimenForCrossing(population, populationSize, citiesNumber);
		population = crossoverPopulation(newParentsSet, crossoverRatio, crossoverType, populationSize);
		population = mutatePopulation(population, mutationRatio, mutationType, populationSize);
	}

	return population;
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

std::vector<std::vector<int>> crossoverPopulation(std::vector<std::vector<int>> population, double crossoverRatio, int crossoverType, int populationSize)
{
	int parentPopulationSize = populationSize / 2;
	int targetPopulationSize = populationSize;
	int currentPopulationSize = 0;

	std::vector<std::vector<int>> newPopulation;

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


			// PMX
			if (crossoverType == 1)
			{
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

			// OX
			if (crossoverType == 2)
			{
				// initial crossing
				for (int i = crossStartPoint; i < crossEndPoint; i++)
				{
					child1[i] = parent2[i];
					child2[i] = parent1[i];
				}

				// fill remaining values
				for (int i = 0; i < citiesNumber; i++)
				{
					if (child1.at(i) == 0)
					{
						for (int j = 0; j < citiesNumber; j++)
						{
							if (std::find(child1.begin(), child1.end(), parent2.at(j)) == child1.end())
							{
								child1[i] = parent2.at(j);
							}
						}
					}

					if (child2.at(i) == 0)
					{
						for (int j = 0; j < citiesNumber; j++)
						{
							if (std::find(child2.begin(), child2.end(), parent1.at(j)) == child2.end())
							{
								child2[i] = parent1.at(j);
							}
						}
					}
				}
			}
		}

		newPopulation.push_back(child1);
		newPopulation.push_back(child2);
		currentPopulationSize += 2;
	}

	return newPopulation;
}

std::vector<std::vector<int>> mutatePopulation(std::vector<std::vector<int>> population, double mutationRatio, int mutationType, int populationSize)
{
	std::mt19937_64 rng;
	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
	rng.seed(ss);
	std::uniform_real_distribution<double> unif(0, 1);

	// Inversion
	if (mutationType == 1)
	{
		for (int i = 0; i < populationSize; i++)
		{
			for (int j = 0; j < citiesNumber; j++)
			{
				double mutationFactor = unif(rng);

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
	}

	// swap
	if (mutationType == 2)
	{
		for (int i = 0; i < populationSize; i++)
		{
			for (int j = 0; j < citiesNumber; j++)
			{
				double mutationFactor = unif(rng);
			
				if (mutationFactor < mutationRatio)
				{
					int mutationStartPoint = j;
					int mutationEndPoint = rand() % citiesNumber;

					std::swap(population[i][mutationStartPoint], population[i][mutationEndPoint]);

					break;
				}
			}
		}
	}
	
	return population;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
// utils

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

std::pair<std::vector<int>, int> getBestSpecimen(std::vector<std::vector<int>> population)
{
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

	return std::make_pair(population[bestCandidate], bestCost);
}

double averageInPopulation(std::vector<std::vector<int>> population)
{
	double avg = 0.0;
	for (int i = 0; i < population.size(); i++)
	{
		avg += getRouteCost(population.at(i), edgesMatrix, citiesNumber);
	}

	return 1.0 * avg / population.size();
}

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

/////////////////////////////////////////////////////////////////////////////////////////////////////
// rading tsp files

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

/////////////////////////////////////////////////////////////////////////////////////////////////////
// benchmarks

void benchmarkGA(std::fstream &file, int correctSolution, int populationSize, int maxPopulations, int crossoverType, double crossoverRatio, int mutationType, double mutationRatio)
{
	std::vector<std::vector<int>> result;
	Stopwatch timer = Stopwatch();
	
	timer.StartCounter();
	result = runGeneticAlgorithm(populationSize, maxPopulations, crossoverType, crossoverRatio, mutationType, mutationRatio);
	double endTime = timer.GetCounter();

	std::pair<std::vector<int>, int> bestSpecimen = getBestSpecimen(result);

	int bestCost = bestSpecimen.second;
	double average = averageInPopulation(result);
	double error = ((bestSpecimen.second - correctSolution) * 100.0 / correctSolution);

	std::cout << "Best cost: " << bestCost << std::endl;
	std::cout << "Error:     " << error << "%" << std::endl;
	std::cout << "Average :  " << average << std::endl;
	std::cout << "Time:      " << endTime << " seconds" << std::endl;

	file << populationSize << "," << maxPopulations << "," << bestCost << "," << error << "," << average << "," << endTime << std::endl;
}

void runInstanceBenchmarkSet(std::fstream &file, int correctSolution)
{
	system("cls");
	std::cout << "Starting for population 50 " << std::endl;

	benchmarkGA(file, correctSolution, 50, 250, 1, 0.7, 1, 0.1);
	benchmarkGA(file, correctSolution, 50, 500, 1, 0.7, 1, 0.1);
	benchmarkGA(file, correctSolution, 50, 750, 1, 0.7, 1, 0.1);
	benchmarkGA(file, correctSolution, 50, 1000, 1, 0.7, 1, 0.1);

	system("cls");
	std::cout << "Starting for population 100 " << std::endl;

	benchmarkGA(file, correctSolution, 100, 250, 1, 0.7, 1, 0.1);
	benchmarkGA(file, correctSolution, 100, 500, 1, 0.7, 1, 0.1);
	benchmarkGA(file, correctSolution, 100, 750, 1, 0.7, 1, 0.1);
	benchmarkGA(file, correctSolution, 100, 1000, 1, 0.7, 1, 0.1);

	system("cls");
	std::cout << "Starting for population 150 " << std::endl;

	benchmarkGA(file, correctSolution, 150, 250, 1, 0.7, 1, 0.1);
	benchmarkGA(file, correctSolution, 150, 500, 1, 0.7, 1, 0.1);
	benchmarkGA(file, correctSolution, 150, 750, 1, 0.7, 1, 0.1);
	benchmarkGA(file, correctSolution, 150, 1000, 1, 0.7, 1, 0.1);

	system("cls");
	std::cout << "Starting for population 200 " << std::endl;

	benchmarkGA(file, correctSolution, 200, 250, 1, 0.7, 1, 0.1);
	benchmarkGA(file, correctSolution, 200, 500, 1, 0.7, 1, 0.1);
	benchmarkGA(file, correctSolution, 200, 750, 1, 0.7, 1, 0.1);
	benchmarkGA(file, correctSolution, 200, 1000, 1, 0.7, 1, 0.1);

	std::cout << "Finished for population 200 " << std::endl;
}

void runCrossingRatioBenchmarkSet(std::fstream &file, int correctSolution)
{
	system("cls");
	std::cout << "Cross ratio = 0.3" << std::endl;
	file << "Cross_ratio,0.3" << std::endl;
	benchmarkGA(file, correctSolution, 100, 500, 1, 0.3, 1, 0.1);

	system("cls");
	std::cout << "Cross ratio = 0.6" << std::endl;
	file << "Cross_ratio,0.6" << std::endl; 
	benchmarkGA(file, correctSolution, 100, 500, 1, 0.6, 1, 0.1);

	system("cls");
	std::cout << "Cross ratio = 0.9" << std::endl;
	file << "Cross_ratio,0.9" << std::endl;
	benchmarkGA(file, correctSolution, 100, 500, 1, 0.9, 1, 0.1);

	system("cls");
	std::cout << "Cross ratio = 1.0" << std::endl;
	file << "Cross_ratio,1.0" << std::endl;
	benchmarkGA(file, correctSolution, 100, 500, 1, 1.0, 1, 0.1);
}

void runCrossingMutationTypesBenchmarkSet(std::fstream &file, int correctSolution)
{
	system("cls");
	std::cout << "Crossover type: PMX" << std::endl;
	file << "Crossover_type,PMX" << std::endl;
	benchmarkGA(file, correctSolution, 100, 250, 1, 0.6, 1, 0.1);
	benchmarkGA(file, correctSolution, 100, 500, 1, 0.6, 1, 0.1);
	benchmarkGA(file, correctSolution, 100, 750, 1, 0.6, 1, 0.1);
	benchmarkGA(file, correctSolution, 100, 1000, 1, 0.6, 1, 0.1);

	system("cls");
	std::cout << "Crossover type: OX" << std::endl;
	file << "Crossover_type,OX" << std::endl;
	benchmarkGA(file, correctSolution, 100, 250, 2, 0.6, 1, 0.1);
	benchmarkGA(file, correctSolution, 100, 500, 2, 0.6, 1, 0.1);
	benchmarkGA(file, correctSolution, 100, 750, 2, 0.6, 1, 0.1);
	benchmarkGA(file, correctSolution, 100, 1000, 2, 0.6, 1, 0.1);

	/////////////////////////////////////////////////////////////////////////////////

	system("cls");
	std::cout << "Mutation type: invert" << std::endl;
	file << "Mutation_type,invert" << std::endl;
	benchmarkGA(file, correctSolution, 100, 250, 1, 0.6, 1, 0.1);
	benchmarkGA(file, correctSolution, 100, 500, 1, 0.6, 1, 0.1);
	benchmarkGA(file, correctSolution, 100, 750, 1, 0.6, 1, 0.1);
	benchmarkGA(file, correctSolution, 100, 1000, 1, 0.6, 1, 0.1);

	system("cls");
	std::cout << "Mutation type: swap" << std::endl;
	file << "Mutation_type,swap" << std::endl;
	benchmarkGA(file, correctSolution, 100, 250, 1, 0.6, 2, 0.1);
	benchmarkGA(file, correctSolution, 100, 500, 1, 0.6, 2, 0.1);
	benchmarkGA(file, correctSolution, 100, 750, 1, 0.6, 2, 0.1);
	benchmarkGA(file, correctSolution, 100, 1000, 1, 0.6, 2, 0.1);
}

void runGAbenchmarks()
{
	std::fstream geneticAlgorithm_ATSP_Instances;
	geneticAlgorithm_ATSP_Instances.open("geneticAlgorithm_benchmark_ATSP_Instances.txt", std::ios::out | std::ios::trunc);
	if (!geneticAlgorithm_ATSP_Instances.good())
	{
		std::cout << "Error opening file..." << std::endl;
		return;
	}
	else
	{
		geneticAlgorithm_ATSP_Instances << "population_size,max_populations,best_cost,error,average_population,time" << std::endl;

		edgesMatrix = readAsymetricTSPData("ATSP_Data/br17.atsp");
		correctSolution = 39;
		runInstanceBenchmarkSet(geneticAlgorithm_ATSP_Instances, correctSolution);
		
		edgesMatrix = readAsymetricTSPData("ATSP_Data/ftv35.atsp");
		correctSolution = 1473;
		runInstanceBenchmarkSet(geneticAlgorithm_ATSP_Instances, correctSolution);

		edgesMatrix = readAsymetricTSPData("ATSP_Data/ftv64.atsp");
		correctSolution = 1839;
		runInstanceBenchmarkSet(geneticAlgorithm_ATSP_Instances, correctSolution);

		edgesMatrix = readAsymetricTSPData("ATSP_Data/kro124p.atsp");
		correctSolution = 36230;
		runInstanceBenchmarkSet(geneticAlgorithm_ATSP_Instances, correctSolution);

		edgesMatrix = readAsymetricTSPData("ATSP_Data/ftv170.atsp");
		correctSolution = 2755;
		runInstanceBenchmarkSet(geneticAlgorithm_ATSP_Instances, correctSolution);

		geneticAlgorithm_ATSP_Instances.close();
	}

	std::fstream geneticAlgorithm_TSP_Instances;
	geneticAlgorithm_TSP_Instances.open("geneticAlgorithm_benchmark_TSP_Instances.txt", std::ios::out | std::ios::trunc);
	if (!geneticAlgorithm_TSP_Instances.good())
	{
		std::cout << "Error opening file..." << std::endl;
		return;
	}
	else
	{
		geneticAlgorithm_TSP_Instances << "population_size,max_populations,best_cost,error,average_population,time" << std::endl;

		edgesMatrix = readSymetricTSPData("TSP_Data/burma14.tsp");
		correctSolution = 3323;
		runInstanceBenchmarkSet(geneticAlgorithm_TSP_Instances, correctSolution);

		edgesMatrix = readSymetricTSPData("TSP_Data/ulysses22.tsp");
		correctSolution = 7013;
		runInstanceBenchmarkSet(geneticAlgorithm_TSP_Instances, correctSolution);

		edgesMatrix = readSymetricTSPData("TSP_Data/gr96.tsp");
		correctSolution = 55209;
		runInstanceBenchmarkSet(geneticAlgorithm_TSP_Instances, correctSolution);

		edgesMatrix = readSymetricTSPData("TSP_Data/gr137.tsp");
		correctSolution = 69853;
		runInstanceBenchmarkSet(geneticAlgorithm_TSP_Instances, correctSolution);

		edgesMatrix = readSymetricTSPData("TSP_Data/gr202.tsp");
		correctSolution = 40160;
		runInstanceBenchmarkSet(geneticAlgorithm_TSP_Instances, correctSolution);

		geneticAlgorithm_TSP_Instances.close();
	}

	std::fstream geneticAlgorithm_ATSP_CrossRatio;
	geneticAlgorithm_ATSP_CrossRatio.open("geneticAlgorithm_benchmark_ATSP_CrossRatio.txt", std::ios::out | std::ios::trunc);
	if (!geneticAlgorithm_ATSP_CrossRatio.good())
	{
		std::cout << "Error opening file..." << std::endl;
		return;
	}
	else
	{
		geneticAlgorithm_ATSP_CrossRatio << "population_size,max_populations,best_cost,error,average_population,time" << std::endl;

		edgesMatrix = readAsymetricTSPData("ATSP_Data/ftv35.atsp");
		correctSolution = 1473;
		runCrossingRatioBenchmarkSet(geneticAlgorithm_ATSP_CrossRatio, correctSolution);

		edgesMatrix = readAsymetricTSPData("ATSP_Data/ftv170.atsp");
		correctSolution = 2755;
		runCrossingRatioBenchmarkSet(geneticAlgorithm_ATSP_CrossRatio, correctSolution);

		geneticAlgorithm_ATSP_CrossRatio.close();
	}

	std::fstream geneticAlgorithm_TSP_CrossRatio;
	geneticAlgorithm_TSP_CrossRatio.open("geneticAlgorithm_benchmark_TSP_CrossRatio.txt", std::ios::out | std::ios::trunc);
	if (!geneticAlgorithm_TSP_CrossRatio.good())
	{
		std::cout << "Error opening file..." << std::endl;
		return;
	}
	else
	{
		geneticAlgorithm_TSP_CrossRatio << "population_size,max_populations,best_cost,error,average_population,time" << std::endl;

		edgesMatrix = readSymetricTSPData("TSP_Data/ulysses22.tsp");
		correctSolution = 7013;
		runInstanceBenchmarkSet(geneticAlgorithm_TSP_CrossRatio, correctSolution);

		edgesMatrix = readSymetricTSPData("TSP_Data/gr202.tsp");
		correctSolution = 40160;
		runInstanceBenchmarkSet(geneticAlgorithm_TSP_CrossRatio, correctSolution);

		geneticAlgorithm_TSP_CrossRatio.close();
	}

	std::fstream geneticAlgorithm_TSP_Parameters;
	geneticAlgorithm_TSP_Parameters.open("geneticAlgorithm_benchmark_TSP_Parameters.txt", std::ios::out | std::ios::trunc);
	if (!geneticAlgorithm_TSP_Parameters.good())
	{
		std::cout << "Error opening file..." << std::endl;
		return;
	}
	else
	{
		geneticAlgorithm_TSP_Parameters << "population_size,max_populations,best_cost,error,average_population,time" << std::endl;

		edgesMatrix = readSymetricTSPData("TSP_Data/gr96.tsp");
		correctSolution = 55209;
		runCrossingMutationTypesBenchmarkSet(geneticAlgorithm_TSP_Parameters, correctSolution);

		geneticAlgorithm_TSP_Parameters.close();
	}
}
