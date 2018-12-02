#include "pch.h"
#include "SimulatedAnnealing.h"
#include <algorithm>
#include <ctime>
#include <random>
#include <iostream>

#define epsilon 0.0001

std::vector<int> getRandomPermutation(int citiesNumber)
{
	std::srand(unsigned(std::time(0)));
	std::vector<int> randomPermutation;
	for (int i = 1; i <= citiesNumber; ++i) randomPermutation.push_back(i);
	std::random_shuffle(randomPermutation.begin(), randomPermutation.end());

	return randomPermutation;
}

std::vector<int> getNextNeighbour(std::vector<int> currentPermutation, int citiesNumber)
{
	int a, b;
	std::vector<int> nextNeighbour = currentPermutation;
	std::random_device rd; // obtain a random number from hardware
	std::mt19937 eng(rd()); // seed the generator
	std::uniform_int_distribution<> distr(0, citiesNumber - 1); // define the range

	do
	{
		a = distr(eng);
		b = distr(eng);

	} while (a == b);
	std::swap(nextNeighbour.at(a), nextNeighbour.at(b));

	return nextNeighbour;
}

int getRouteCost(std::vector<int> permutation, int** edgesMatrix, int citiesNumber)
{
	int totalCost = 0;
	int startingCityIndex = permutation.at(0) - 1;
	int lastCityIndex = permutation.at(permutation.size() - 1) - 1;

	for (int i = 0; i < permutation.size() - 1; i++)
	{
		int currentCity = permutation.at(i) - 1;
		int nextOnRoute = permutation.at(i + 1) - 1;

		totalCost += edgesMatrix[currentCity][nextOnRoute];
	}
	totalCost += edgesMatrix[lastCityIndex][startingCityIndex];

	return totalCost;
}

std::vector<int> simulatedAnnealing(double maxTemperature, double tempChange, int citiesNumber, int** edgesMatrix)
{
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dist(0.0, 1.0);

	std::vector<int> currentSolution = getRandomPermutation(citiesNumber);
	int currentCost = getRouteCost(currentSolution, edgesMatrix, citiesNumber);

	std::vector<int> bestSolution = currentSolution;
	int bestCost;

	int iteration = 1;
	double temp = maxTemperature;
	std::vector<int> nextNeighbour;
	while (temp > epsilon)
	{	
		currentCost = getRouteCost(currentSolution, edgesMatrix, citiesNumber);
		bestCost = getRouteCost(bestSolution, edgesMatrix, citiesNumber);

		nextNeighbour = getNextNeighbour(currentSolution, citiesNumber);
		int neighbourCost = getRouteCost(nextNeighbour, edgesMatrix, citiesNumber);

		if (neighbourCost <= currentCost)
		{
			currentSolution = nextNeighbour;
			if (neighbourCost < bestCost)
			{
				bestSolution = nextNeighbour;
			}
		}
		else
		{
			double delta = currentCost - neighbourCost;
			double u = dist(mt);

			if (u <= std::exp(1 / temp * delta))
			{
				currentSolution = nextNeighbour;
			}
		}

		temp *= tempChange;
		iteration++;
	}

	return bestSolution;
}
