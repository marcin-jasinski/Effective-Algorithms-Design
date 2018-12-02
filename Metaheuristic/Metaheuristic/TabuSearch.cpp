#include "pch.h"
#include <algorithm>
#include <ctime>
#include <random>
#include <iostream>
#include <iomanip>
#include "TabuSearch.h"

#define epsilon 0.0001

std::vector<int> getRandomPermutationTabu(int citiesNumber)
{
	std::srand(unsigned(std::time(0)));
	std::vector<int> randomPermutation;
	for (int i = 1; i <= citiesNumber; ++i) randomPermutation.push_back(i);
	std::random_shuffle(randomPermutation.begin(), randomPermutation.end());

	return randomPermutation;
}

std::vector<std::vector<int>> getNeighbourhood(std::vector<int> currentPermutation, int citiesNumber)
{
	int neighboursCount = (citiesNumber * citiesNumber) / 2;
	int i = 1;
	int begin = 0;
	std::vector<std::vector<int>> neighbourhood;
	
	while(begin != citiesNumber - 1)
	{
		for (int currentSwap = begin + 1; currentSwap < citiesNumber; currentSwap++)
		{
			std::vector<int> nextNeighbour = currentPermutation;
			std::swap(nextNeighbour.at(begin), nextNeighbour.at(currentSwap));
			neighbourhood.push_back(nextNeighbour);
		}
		
		begin++;
	}

	return neighbourhood;
}

int fitness(std::vector<int> permutation, int** edgesMatrix, int citiesNumber)
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

std::vector<int> tabuSearch(int citiesNumber, int** edgesMatrix)
{
	std::vector<int> s0 = getRandomPermutationTabu(citiesNumber);
	std::vector<int> sBest = s0;
	std::vector<int> bestCandidate = s0;

	std::vector<std::vector<int>> tabuList;
	tabuList.push_back(s0);

	int iteration = 1;
	while (iteration <= 100)
	{
		std::vector<std::vector<int>> sNeighbourhood = getNeighbourhood(bestCandidate, citiesNumber);
		for (int i = 0; i < sNeighbourhood.size(); i++)
		{
			std::vector<int> sCandidate = sNeighbourhood.at(i);
			if (std::find(tabuList.begin(), tabuList.end(), sCandidate) == tabuList.end() //
				&& fitness(sCandidate, edgesMatrix, citiesNumber) < fitness(bestCandidate, edgesMatrix, citiesNumber)) {

				bestCandidate = sCandidate;
			}
		}

		if(fitness(bestCandidate, edgesMatrix, citiesNumber) < fitness(sBest, edgesMatrix, citiesNumber))
		{
			sBest = bestCandidate;
		}

		tabuList.push_back(bestCandidate);
		if (tabuList.size() > 7) tabuList.erase(tabuList.begin());

		iteration++;
	}

	return sBest;
}