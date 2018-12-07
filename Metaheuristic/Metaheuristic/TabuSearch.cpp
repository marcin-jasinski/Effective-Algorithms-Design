#include "pch.h"
#include <algorithm>
#include <ctime>
#include <random>
#include <iostream>
#include <iomanip>
#include "TabuSearch.h"

#define epsilon 0.0001
std::vector<std::vector<int>> moves;

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
	int i = 1;
	int begin = 0;
	std::vector<std::vector<int>> neighbourhood;
	
	while(begin != citiesNumber - 1)
	{
		for (int currentSwap = begin + 1; currentSwap < citiesNumber; currentSwap++)
		{
			std::vector<int> move;
			move.push_back(begin);
			move.push_back(currentSwap);

			std::vector<int> nextNeighbour = currentPermutation;
			std::swap(nextNeighbour.at(begin), nextNeighbour.at(currentSwap));

			neighbourhood.push_back(nextNeighbour);
			moves.push_back(move);
		}
		
		begin++;
	}

	return neighbourhood;
}

int pathCost(std::vector<int> permutation, int** edgesMatrix, int citiesNumber)
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

std::vector<int> tabuSearch(int maxIterations, int tabuMultiplier, int citiesNumber, int** edgesMatrix)
{
	std::vector<int> s0 = getRandomPermutationTabu(citiesNumber);
	std::vector<int> sBest = s0;
	std::vector<int> bestCandidate = s0;

	std::vector<int> sCandidate;
	std::vector<int> moveToCandidate;

	std::vector<std::vector<int>> sNeighbourhood;
	std::vector<std::vector<int>> tabuList;
	std::vector<int> moveToBePutOnTabuList; // czyNieZarezerwowanoCa³egoPrzedzia³uDlaPsa()

	int iterationThreshold = 25;
	int iteration = 1;
	while (iteration <= maxIterations)
	{
		bool noAcceptedCandidates = true;
		sNeighbourhood = getNeighbourhood(bestCandidate, citiesNumber);
		for (int i = 0; i < sNeighbourhood.size(); i++)
		{
			sCandidate = sNeighbourhood.at(i);
			moveToCandidate = moves.at(i);

			if (std::find(tabuList.begin(), tabuList.end(), moveToCandidate) != tabuList.end())
			{
				if (pathCost(sCandidate, edgesMatrix, citiesNumber) < pathCost(sBest, edgesMatrix, citiesNumber))
				{
					bestCandidate = sCandidate;
					moveToBePutOnTabuList = moveToCandidate;
					noAcceptedCandidates = false;
				}
			}
			else
			{
				if (pathCost(sCandidate, edgesMatrix, citiesNumber) < pathCost(bestCandidate, edgesMatrix, citiesNumber))
				{
					noAcceptedCandidates = false;
					bestCandidate = sCandidate;
					moveToBePutOnTabuList = moveToCandidate;
				}
			}
		}

		if (noAcceptedCandidates) {
			iterationThreshold--;
		}

		if (iterationThreshold <= 0)
		{
			sNeighbourhood.clear();
			iterationThreshold = 25;
			bestCandidate = getRandomPermutationTabu(citiesNumber);
		}
		else
		{
			if (pathCost(bestCandidate, edgesMatrix, citiesNumber) < pathCost(sBest, edgesMatrix, citiesNumber))
			{
				sBest = bestCandidate;
			}

			tabuList.push_back(moveToBePutOnTabuList);
			if (tabuList.size() > tabuMultiplier * citiesNumber) tabuList.erase(tabuList.begin());

			iteration++;
		}

		moves.clear();
	}

	return sBest;
}