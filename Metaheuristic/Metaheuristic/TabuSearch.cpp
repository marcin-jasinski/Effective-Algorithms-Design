#include "pch.h"
#include <algorithm>
#include <ctime>
#include <random>
#include <iostream>
#include <iomanip>
#include "TabuSearch.h"

#define epsilon 0.0001
std::vector<std::vector<int>> moves;

std::vector<int> getInitialSolution(int citiesNumber, int** edgesMatrix) 
{
	int startNode = 0;
	int foundNeighbours = 0;

	std::vector<int> initialSolution;
	initialSolution.push_back(startNode);

	while (foundNeighbours < citiesNumber - 1)
	{
		int bestNeighbourIndex = -1;
		int bestNeighbourCost = INT_MAX;

		for (int i = 0; i < citiesNumber; i++)
		{
			if (i != startNode && edgesMatrix[startNode][i] < bestNeighbourCost &&
				std::find(initialSolution.begin(), initialSolution.end(), i) == initialSolution.end())
			{
				bestNeighbourIndex = i;
				bestNeighbourCost = edgesMatrix[startNode][i];
			}
		}

		initialSolution.push_back(bestNeighbourIndex);
		startNode = bestNeighbourIndex;
		foundNeighbours++;
	}

	for (int i = 0; i < citiesNumber; i++) initialSolution.at(i) += 1;

	return initialSolution;
}

std::vector<int> getRandomPermutationTabu(int citiesNumber)
{
	std::srand(unsigned(std::time(0)));
	std::vector<int> randomPermutation;
	for (int i = 0; i < citiesNumber; ++i) randomPermutation.push_back(i);
	std::random_shuffle(randomPermutation.begin(), randomPermutation.end());

	for (int i = 0; i < citiesNumber; i++) randomPermutation.at(i) += 1;
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

int pathCostTabu(std::vector<int> permutation, int** edgesMatrix, int citiesNumber)
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

std::vector<int> tabuSearch(int maxIterations, int citiesNumber, int** edgesMatrix)
{
	int tabuListSize = citiesNumber * 3;
	int iterationThreshold = tabuListSize * 3;

	std::vector<int> s0 = getInitialSolution(citiesNumber, edgesMatrix);
	std::vector<int> sBest = s0;
	std::vector<int> bestCandidate = s0;

	std::vector<int> sCandidate;
	std::vector<int> moveToCandidate;

	std::vector<std::vector<int>> sNeighbourhood;
	std::vector<std::vector<int>> tabuList;
	std::vector<int> moveToBePutOnTabuList; 

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
				if (pathCostTabu(sCandidate, edgesMatrix, citiesNumber) < pathCostTabu(sBest, edgesMatrix, citiesNumber))
				{
					bestCandidate = sCandidate;
					moveToBePutOnTabuList = moveToCandidate;
					noAcceptedCandidates = false;
				}
			}
			else
			{
				if (pathCostTabu(sCandidate, edgesMatrix, citiesNumber) < pathCostTabu(bestCandidate, edgesMatrix, citiesNumber))
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
			moves.clear();
			tabuList.clear();
			sNeighbourhood.clear();
			iterationThreshold = tabuListSize * 3;

			bestCandidate = getRandomPermutationTabu(citiesNumber);
		}
		else
		{
			if (pathCostTabu(bestCandidate, edgesMatrix, citiesNumber) < pathCostTabu(sBest, edgesMatrix, citiesNumber))
			{
				sBest = bestCandidate;
			}

			tabuList.push_back(moveToBePutOnTabuList);
			if (tabuList.size() > tabuListSize) tabuList.erase(tabuList.begin());

			iteration++;
		}

		moves.clear();
	}

	return sBest;
}