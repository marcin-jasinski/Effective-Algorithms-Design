#include "pch.h"
#include "BranchAndBound.h"
#include "Stopwatch.h"
#include "Node.h"

#include <iostream>
#include <algorithm>
#include <vector>
#include <queue>
#include <utility>
#include <cstdlib>
#include <ctime>
#include <windows.h>

struct nodeComparator
{
	bool operator()(const Node& firstNode, const Node& secondNode) const
	{
		return firstNode.estimatedLowerBound > secondNode.estimatedLowerBound;
	}
};

void branchAndBound(int citiesNumber, int** edgesMatrix)
{
	std::priority_queue<Node, std::vector<Node>, nodeComparator> nodesQueue;
	std::vector<std::pair<int, int>> v;

	Node rootNode = Node(edgesMatrix, 0, -1, 0, citiesNumber);

	rootNode.estimatedLowerBound= calculateLowerBound(edgesMatrix, citiesNumber);

	nodesQueue.push(rootNode);

	double endTime = 0;
	Stopwatch timer = Stopwatch();
	timer.StartCounter();
	timer.StartCounter();
	while (!nodesQueue.empty())
	{
		Node minimalCostNode = nodesQueue.top();
		nodesQueue.pop();

		int currentCityNumber = minimalCostNode.getCityNumber();
		if (minimalCostNode.getTreeLevel() == citiesNumber - 1)
		{
			double endTime = timer.GetCounter();
			std::cout << "\n Elapsed time: " << endTime << " seconds." << std::endl;
			std::cout << " Best route cost: " << minimalCostNode.estimatedLowerBound << std::endl;
			return;
		}

		for (int j = 0; j < citiesNumber; j++)
		{
			if (minimalCostNode.reducedNodeMatrix[currentCityNumber][j] != INT_MAX)
			{
				Node childNode = Node(minimalCostNode.reducedNodeMatrix, minimalCostNode.getTreeLevel() + 1, currentCityNumber, j, citiesNumber);

				int childNodeLowerBound = minimalCostNode.getEstimatedLowerBound() + minimalCostNode.reducedNodeMatrix[currentCityNumber][j] + calculateLowerBound(childNode.reducedNodeMatrix, citiesNumber);
				childNode.estimatedLowerBound = childNodeLowerBound;
				nodesQueue.push(childNode);
			}
		}
	}
}

int* reduceMatrixRow(int** nodeEdgesMatrix, int citiesNumber)
{
	int* reductionRow = new int[citiesNumber];
	for (int i = 0; i < citiesNumber; i++)
	{
		reductionRow[i] = INT_MAX;
	}

	for (int i = 0; i < citiesNumber; i++)
	{
		for (int j = 0; j < citiesNumber; j++)
		{
			if (nodeEdgesMatrix[i][j] < reductionRow[i]) reductionRow[i] = nodeEdgesMatrix[i][j];
		}
	}

	for (int i = 0; i < citiesNumber; i++)
	{
		for (int j = 0; j < citiesNumber; j++)
		{
			if (nodeEdgesMatrix[i][j] != INT_MAX && reductionRow[i] != INT_MAX) nodeEdgesMatrix[i][j] -= reductionRow[i];
		}
	}

	return reductionRow;
}

int* reduceMatrixCol(int** nodeEdgesMatrix, int citiesNumber)
{
	int* reductionRow = new int[citiesNumber];
	for (int i = 0; i < citiesNumber; i++)
	{
		reductionRow[i] = INT_MAX;
	}

	for (int i = 0; i < citiesNumber; i++)
	{
		for (int j = 0; j < citiesNumber; j++)
		{
			if (nodeEdgesMatrix[i][j] < reductionRow[j]) reductionRow[j] = nodeEdgesMatrix[i][j];
		}
	}

	for (int i = 0; i < citiesNumber; i++)
	{
		for (int j = 0; j < citiesNumber; j++)
		{
			if (nodeEdgesMatrix[i][j] != INT_MAX && reductionRow[j] != INT_MAX) nodeEdgesMatrix[i][j] -= reductionRow[j];
		}
	}

	return reductionRow;
}

int calculateLowerBound(int** nodeEdgesMatrix, int citiesNumber)
{
	int lowerBound = 0;
	int* rowReductor;
	int* colReductor;

	rowReductor = reduceMatrixRow(nodeEdgesMatrix, citiesNumber);
	colReductor = reduceMatrixCol(nodeEdgesMatrix, citiesNumber);

	for (int i = 0; i < citiesNumber; i++)
	{
		if (rowReductor[i] != INT_MAX) lowerBound += rowReductor[i];
		if (colReductor[i] != INT_MAX) lowerBound += colReductor[i];
	}

	return lowerBound;
}