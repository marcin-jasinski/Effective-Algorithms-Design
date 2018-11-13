#include "pch.h"
#include "BranchAndBound.h"
#include "Stopwatch.h"
#include "Node.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <windows.h>
#include <vector>
#include <queue>
#include <utility>
#include <cstring>
#include <climits>

// This method starts the main B&B algorithm execution.
int branchAndBound(int citiesNumber, int** costMatrix)
{
	double endTime;
	Stopwatch timer = Stopwatch();
	timer.StartCounter();

	// Creating priority queue for nodes
	std::priority_queue<Node*, std::vector<Node*>, comp> nodesPriorityQueue;

	std::vector<int> pathToNode;
	Node* startingNode = new Node(costMatrix, pathToNode, 0, 0, 0, citiesNumber);

	// Calculating lower bound and pushing node to priority queue
	startingNode->lowerBound = calculateLowerBound(startingNode->nodeMatrix, citiesNumber);
	nodesPriorityQueue.push(startingNode);

	while (!nodesPriorityQueue.empty())
	{
		Node* minimalCostNode = nodesPriorityQueue.top();
		nodesPriorityQueue.pop();

		int i = minimalCostNode->currentCityNumber;

		// If the end of the graph has been reached
		if (minimalCostNode->graphLevel == citiesNumber - 1)
		{
			endTime = timer.GetCounter(); // stop timer

			minimalCostNode->pathToVertex.push_back(0);

			std::cout << "\n=== Branch & Bound ===" << std::endl;
			std::cout << " Best route cost: " << minimalCostNode->lowerBound << std::endl;
			std::cout << " Optimal path: ";

			std::cout << "[";
			for (int i = minimalCostNode->pathToVertex.size() - 1; i >= 0; i--)
			{
				std::cout << minimalCostNode->pathToVertex[i] << ",";
			}
			std::cout << minimalCostNode->pathToVertex[minimalCostNode->pathToVertex.size() - 1] << "]" << std::endl;
			std::cout << " \nElapsed time: " << endTime << " seconds." << std::endl;

			return minimalCostNode->lowerBound;
		}

		// Creating child nodes for the current node.
		for (int j = 0; j < citiesNumber; j++)
		{
			if (minimalCostNode->nodeMatrix[i][j] != INT_MAX)
			{
				Node* child = new Node(minimalCostNode->nodeMatrix, minimalCostNode->pathToVertex, minimalCostNode->graphLevel + 1, i, j, citiesNumber);

				// Calculating lower bound for child node
				child->lowerBound = minimalCostNode->lowerBound + minimalCostNode->nodeMatrix[i][j] + calculateLowerBound(child->nodeMatrix, citiesNumber);
				nodesPriorityQueue.push(child);
			}
		}

		delete minimalCostNode;
	}

	return -1;
}

// This method performs matrix reduction using smallest values in each row
int* reduceMatrixRow(int** reducedMatrix, int citiesNumber)
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
			if (reducedMatrix[i][j] < reductionRow[i]) reductionRow[i] = reducedMatrix[i][j];
		}
	}

	for (int i = 0; i < citiesNumber; i++)
	{
		for (int j = 0; j < citiesNumber; j++)
		{
			if (reducedMatrix[i][j] != INT_MAX && reductionRow[i] != INT_MAX) reducedMatrix[i][j] -= reductionRow[i];
		}
	}

	return reductionRow;
}

// This method performs matrix reduction using smallest values in each column
int* reduceMatrixCol(int** reducedMatrix, int citiesNumber)
{
	int* reductionCol = new int[citiesNumber];
	for (int i = 0; i < citiesNumber; i++)
	{
		reductionCol[i] = INT_MAX;
	}

	for (int i = 0; i < citiesNumber; i++)
	{
		for (int j = 0; j < citiesNumber; j++)
		{
			if (reducedMatrix[i][j] < reductionCol[j]) reductionCol[j] = reducedMatrix[i][j];
		}
	}
		
	for (int i = 0; i < citiesNumber; i++)
	{
		for (int j = 0; j < citiesNumber; j++)
		{
			if (reducedMatrix[i][j] != INT_MAX && reductionCol[j] != INT_MAX) reducedMatrix[i][j] -= reductionCol[j];
		}
	}

	return reductionCol;
}

// This method calculates lower bound for a given node
// using values from the reduction methods.
int calculateLowerBound(int** reducedMatrix, int citiesNumber)
{
	int cost = 0;
	int* reductionRow = reduceMatrixRow(reducedMatrix,citiesNumber);
	int* reductionCol = reduceMatrixCol(reducedMatrix, citiesNumber);

	for (int i = 0; i < citiesNumber; i++)
	{
		cost += (reductionRow[i] != INT_MAX) ? reductionRow[i] : 0,
		cost += (reductionCol[i] != INT_MAX) ? reductionCol[i] : 0;
	}

	delete[] reductionRow;
	delete[] reductionCol;

	return cost;
}


