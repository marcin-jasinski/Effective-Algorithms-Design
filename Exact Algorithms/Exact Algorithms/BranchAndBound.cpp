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
	std::priority_queue<Node*, std::vector<Node*>, comp> nodesPriorityQueue;
	std::vector<int> pathToNode;

	Node* startingNode = new Node(costMatrix, pathToNode, 0, -1, 0, citiesNumber);

	startingNode->lowerBound = calculateLowerBound(startingNode->nodeMatrix, citiesNumber);

	nodesPriorityQueue.push(startingNode);

	while (!nodesPriorityQueue.empty())
	{
		Node* min = nodesPriorityQueue.top();
		nodesPriorityQueue.pop();

		int i = min->currentCityNumber;

		if (min->graphLevel == citiesNumber - 1)
		{
			min->pathToVertex.push_back(0);

			std::cout << "=== Branch & Bound ===" << std::endl;
			std::cout << " Best route cost: " << min->lowerBound << std::endl;
			std::cout << " Optimal path: ";

			for (int i = min->pathToVertex.size() - 1; i >= 0; i--)
			{
				std::cout << min->pathToVertex[i] << " -> ";
			}
			std::cout << min->pathToVertex[min->pathToVertex.size() - 1] << std::endl;

			nodesPriorityQueue.empty();
			return min->lowerBound;
		}

		for (int j = 0; j < citiesNumber; j++)
		{
			if (min->nodeMatrix[i][j] != INT_MAX)
			{
				Node* child = new Node(min->nodeMatrix, min->pathToVertex, min->graphLevel + 1, i, j, citiesNumber);
				child->lowerBound = min->lowerBound + min->nodeMatrix[i][j] + calculateLowerBound(child->nodeMatrix, citiesNumber);
				nodesPriorityQueue.push(child);
			}
		}

		delete min;
	}

	return -1;
}

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


