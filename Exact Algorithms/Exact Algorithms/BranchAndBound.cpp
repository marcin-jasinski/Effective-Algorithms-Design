#include "pch.h"
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

int* reduceMatrixRow(int** reducedMatrix, int citiesNumber)
{
	int* row = new int[citiesNumber];
	for (int i = 0; i < citiesNumber; i++)
	{
		row[i] = INT_MAX;
	}

	// row[i] contains minimum in row i
	for (int i = 0; i < citiesNumber; i++)
		for (int j = 0; j < citiesNumber; j++)
			if (reducedMatrix[i][j] < row[i])
				row[i] = reducedMatrix[i][j];

	// reduce the minimum value from each element in each row
	for (int i = 0; i < citiesNumber; i++)
		for (int j = 0; j < citiesNumber; j++)
			if (reducedMatrix[i][j] != INT_MAX && row[i] != INT_MAX)
				reducedMatrix[i][j] -= row[i];

	return row;
}

int* reduceMatrixCol(int** reducedMatrix, int citiesNumber)
{
	int* col = new int[citiesNumber];
	for (int i = 0; i < citiesNumber; i++)
	{
		col[i] = INT_MAX;
	}

	// col[j] contains minimum in col j
	for (int i = 0; i < citiesNumber; i++)
		for (int j = 0; j < citiesNumber; j++)
			if (reducedMatrix[i][j] < col[j])
				col[j] = reducedMatrix[i][j];

	// reduce the minimum value from each element in each column
	for (int i = 0; i < citiesNumber; i++)
		for (int j = 0; j < citiesNumber; j++)
			if (reducedMatrix[i][j] != INT_MAX && col[j] != INT_MAX)
				reducedMatrix[i][j] -= col[j];

	return col;
}

// Function to get the lower bound on
// on the path starting at current min node
int calculateCost(int** reducedMatrix, int citiesNumber)
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

struct comp {
	bool operator()(const Node* nodeA, const Node* nodeB) const
	{
		return nodeA->lowerBound > nodeB->lowerBound;
	}
};

// Function to solve Traveling Salesman Problem using Branch and Bound
int branchAndBound(int citiesNumber, int** costMatrix)
{
	std::priority_queue<Node*, std::vector<Node*>, comp> pq;
	std::vector<int> v;

	Node* root = new Node(costMatrix, v, 0, -1, 0, citiesNumber);

	root->lowerBound = calculateCost(root->nodeMatrix, citiesNumber);

	// Add root to list of live nodes;
	pq.push(root);

	// Finds a live node with least cost, add its children to list of
	// live nodes and finally deletes it from the list
	while (!pq.empty())
	{
		// Find a live node with least estimated cost
		Node* min = pq.top();

		// The found node is deleted from the list of live nodes
		pq.pop();

		// i stores current city number
		int i = min->currentCityNumber;

		// if all cities are visited
		if (min->graphLevel == citiesNumber - 1)
		{
			// return to starting city
			min->pathToVertex.push_back(0);

			// return optimal cost
			std::cout << "=== Branch & Bound ===" << std::endl;
			std::cout << " Best route cost: " << min->lowerBound << std::endl;
			std::cout << " Optimal path: ";

			for (int i = min->pathToVertex.size() - 1; i >= 0; i--)
			{
				std::cout << min->pathToVertex[i] << " -> ";
			}
			std::cout << min->pathToVertex[min->pathToVertex.size() - 1] << std::endl;
			
			pq.empty();
			return min->lowerBound;
		}

		for (int j = 0; j < citiesNumber; j++)
		{
			if (min->nodeMatrix[i][j] != INT_MAX)
			{
				Node* child = new Node(min->nodeMatrix, min->pathToVertex, min->graphLevel+ 1, i, j, citiesNumber);
				child->lowerBound = min->lowerBound + min->nodeMatrix[i][j] + calculateCost(child->nodeMatrix, citiesNumber);
				pq.push(child);
			}
		}

		delete min;
	}

	return -1;
}
