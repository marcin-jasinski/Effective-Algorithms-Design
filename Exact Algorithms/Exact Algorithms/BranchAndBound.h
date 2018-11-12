#pragma once
#include "Node.h"

// This structure is used as a comparator for nodes priority queue.
struct comp {
	bool operator()(const Node* nodeA, const Node* nodeB) const
	{
		return nodeA->lowerBound > nodeB->lowerBound;
	}
};

int branchAndBound(int citiesNumber, int** costMatrix);
int* reduceMatrixRow(int** reducedMatrix, int citiesNumber);
int* reduceMatrixCol(int** reducedMatrix, int citiesNumber);
int calculateLowerBound(int** reducedMatrix, int citiesNumber);