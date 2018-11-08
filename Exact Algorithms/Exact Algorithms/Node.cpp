#include "pch.h"
#include "Node.h"
#include <climits>

Node::Node()
{
}

Node::~Node()
{
}

Node::Node(int** parentNodeMatrix, int treeLevel, int i, int j, int citiesNumber)
{
	this->reducedNodeMatrix = parentNodeMatrix;

	for (int x = 0; x < citiesNumber && this->treeLevel != 0; x++)
	{
		this->reducedNodeMatrix[i][x] = INT_MAX;
		this->reducedNodeMatrix[x][j] = INT_MAX;
	}
	this->reducedNodeMatrix[j][0] = INT_MAX;

	this->treeLevel = treeLevel;
	this->cityNumber = j;
}

int Node::getTreeLevel()
{
	return this->treeLevel;
}

int Node::getCityNumber()
{
	return this->cityNumber;
}

int Node::getEstimatedLowerBound()
{
	return this->estimatedLowerBound;
}

int** Node::getReducedNodeMatrix()
{
	return this->reducedNodeMatrix;
}
