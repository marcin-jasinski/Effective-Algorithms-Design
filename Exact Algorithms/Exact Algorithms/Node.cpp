#include "pch.h"
#include "Node.h"
#include <climits>
#include <iostream>
#include <cstdint>
#include <cstring>

Node::Node()
{

}

Node::Node(int** parentMatrix, std::vector<int> const & path, int level, int i, int j, int citiesNumber)
{
	this->pathToVertex = path;
	if (level != 0) this->pathToVertex.push_back(j);

	this->nodeMatrix = new int*[citiesNumber];
	for (int i = 0; i < citiesNumber; i++)
	{
		this->nodeMatrix[i] = new int[citiesNumber];
	}

	for (int i = 0; i < citiesNumber; i++)
	{
		for (int j = 0; j < citiesNumber; j++)
		{
			this->nodeMatrix[i][j] = 0;
		}
	}

	for (int i = 0; i < citiesNumber; i++)
	{
		for (int j = 0; j < citiesNumber; j++)
		{
			this->nodeMatrix[i][j] = parentMatrix[i][j];
		}
	}

	for (int x = 0; level != 0 && x < citiesNumber; x++)
	{
		this->nodeMatrix[i][x] = INT_MAX;
		this->nodeMatrix[x][j] = INT_MAX;
	}

	this->nodeMatrix[j][0] = INT_MAX;

	this->graphLevel = level;
	this->currentCityNumber = j;
	this->citiesNumber = citiesNumber;
}

Node::~Node()
{
	//Free each sub-array
	for (int i = 0; i < this->citiesNumber; ++i) {
		delete[] this->nodeMatrix[i];
	}
	//Free the array of pointers
	delete[] this->nodeMatrix;
}
