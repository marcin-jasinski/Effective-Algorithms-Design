#pragma once
#include <vector>

// This class represents a node in a graph used in Branch And Bound calculations.
class Node
{
public:
	std::vector<int> pathToVertex;
	int** nodeMatrix;
	int lowerBound;
	int currentCityNumber;
	int graphLevel;
	int citiesNumber;

	Node();
	~Node();

	Node(int **parentMatrix, std::vector<int> const &path, int level, int i, int j, int citiesNumber);
};

