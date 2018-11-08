#pragma once

class Node
{
private:

	int treeLevel;
	int cityNumber;
	int estimatedLowerBound;

	int** reducedNodeMatrix;

public:
	Node();
	Node(int**, int, int, int, int);
	~Node();

	int getTreeLevel();
	int getCityNumber();
	int getEstimatedLowerBound();
	int** getReducedNodeMatrix();
};

