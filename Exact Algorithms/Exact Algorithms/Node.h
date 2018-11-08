#pragma once

class Node
{
private:

	int treeLevel;
	int cityNumber;
	
public:

	int estimatedLowerBound;
	int** reducedNodeMatrix;

	Node();
	Node(int**, int, int, int, int);
	~Node();

	int getTreeLevel();
	int getCityNumber();
	int getEstimatedLowerBound();
	int** getReducedNodeMatrix();

	void setEstimatedLowerBound(int);
};

