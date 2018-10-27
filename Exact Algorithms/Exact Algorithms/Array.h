#pragma once
#include <iostream>

class Array
{
private:
	int* _headPtr;	// pointer to dynamic array of integers
	int arraySize;  // current number of elements in the array

public:
	Array();
	Array(int);
	~Array();

	int getSize() const;
	int* getHeadPtr();

	void readDataFromFile();
	void readDataFromKeyboard();

	void pushFront(int);
	void pushBack(int);
	void popFront();
	void popBack();

	void insertValueOnIndex(int, int);
	void deleteValueFromIndex(int);
	void replaceValueOnIndex(int, int);
	void findValue(int);

	int get(int);
	int operator[](int) const;

	friend std::ostream& operator << (std::ostream&, Array&);
};

