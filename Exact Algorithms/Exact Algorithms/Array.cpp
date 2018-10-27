#include "pch.h"
#include "Array.h"
#include <fstream>
#include <string>
#include <iostream>
#include <conio.h>

// default constructor of an array
Array::Array()
{
	this->arraySize = 0;
	this->_headPtr = nullptr;
}

// constructor creating an array of a specified size
// assigns new dynamicly allocated array of integers to _headPtr
Array::Array(int arraySize)
{
	this->arraySize = arraySize;
	this->_headPtr = new int[arraySize];
}

// default destructor
// releasing memory 
Array::~Array()
{
	this->_headPtr = nullptr;
	this->arraySize = 0;
}

// returns current size (number of elements in the array)
int Array::getSize() const
{
	return this->arraySize;
}

// returns pointer to the first element in an array (or nullptr if array was created with default constructor)
int* Array::getHeadPtr()
{
	return this->_headPtr;
}

// reading data from a text file "testData.txt"
// first line sets array size (number of elements)
void Array::readDataFromFile()
{
	std::fstream file;
	file.open("testData.txt", std::ios::in);
	if (file.good() == true)
	{	
		std::cout << "\nFile acces granted." << std::endl;
		std::string input;
		getline(file, input);

		this->arraySize = std::stoi(input);
		this->_headPtr = new int[arraySize];	// creating new integer array

		for (int i = 0; i < this->arraySize && !file.eof(); i++) {
			input.clear();
			std::getline(file, input);
			// std::cout << "Input " << i << " = " << input << std::endl;
			this->_headPtr[i] = std::stoi(input);
		}

		file.close();
	}

	else std::cout << "Error opening file!!!" << std::endl;
}

// reading data from keyboard input
void Array::readDataFromKeyboard()
{
	std::cout << "\nSet array size: ";
	int userSize;
	std::cin >> userSize;
	this->arraySize = userSize;
	this->_headPtr = new int[arraySize];	// creating new integer array

	int userInput;
	for (int i = 0; i < userSize; i++) {
		std::cout << "Value at index [" << i << "] : ";
		std::cin >> userInput;
		this->_headPtr[i] = userInput;
	}
}

// inserts en element on the beginning of the array
void Array::pushFront(int element)
{
	//if the array is empty
	if (this->arraySize == 0) {
		this->_headPtr = new int[1];	// create new one-element array
		this->_headPtr[0] = element;
		this->arraySize++;
	}
	else {
		int* _tempPtr = new int[this->arraySize + 1];						// temporal "buffer" array for holding already existing elements
		memcpy(_tempPtr + 1, this->_headPtr, arraySize * sizeof(int));		// copy elements from existing array to buffer shifted by one index up
		delete[] _headPtr;													// free memory currently occupied by old array (elements are safely copied to buffer)
		_tempPtr[0] = element;												// place new element on the beginning of buffer array
		this->_headPtr = _tempPtr;											// assign head pointer to buffer array
		this->arraySize++;		
	}
}

// inserts an element on the end of the array
void Array::pushBack(int element)
{
	// if the array is empty
	if (this->arraySize == 0) {
		this->_headPtr = new int[1];	// create new one-element array
		*_headPtr = element;
		this->arraySize++;
	}
	else {
		int* _tempPtr = new int[arraySize + 1];								// same thing as above, but this time elements is inserted on the end of array
		memcpy(_tempPtr, _headPtr, arraySize * sizeof(int));
		delete[] _headPtr;
		_tempPtr[arraySize] = element;
		_headPtr = _tempPtr;
		this->arraySize++;
	}
}
 
// deletes first array element
void Array::popFront()
{
	if (this->arraySize == 0) {
		std::cout << "Array is empty." << std::endl;
		return;
	}

	deleteValueFromIndex(0);
}

// deletes last array element
void Array::popBack()
{
	if (this->arraySize == 0) {
		std::cout << "Array is empty." << std::endl;
		return;
	}

	deleteValueFromIndex(this->arraySize - 1);
}

// relocates array with a new size and places new element on a specified index
// elements originally placed after selected index are shifted by one index number up
void Array::insertValueOnIndex(int index, int element)
{
	// Index out of bounds 
	if (index < 0 || index >= arraySize)
	{
		std::cout << "Index out of bounds." << std::endl;
		return;
	}

	// there is no point of inserting a value to a non-existant array without explicit intention
	if (this->arraySize == 0) return;
	
	int* _tempPtr = new int[this->arraySize + 1];				 // temporal "buffer" array for holding already existing elements
	memcpy(_tempPtr, this->_headPtr, index * sizeof(int));		 // copying to buffer only elements originally being "above" new element (thus index*sizeof(int))
	_tempPtr[index] = element;
	memcpy(_tempPtr + index + 1, this->_headPtr + index, (this->arraySize - index) * sizeof(int)); // placing the rest of the elements back with one index up
	delete[] _headPtr;
	this->_headPtr = _tempPtr;	// assigning head pointer back to point on array
	this->arraySize++;
}

// deletes element from the specified index and relocates array with a new size
// elements originally placed after selected index are shifted by one index number down
void Array::deleteValueFromIndex(int index)
{
	// Index out of bounds 
	if (index < 0 || index >= arraySize)
	{
		std::cout << "Index out of bounds." << std::endl;
		return;
	}

	// there is no point of deleting a value from a non-existant array
	if (arraySize == 0) {
		std::cout << "Array is empty." << std::endl;
		return;
	}

	int* _tempPtr = new int[arraySize - 1];						// temporal "buffer" array for holding already existing elements
	memcpy(_tempPtr, _headPtr, index * sizeof(int));
	this->arraySize--;
	memcpy(_tempPtr + index, this->_headPtr + index + 1, (this->arraySize - index) * sizeof(int));
	delete[] _headPtr;
	this->_headPtr = _tempPtr;
}

// replacing an old value at index "index" with a new one "element"
void Array::replaceValueOnIndex(int index, int element)
{
	// Index out of bounds 
	if (index < 0 || index >= arraySize)
	{
		std::cout << "Index out of bounds." << std::endl;
		return;
	}

	this->_headPtr[index] = element;
}

// returns true if array contains specified value
void Array::findValue(int element)
{
	for (int i = 0; i < this->arraySize; i++) {
		if (*(this->_headPtr + i) == element) {			// if determined value matches with the one in the array
			std::cout << "\nElement found on position " << i << std::endl;
			return;
		}
	}
	std::cout << "\nElement not found." << std::endl;
	return;
}

int Array::get(int index)
{
	return *(_headPtr + index);
}

// overloaded [] operator for array-like element acces
int Array::operator[](int index) const
{
	return *(_headPtr+index);
}

// overloaded operator for writing array contents to the output stream
std::ostream & operator<<(std::ostream& out, Array& array)
{
	if (array.arraySize == 0) {
		out << "[empty]" << std::endl;
		return out;
	}

	out << "[";
	for (int i = 0; i < array.getSize(); i++) {
		out << array[i];
		if (i == array.getSize() - 1) out << "]\n";
		else out << ",";
	}
	return out;
}
