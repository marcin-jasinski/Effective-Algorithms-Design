#pragma once
#include <cstdlib>
#include <ctime>
#include <windows.h>

class Stopwatch
{
private:
	LARGE_INTEGER li;
	double PCFreq = double(li.QuadPart) / 1000000.0;
	__int64 CounterStart = 0;

public:
	Stopwatch();
	~Stopwatch();

	void StartCounter();
	double GetCounter();
};