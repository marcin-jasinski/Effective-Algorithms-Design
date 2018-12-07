#include "pch.h"
#include "Stopwatch.h"
#include <iostream>

Stopwatch::Stopwatch()
{
}

Stopwatch::~Stopwatch()
{
}

void Stopwatch::StartCounter()
{
	LARGE_INTEGER li;
	if (!QueryPerformanceFrequency(&li)) std::cout << "QueryPerformanceFrequency failed!\n";

	PCFreq = double(li.QuadPart);

	QueryPerformanceCounter(&li);
	CounterStart = li.QuadPart;
}

double Stopwatch::GetCounter()
{
	LARGE_INTEGER li;
	QueryPerformanceCounter(&li);
	return double(li.QuadPart - CounterStart) / PCFreq;
}
