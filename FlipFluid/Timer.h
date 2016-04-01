#ifndef TIMER_H
#define TIMER_H
//------------------------------------------------------------------------------
//My Timer. Use for calculate time usage for per frame.
//
//------------------------------------------------------------------------------


#include <iostream>
#include <windows.h>
#include <winnt.h>
#include <vector>
#include <map>
#include <string>

using namespace std;
class Timer
{
private:
	static __int64 staTime;
	static __int64 preTime;
	static __int64 curTime;
	static __int64 totalTime;
	static double secsPerCnt;
	static __int64 numFrames;
	static map<string, __int64> timingBreakdown;

public:
	static void start()
	{
		__int64 cntsPerSec = 0;
		QueryPerformanceFrequency((LARGE_INTEGER*)&cntsPerSec);
		secsPerCnt = 1.0f / cntsPerSec;
		QueryPerformanceCounter((LARGE_INTEGER*)&staTime);
		numFrames = 0;
	}
	inline static void AddFrame()
	{
		numFrames++;
	}

	//record the time now.
	inline static void tic()													//tic() and toc() should be used together
	{																			//like
		QueryPerformanceCounter((LARGE_INTEGER*)&preTime);						// tic();
	}																			//     do something
	//record the time used by input parameter "name".							// toc("something");
	inline static void toc(string name)
	{
		QueryPerformanceCounter((LARGE_INTEGER*)&curTime);
		timingBreakdown[name] += (curTime - preTime);
	}
	// print the time usage of your record.
	static void printTime()
	{
		double sum = 0.0;
		cout << "========================================================================================" << endl;
		for (map<string, __int64>::iterator i = timingBreakdown.begin(); i != timingBreakdown.end(); i++)
		{
			cout << i->first << "                         " << (double)i->second * secsPerCnt / numFrames << "s / frame" <<endl;
			sum += (double)i->second / totalTime;
			cout << "Persent: " << 100.0 * (double)i->second / totalTime << endl;
		}
		cout << "mic: " << endl;
		cout << "Persent: " << 100.0 - 100.0 * sum << endl;
		cout << "========================================================================================" << endl;
	}
	//get fps
	inline static double fps()
	{
		QueryPerformanceCounter((LARGE_INTEGER*)&curTime);
		totalTime = curTime - staTime;
		return numFrames / (totalTime * secsPerCnt);
	}
};


#endif