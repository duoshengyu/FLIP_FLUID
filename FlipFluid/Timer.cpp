#include "Timer.h"

__int64 Timer :: staTime;
__int64 Timer :: preTime;
__int64 Timer :: curTime;
__int64 Timer :: totalTime;
double  Timer :: secsPerCnt;
__int64 Timer :: numFrames;
map<string, __int64> Timer :: timingBreakdown;