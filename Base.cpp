#include <fstream>
#include <iostream>
#include <Windows.h>
#include "Base.h"

long Base::GetElapsedProcessTime()
{
	FILETIME createTime;
	FILETIME exitTime;
	FILETIME kernelTime;
	FILETIME userTime;

	long ElapsedTime;
	if (GetProcessTimes(GetCurrentProcess(), &createTime, &exitTime, &kernelTime, &userTime) != 0)//����GetProcessTimesϵͳ����
	{
		SYSTEMTIME userSystemTime;
		if (FileTimeToSystemTime(&userTime, &userSystemTime) != -1)//FileTimeToSystemTimeʱ���ʽת��
			ElapsedTime = (userSystemTime.wDay - 1) * 24 * 3600 * 1000
			+ userSystemTime.wHour * 3600 * 1000 +
			userSystemTime.wMinute * 60 * 1000 +
			userSystemTime.wSecond * 1000 +
			userSystemTime.wMilliseconds;
		else
			ElapsedTime = 0;
	}
	else
		ElapsedTime = 0;
	return ElapsedTime;
}


