#include <iostream>
#include <windows.h>
using namespace std;
const int SIZE_ = 10;
const int RE = 50000000;
int a[SIZE_];
int A[SIZE_];
int main()
{
    cout << "数组大小: " << SIZE_ << " 重复次数: " << RE << endl;
    LARGE_INTEGER nFreq;
    LARGE_INTEGER nBeginTime;
    LARGE_INTEGER nEndTime;
    double time, time_op, time_op2;
    QueryPerformanceFrequency(&nFreq);
    //Initializing
    for(int i = 0; i < SIZE_; i++){
        a[i] = i;
        A[i] = i;
    }


    //Common algorithm:
    QueryPerformanceCounter(&nBeginTime);
    for(int k = 0; k < RE; k++){
        int sum = 0;
        for (int i = 0; i < SIZE_; i++){
            sum += a[i];
        }
    }
    QueryPerformanceCounter(&nEndTime);
    time = (double)(nEndTime.QuadPart - nBeginTime.QuadPart) / (double)nFreq.QuadPart;
    cout << "General algorithm time consumption:   ";
    printf("%f\n", time);

    //Parallel algorithm(dual-link):
    QueryPerformanceCounter(&nBeginTime);
    for(int k = 0; k < RE; k++){
        int sum1 = 0, sum2 = 0, sum = 0;
        for (int i = 0; i < SIZE_; i += 2) {
            sum1 += A[i];
            sum2 += A[i + 1];
        }sum = sum1 + sum2;
    }

    QueryPerformanceCounter(&nEndTime);
    time_op = (double)(nEndTime.QuadPart - nBeginTime.QuadPart) / (double)nFreq.QuadPart;
    cout << "Optimized algorithm time consumption:   ";
    printf("%f\n", time_op);

    //Parallel algorithm(recursion):
    QueryPerformanceCounter(&nBeginTime);
    for(int k = 0; k < RE; k++){
        for (int m = SIZE_; m > 1; m /= 2){
            if(m % 2 != 0) {a[m] = 0;m ++;}
            for (int i = 0; i < m / 2; i++){
                a[i] = a[i * 2] + a[i * 2 + 1];
            }
        }
    }

    QueryPerformanceCounter(&nEndTime);
    time_op2 = (double)(nEndTime.QuadPart - nBeginTime.QuadPart) / (double)nFreq.QuadPart;
    cout << "Optimized algorithm2 time consumption:   ";
    printf("%f\n", time_op2);

}
