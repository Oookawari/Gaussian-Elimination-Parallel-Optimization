#include <arm_neon.h>
#include <iostream>
#include <iomanip>
#include <random>
#include <stdio.h>
#include <time.h>

using namespace std;
const int ORDER = 100;
const int TIME = 100;
float matrix[ORDER][ORDER + 1];
float temp[ORDER][ORDER + 1];
int main() {
	unsigned seed = 47;
	srand(seed);
	for (int i = 0; i < ORDER; i++)
		for (int j = 0; j <= ORDER; j++)
			temp[i][j] = rand() % 10 + rand() % 10000 * 0.0001;
	struct timespec sts, ets;
	double time_sum = 0.0;
	for (int time = 0; time < TIME; time++) {
		//初始化/恢复矩阵 不计入时间
		for (int i = 0; i < ORDER; i++)
			for (int j = 0; j <= ORDER; j++)
				matrix[i][j] = temp[i][j];

		timespec_get(&sts, TIME_UTC);

		for (int k = 0; k < ORDER; k++) {
			float32x4_t vector_k_k = vdupq_n_f32(matrix[k][k]);
			int j;
			for (j = k + 1; j + 4 <= ORDER + 1; j += 4) {
				float32x4_t vector_k_j = vld1q_f32(&matrix[k][j]);
				vector_k_j = vdivq_f32(vector_k_j, vector_k_k);
				vst1q_f32(&matrix[k][j], vector_k_j);
			}
			for (; j <= ORDER; j++) {
				matrix[k][j] = matrix[k][j] / matrix[k][k];
			}
			matrix[k][k] = 1.0;
			for (int i = k + 1; i < ORDER; i++) {
				float32x4_t vector_i_k = vdupq_n_f32(matrix[i][k]);
				for (j = k + 1; j + 4 <= ORDER + 1; j += 4) {
					float32x4_t vector_k_j = vld1q_f32(&matrix[k][j]);
					float32x4_t vector_i_j = vld1q_f32(&matrix[i][j]);
					float32x4_t vector_temp = vmulq_f32(vector_i_k, vector_k_j);
					vector_i_j = vsubq_f32(vector_i_j, vector_temp);
					vst1q_f32(&matrix[i][j], vector_i_j);
				}
				for (; j <= ORDER; j++) {
					matrix[i][j] -= matrix[i][k] * matrix[k][j];
				}
				matrix[i][k] = 0.0;
			}
		}
		timespec_get(&ets, TIME_UTC);
		time_sum += (ets.tv_nsec - sts.tv_nsec) * 0.000001 + (ets.tv_sec - sts.tv_sec) * 1000;
	}

	
	cout << fixed << setprecision(4) << time_sum << endl;
	/*
	for (int i = 0; i < ORDER; i++) {
		for (int j = 0; j <= ORDER; j++) {
			cout << fixed << setprecision(4) << matrix[i][j] << " ";
		}
		cout << endl;
	}
	*/
	return 0;
}
