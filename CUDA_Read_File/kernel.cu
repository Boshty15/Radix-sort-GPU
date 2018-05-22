#pragma once
#ifdef __INTELLISENSE__
void __syncthreads();
#endif

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <device_functions.h>
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <stdio.h>
#include "kernel.h"
#include <iostream>     
#include <fstream>
#include <vector>

#define SIZE 32 //512
#define BIT 32

using namespace std;

__device__ unsigned int device_data[SIZE];
__device__ int ddata_s[SIZE];

//Test v gpu in nazaj
__global__ void gpucopy(int* src, int* dst)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;;
	dst[i] = src[i];
}

__global__ void radixSort(int* sort_tmp, const int num_lists, const int num_element, int* sort_tmp_0, int* sort_tmp_1) {
	int tid = (blockIdx.x * blockDim.x) + threadIdx.x;
	//printf("\n%d num list: %d" ,tid, num_lists );
	//printf("\n num element: %d" , num_element);
	if (sort_tmp == NULL) {
		printf("\nsort_tmp is null");
	}
	
	
	for (int bit = 0; bit < BIT; bit++) {
		int base_cnt_0 = 0, base_cnt_1 = 0;
		//printf("\ntest1 + %d",tid );
		for (int i = 0; i < num_element; i+=num_lists) {
			int elem = sort_tmp[i + tid];
			const int bit_mask = (1 << bit);
			//printf("\ntest1 + %d", elem);
			if ((elem & bit_mask) > 0) {
				sort_tmp_1[base_cnt_1 + tid] = elem;
				base_cnt_1 += num_lists;
			}
			else {
				sort_tmp_0[base_cnt_0 + tid] = elem;
				base_cnt_0+=num_lists;
			}
		}
		//copy data back to source - first the zero list
		for (int i = 0; i < base_cnt_0; i += num_lists) {
			sort_tmp[i + tid] = sort_tmp_0[i + tid];
		}
		//copy data back to source - then the one list
		for (int i = 0; i < base_cnt_1; i += num_lists) {
			sort_tmp[base_cnt_0 + i + tid] = sort_tmp_1[i + tid];
		}
	}
	
	__syncthreads();
	/*for (int i = 0; i < num_element; i++) {
		printf("\n%d", sort_tmp[i]);
	}*/
}


void WriteFile(vector<double> vector) {
	ofstream myfile;
	myfile.open("out.txt");
	for each (double var in vector)
	{
		myfile << var << " ";
	}
	myfile.close();
}
vector<int> ReadFile(string filepath) {
	vector<int> vector;
	ifstream inputFile(filepath);
	double value;

	// read the elements in the file into a vector  
	while (inputFile >> value) {
		vector.push_back(value * 100);
	}
	for each (int var in vector)
	{
		cout << var << " ";
	}
	inputFile.close();
	return vector;
}

int main(int argc, char * argv[])
{
	cout << argv[1] << endl;

	std::ifstream fin(argv[1]);
	
	std::vector<int> host_double(7);
	std::vector<int> host_double2(7);

	float totalTime = 0;
	float milliseconds = 0;

	//Time
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	if (!fin) {
		cerr << "Datoteka ne obstaja " << argv[1] << endl;
	}
	else {
		cout << "Datoteka obstaja " << argv[1] << endl;
		host_double = ReadFile(argv[1]);

		//WriteFile(host_double);
	}
	int size = host_double.size();
	
	int* device = NULL;
	int* tmp_0 = NULL;
	int* tmp_1 = NULL;
	
	cudaMalloc((void**)&device, size * sizeof(int));
	cudaMalloc((void**)&tmp_0, size * sizeof(int));
	cudaMalloc((void**)&tmp_1, size * sizeof(int));

	int* a = host_double.data();
	

	//num_list 32 

	cudaMemcpy(device,&*host_double.data(), size * sizeof(int), cudaMemcpyHostToDevice);
	
	radixSort << < 1,SIZE >> >(device, SIZE, size, tmp_0, tmp_1);
	//ParalelRadixSort << <(size + size - 1) / size, size >> >();
	//gpucopy << <1,SIZE>> >(device, device2);
	
	cudaDeviceSynchronize();

	cudaMemcpy(&*host_double.data(), device, size * sizeof(int), cudaMemcpyDeviceToHost);

	cout << endl << "host data" << endl;
	for each (int var in host_double)
	{
		double d = (double)var / 100;
		cout << d << " " << endl;
	}

	cudaFree(device);
	cudaFree(tmp_0);
	cudaFree(tmp_1);

	//std::vector<int> v(10);

	//unsigned int hdata[SIZE];


	//	int size = host_double.size();
	//	//vector.resize(size);

	//	std::copy(host_double.begin(), host_double.end(), hdata);
	//	//int *a = vector.data();

	//	// Copy data from host to device
	//	cudaMemcpyToSymbol(device_data, hdata, SIZE * sizeof(unsigned int));

	//	// Execution time measurement, that point starts the clock
	//	cudaEventRecord(start);
	//	ParalelRadixSort << <1, SIZE>> >();

	//	// Execution time measurement, that point stops the clock
	//	cudaEventRecord(stop);

	//	// Make kernel function synchronous
	//	cudaDeviceSynchronize();

	//	cudaEventElapsedTime(&milliseconds, start, stop);
	//	totalTime += milliseconds;
	//	//milliseconds = 0;

	//	// Copy data from device to host
	//	cudaMemcpyFromSymbol(hdata, device_data, SIZE * sizeof(unsigned int));

	//	/*std::vector<int> v(hdata, hdata + sizeof hdata / sizeof hdata[0]);
	//	v.resize(size);*/
	//	milliseconds = 0;
	//

	///*printf("TIME %4.2fs", milliseconds);
	//printf("\n");
	//printf("Effective Bandwidth (GB/s): %fn", SIZE * 4 * 3 / milliseconds / 1e6);
	//printf("\n");*/
	//cout << endl;
	//cout << "Host Data" << endl;
	///*for each (int var in v)
	//{
	//	double d = (double)var / 100;
	//	cout << d << " " << endl;
	//}*/
	//for (int i = 0; i < SIZE; i++) {
	//	double d = (double)hdata[i] / 100;
	//	cout << d << endl;
	//}
	//
	//cout << endl;
	//cout << "TIME %4.2fs " << totalTime << endl;

    return 0;
}

