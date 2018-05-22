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
#include <bitset>
//
//#define SIZE 5
//#define BIT 5
//
using namespace std;
//
//__device__ unsigned int device_data[SIZE];
//__device__ int ddata_s[SIZE];
//
////Test v gpu in nazaj
//__global__ void gpucopy(int* src, int* dst)
//{
//	int i = (blockIdx.x * blockDim.x) + threadIdx.x;
//	dst[i] = src[i];
//}
//
//__global__ void radixSort(int* sort_tmp, const int num_lists, const int num_element, int* sort_tmp_1) {
//	int tid = (blockIdx.x * blockDim.x) + threadIdx.x;
//	//printf("\n%d num list: %d" ,tid, num_lists );
//	//printf("\n num element: %d" , num_element);
//	if (sort_tmp == NULL) {
//		printf("\nsort_tmp is null");
//	}
//	
//	
//	for (int bit = 0; bit < BIT; bit++) {
//		int base_cnt_0 = 0, base_cnt_1 = 0;
//		//printf("\ntest1 + %d",tid );
//		for (int i = 0; i < num_element; i+=num_lists) {
//			int elem = sort_tmp[i + tid];
//			const int bit_mask = (1 << bit);
//			//printf("\ntest %d + %d bit %d",tid, elem, bit_mask);
//			if ((elem & bit_mask) > 0) {
//				sort_tmp_1[base_cnt_1 + tid] = elem;
//				base_cnt_1 += num_lists;
//				//printf("\n element 1 %d = %d",i, elem);
//			}
//			else {
//				sort_tmp[base_cnt_0 + tid] = elem;
//				base_cnt_0+=num_lists;
//				//printf("\n element 0 %d = %d",i, elem);
//			}
//		}
//		//copy data back to source - first the zero list
//		for (int i = 0; i < base_cnt_1; i += num_lists) {
//			sort_tmp[base_cnt_0 + i + tid] = sort_tmp_1[i + tid];
//			//printf("\nsort tmp: %d", sort_tmp_0[i + tid]);
//		}
//	}
//	
//	__syncthreads();
//
//	for (int i = 0; i < num_element; i++) {
//		printf("\n%d", sort_tmp[i]);
//	}
//}
//
//
//void WriteFile(vector<double> vector) {
//	ofstream myfile;
//	myfile.open("out.txt");
//	for each (double var in vector)
//	{
//		myfile << var << " ";
//	}
//	myfile.close();
//}
vector<unsigned int> ReadFile(string filepath) {
	vector<unsigned int> vector;
	ifstream inputFile(filepath);
	double value;

	// read the elements in the file into a vector  
	while (inputFile >> value) {
		//string binary = bitset<16>(value * 100).to_string();
		vector.push_back(value * 100);
	}
	for each (unsigned int var in vector)
	{
		cout << var << " " << endl;
	}
	inputFile.close();
	return vector;
}
//
//int main(int argc, char * argv[])
//{
//	cout << argv[1] << endl;
//
//	std::ifstream fin(argv[1]);
//	
//	std::vector<string> host_double(7);
//	std::vector<int> host_double2(7);
//
//	float totalTime = 0;
//	float milliseconds = 0;
//
//	//Time
//	cudaEvent_t start, stop;
//	cudaEventCreate(&start);
//	cudaEventCreate(&stop);
//
//	if (!fin) {
//		cerr << "Datoteka ne obstaja " << argv[1] << endl;
//	}
//	else {
//		cout << "Datoteka obstaja " << argv[1] << endl;
//		host_double = ReadFile(argv[1]);
//
//		//WriteFile(host_double);
//	}
//	int size = host_double.size();
//	
//	int* device = NULL;
//	int* tmp_0 = NULL;
//	int* tmp_1 = NULL;
//	
//	cudaMalloc((void**)&device, size * sizeof(int));
//	cudaMalloc((void**)&tmp_0, size * sizeof(int));
//	cudaMalloc((void**)&tmp_1, size * sizeof(int));
//
//	string s = host_double[0];
//	cout << "test: " << s[0] << endl;
//	
//	
//	//string binary = bitset<16>(222).to_string(); //to binary
//	//cout << binary << "\n";	
//
//	//unsigned long decimal = bitset<16>(binary).to_ulong();
//	//cout << decimal << "\n";
//
//	//num_list 32 
//
//	cudaMemcpy(device,&*host_double.data(), size * sizeof(int), cudaMemcpyHostToDevice);
//	
//	//radixSort << < 1,5>> >(device, SIZE, size, tmp_1);
//	//ParalelRadixSort << <(size + size - 1) / size, size >> >();
//	//gpucopy << <1,SIZE>> >(device, device2);
//	
//	cudaDeviceSynchronize();
//
//	cudaMemcpy(&*host_double.data(), device, size * sizeof(int), cudaMemcpyDeviceToHost);
//
//	//cout << endl << "host data" << endl;
//	//for each (string var in host_double)
//	//{
//	//	double d = (double)var / 100;
//	//	cout << d << " " << endl;
//	//}
//
//	cudaFree(device);
//	cudaFree(tmp_0);
//	cudaFree(tmp_1);
//
//	//std::vector<int> v(10);
//
//	//unsigned int hdata[SIZE];
//
//
//	//	int size = host_double.size();
//	//	//vector.resize(size);
//
//	//	std::copy(host_double.begin(), host_double.end(), hdata);
//	//	//int *a = vector.data();
//
//	//	// Copy data from host to device
//	//	cudaMemcpyToSymbol(device_data, hdata, SIZE * sizeof(unsigned int));
//
//	//	// Execution time measurement, that point starts the clock
//	//	cudaEventRecord(start);
//	//	ParalelRadixSort << <1, SIZE>> >();
//
//	//	// Execution time measurement, that point stops the clock
//	//	cudaEventRecord(stop);
//
//	//	// Make kernel function synchronous
//	//	cudaDeviceSynchronize();
//
//	//	cudaEventElapsedTime(&milliseconds, start, stop);
//	//	totalTime += milliseconds;
//	//	//milliseconds = 0;
//
//	//	// Copy data from device to host
//	//	cudaMemcpyFromSymbol(hdata, device_data, SIZE * sizeof(unsigned int));
//
//	//	/*std::vector<int> v(hdata, hdata + sizeof hdata / sizeof hdata[0]);
//	//	v.resize(size);*/
//	//	milliseconds = 0;
//	//
//
//	///*printf("TIME %4.2fs", milliseconds);
//	//printf("\n");
//	//printf("Effective Bandwidth (GB/s): %fn", SIZE * 4 * 3 / milliseconds / 1e6);
//	//printf("\n");*/
//	//cout << endl;
//	//cout << "Host Data" << endl;
//	///*for each (int var in v)
//	//{
//	//	double d = (double)var / 100;
//	//	cout << d << " " << endl;
//	//}*/
//	//for (int i = 0; i < SIZE; i++) {
//	//	double d = (double)hdata[i] / 100;
//	//	cout << d << endl;
//	//}
//	//
//	//cout << endl;
//	//cout << "TIME %4.2fs " << totalTime << endl;
//
//    return 0;
//}
//


#define WSIZE 32
#define LOOPS 1
#define UPPER_BIT 31
#define LOWER_BIT 0

__device__ unsigned int ddata[WSIZE];

// naive warp-level bitwise radix sort

__global__ void mykernel() {
	__shared__ volatile unsigned int sdata[WSIZE * 2];
	// load from global into shared variable
	int tid = (blockIdx.x * blockDim.x) + threadIdx.x;
	sdata[tid] = ddata[tid];
	unsigned int bitmask = 1 << LOWER_BIT;
	unsigned int offset = 0;
	unsigned int thrmask = 0xFFFFFFFFU << tid;
	unsigned int mypos;
	//  for each LSB to MSB
	for (int i = LOWER_BIT; i <= UPPER_BIT; i++) {
		unsigned int mydata = sdata[((WSIZE - 1) - tid) + offset];
		unsigned int mybit = mydata&bitmask;
		// get population of ones and zeroes (cc 2.0 ballot)
		unsigned int ones = __ballot(mybit); // cc 2.0
		unsigned int zeroes = ~ones;
		offset ^= WSIZE; // switch ping-pong buffers
						 // do zeroes, then ones
		if (!mybit) // threads with a zero bit
					// get my position in ping-pong buffer
			mypos = __popc(zeroes&thrmask);
		else        // threads with a one bit
					// get my position in ping-pong buffer
			mypos = __popc(zeroes) + __popc(ones&thrmask);
		// move to buffer  (or use shfl for cc 3.0)
		sdata[mypos - 1 + offset] = mydata;
		// repeat for next bit
		bitmask <<= 1;
	}
	// save results to global
	ddata[tid] = sdata[tid + offset];
}

int main(int argc, char * argv[]) {
	std::ifstream fin(argv[1]);
	std::vector<unsigned int> host_double(WSIZE);

	unsigned int hdata[WSIZE];

	

	//if (!fin) {
	//	cerr << "Datoteka ne obstaja " << argv[1] << endl;
	//}
	//else {
	//	cout << "Datoteka obstaja " << argv[1] << endl;
	//	host_double = ReadFile(argv[1]);
	//	
	//	//WriteFile(host_double);
	//}
	//std::copy(host_double.begin(), host_double.end(), hdata);

	for (int lcount = 0; lcount < LOOPS; lcount++) {
		unsigned int range = 1U << UPPER_BIT;
		for (int i = 0; i < WSIZE; i++) hdata[i] = rand() % range;
		cudaMemcpyToSymbol(ddata, hdata, WSIZE * sizeof(unsigned int));
		mykernel << <1, 32 >> >();
		cudaMemcpyFromSymbol(hdata, ddata, WSIZE * sizeof(unsigned int));
		for (int i = 0; i < WSIZE - 1; i++) if (hdata[i] > hdata[i + 1]) { printf("sort error at loop %d, hdata[%d] = %d, hdata[%d] = %d\n", lcount, i, hdata[i], i + 1, hdata[i + 1]); return 1; }
		 printf("sorted data:\n");
		for (int i = 0; i < WSIZE; i++) cout << (double)hdata[i] / 100 << endl;
	}
	printf("Success!\n");
	return 0;
}
