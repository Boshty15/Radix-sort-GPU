


//#pragma once
//#ifdef __INTELLISENSE__
//void __syncthreads();
//#endif


#include "cuda_runtime.h"
#include "device_launch_parameters.h"

//#ifndef __CUDACC_RTC__ 
//#define __CUDACC_RTC__
//#endif
//#pragma once
//#ifdef __INTELLISENSE__
//void __syncthreads();
//
//#endif

#include <cuda.h>
#include <cuda_runtime_api.h>
#include <stdio.h>
#include "kernel.h"
#include <iostream>     
#include <fstream>
#include <vector>
#include <bitset>
#include <time.h>
#include <chrono>
#include "device_functions.h"
#include <math.h> 


//#define SIZE 32
#define BIT 32
//
using namespace std;
using namespace chrono;
//
//__device__ unsigned int device_data[SIZE];
//__device__ int ddata_s[SIZE];
//
////Test v gpu in nazaj
//__global__ void gpucopy(int* src, int* dst)
//{
//	int index = blockIdx.x * blockDim.x + threadIdx.x;
//	int stride = blockDim.x * gridDim.x;
//	for (int i = index; i < SIZE; i += stride) {
//		if(src[i] > 50)
//			src[i] = dst[i] + 100;
//		else 
//			src[i] = dst[i] - 100;
//	}
//		
//}
//

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

//
//int main(int argc, char * argv[])
//{
//	cout << argv[1] << endl;
//
//	std::ifstream fin(argv[1]);
//	
//	std::vector<unsigned int> host_double(7);
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
//	int blockSize = 4;
//	int numBlocks = (SIZE + blockSize - 1) / blockSize;
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
//	//radixSort << <numBlocks, blockSize>> >(device, SIZE, size, tmp_1);
//	//ParalelRadixSort << <(size + size - 1) / size, size >> >();
//	gpucopy << <numBlocks, blockSize>> >(device, device);
//	
//	cudaDeviceSynchronize();
//
//	cudaMemcpy(&*host_double.data(), device, size * sizeof(int), cudaMemcpyDeviceToHost);
//
//	cout << endl << "host data" << endl;
//	for each (unsigned int var in host_double)
//	{
//		double d = (double)var / 100;
//		cout << d << " " << endl;
//	}
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



#define SIZE (33 * 1024)
#define LOOPS 1
#define UPPER_BIT 31
#define LOWER_BIT 0
int getMax(unsigned int* arr, int n)
{
	int mx = arr[0];
	for (int i = 1; i < n; i++)
		if (arr[i] > mx)
			mx = arr[i];
	return mx;
}

// A function to do counting sort of arr[] according to
// the digit represented by exp.
void countSort(unsigned int* arr, int n, int exp)
{
// output array
	vector<unsigned int> ve(n);
	unsigned int* output = ve.data();
	int i, count[10] = { 0 };

	// Store count of occurrences in count[]
	for (i = 0; i < n; i++)
		count[(arr[i] / exp) % 10]++;

	// Change count[i] so that count[i] now contains actual
	//  position of this digit in output[]
	for (i = 1; i < 10; i++)
		count[i] += count[i - 1];

	// Build the output array
	for (i = n - 1; i >= 0; i--)
	{
		output[count[(arr[i] / exp) % 10] - 1] = arr[i];
		count[(arr[i] / exp) % 10]--;
	}

	// Copy the output array to arr[], so that arr[] now
	// contains sorted numbers according to current digit
	for (i = 0; i < n; i++)
		arr[i] = output[i];
}

// The main function to that sorts arr[] of size n using 
// Radix Sort
void serialRadixSort(unsigned int* arr, int n)
{
	// Find the maximum number to know number of digits
	int m = getMax(arr, n);

	// Do counting sort for every digit. Note that instead
	// of passing digit number, exp is passed. exp is 10^i
	// where i is current digit number
	for (int exp = 1; m / exp > 0; exp *= 10)
		countSort(arr, n, exp);
}

// A utility function to print an array


__device__ unsigned int device_data[SIZE];

// naive warp-level bitwise radix sort
__global__ void radixSortParralel(unsigned int* sort_tmp, unsigned int* sort_tmp_1) {
	const int num_lists = SIZE;
	const int num_element = 10;

	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int stride = blockDim.x * gridDim.x;   
	//printf("\n%d num list: %d" ,tid, num_lists );
	//printf("\n num element: %d" , num_element);
	if (sort_tmp == NULL) {
		printf("\nsort_tmp is null");
	}
	

    for (int bit = 0; bit < BIT; bit++) {
		const int bit_mask = (1 << bit);
		int base_cnt_0 = 0, base_cnt_1 = 0;
		//printf("\ntest1 + %d",tid );
		for (int i = 0; i < num_element; i += num_lists) {
			int elem = sort_tmp[i + tid];
			
			//printf("\ntest %d + %d bit %d",tid, elem, bit_mask);
			if ((elem & bit_mask) > 0) {
				sort_tmp_1[base_cnt_1 + tid] = elem;
				base_cnt_1 += num_lists;
				//printf("\n element 1 %d = %d",i, elem);
			}
			else {
				sort_tmp[base_cnt_0 + tid] = elem;
				base_cnt_0 += num_lists;
				//printf("\n element 0 %d = %d",i, elem);
			}
		}
		//copy data back to source - first the zero list
		for (int i = 0; i < base_cnt_1; i += num_lists) {
			sort_tmp[base_cnt_0 + i + tid] = sort_tmp_1[i + tid];
			//printf("\nsort tmp: %d", sort_tmp_0[i + tid]);
		}
	}
	//cudaDeviceSynchronize();
	//__syncthreads();

	//for (int i = 0; i < num_element; i++) {
	//	printf("\n%d", sort_tmp[i]);
	//}
}

__global__ void radixSort(unsigned int* data, unsigned int* data_tmp) {
	//__shared__ volatile unsigned int shared_data[SIZE * 2];
	// load from global into shared variable
	int tid = threadIdx.x + blockIdx.x * blockDim.x;

	while (tid < SIZE) {
		//data[tid] = data[tid] + 1000;



		tid += blockDim.x * gridDim.x;
	}

}
void print(unsigned int* arr, int n)
{
	for (int i = 0; i < n; i++)
		cout << (double)arr[i] / 100 << endl;
}

void menu() {
	cout << "Izberi: " << endl;
	cout << "0 Exit " << endl;
	cout << " 1 Serial sort" << endl;
	cout << " 2 Parrallel sort" << endl;
}

vector<unsigned int> ReadFile(string filepath) {
	vector<unsigned int> vector;
	ifstream inputFile(filepath);
	double value;

	//read the elements in the file into a vector  
	while (inputFile >> value) {
		//bitset<32>tmp(value *100);
		//cout << tmp << " ";
		vector.push_back(value * 100);
		//cout << (int)(tmp.to_ulong()) << " "; // konvert nazaj to int
	}
	/*for each (unsigned int var in vector)
	{
	cout << var << " " << endl;
	}*/
	inputFile.close();
	return vector;
}

int main(int argc, char * argv[]) {

	int blockSize = 32;
	int numBlocks = (SIZE + blockSize - 1) / blockSize;

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	float totalTime = 0;

	std::ifstream fin(argv[1]);
	std::vector<unsigned int> host_double(SIZE);

	unsigned int host_data[SIZE];

	if (!fin) {
		cerr << "Datoteka ne obstaja " << argv[1] << endl;
	}
	else {
		cout << "Datoteka obstaja " << argv[1] << endl;
		host_double = ReadFile(argv[1]);
		
		//WriteFile(host_double);
	}
	int size = host_double.size();

	int izbiraAlg;
	menu();
	cin >> izbiraAlg;
	if (izbiraAlg == 0) {

	}else if (izbiraAlg == 1) {

		//Serial radix sort

		unsigned int* da = host_double.data();
		//int arr[] = { 170, 45, 75, 90, 802, 24, 2, 66, 170, 45, 75, 90, 802, 24, 2, 66 };
		int n = sizeof(host_double) / sizeof(host_double[0]);
		//auto t1 = high_resolution_clock::now();
		high_resolution_clock::time_point t1 = high_resolution_clock::now();
		serialRadixSort(da, size);
		high_resolution_clock::time_point t2 = high_resolution_clock::now();
		auto duration = duration_cast <milliseconds> (t2 - t1).count();
		//auto t2 = high_resolution_clock::now();
		//auto diff = duration_cast<duration<double>>(t2 - t1);
		// now elapsed time, in seconds, as a double can be found in diff.count()
		//long ms = (long)(1000 * diff.count());
		cout << endl << "Serial radix sort " << endl;
		//print(da, size);
		for (int i = 0; i < size - 1; i++) {
			if (host_double[i] > host_double[i + 1]) {
				printf("sort error at, hdata[%d] = %d, hdata[%d] = %d\n", i, host_double[i], i + 1, host_double[i + 1]);
				return 1;
			}
		}
		cout << "Total time: " << duration << "ms" << endl;
		cout << "Success" << endl;

	}
	else {

		//Parrallel radix sort
		cout << endl;
		//unsigned int t = decToBinary(10);
		//unsigned int tmp22;
		/*bitset<32>tmp(65.55);
		cout << tmp << endl;*/
			
		cout << "Parrallel sort " << endl;
		std::copy(host_double.begin(), host_double.end(), host_data);
		cout << "Not sorted!" << endl;
		print(host_data, size);

		unsigned int* ddata;
		unsigned int* ddata_tmp;
		//

		//for (int lcount = 0; lcount < LOOPS; lcount++) {
		//	unsigned int range = 1U << UPPER_BIT;
		//	//for (int i = 0; i < SIZE; i++) host_data[i] = rand() % range;
			//cudaMemcpyToSymbol(device_data, host_data, SIZE * sizeof(unsigned int));
			cudaMalloc((void**)&ddata_tmp, SIZE * sizeof(unsigned int));
			cudaMalloc((void**)&ddata, SIZE * sizeof(unsigned int));
			cudaMemcpy(ddata, host_data, SIZE * sizeof(unsigned int), cudaMemcpyHostToDevice);			

			high_resolution_clock::time_point start = high_resolution_clock::now();
			radixSortParralel << <1, 10 >> >(ddata, ddata_tmp);
			//radixSort << <(numBlocks, blockSize) >> > (ddata, ddata_tmp);
			high_resolution_clock::time_point stop = high_resolution_clock::now();

			cudaMemcpy(host_data, ddata, SIZE * sizeof(unsigned int), cudaMemcpyDeviceToHost);

			//cudaMemcpyFromSymbol(host_data, device_data, SIZE * sizeof(unsigned int));
			auto durationParralel = duration_cast <milliseconds> (start - stop).count();
		//	/*for (int i = 0; i < SIZE - 1; i++) {
		//	if (host_data[i] > host_data[i + 1]) {
		//	printf("sort error at loop %d, hdata[%d] = %d, hdata[%d] = %d\n", lcount, i, host_data[i], i + 1, host_data[i + 1]);
		//	return 1;
		//	}
		//	}*/
		//	cout << " Sorted data: " << endl;
		//	for (int i = 0; i < size; i++) {
		//		cout << i << "     " <<(double)host_data[i] / 100 << endl;
		//	}
		//	cout << " Time parralel sort: " << durationParralel << endl;
		//}
			cout << endl << "Sorted" << endl;
			print(host_data, size);
		
	}
	
	
	return 0;
}
