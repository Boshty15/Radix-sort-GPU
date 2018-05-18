
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include "kernel.h"
#include <iostream>     
#include <fstream>
#include <vector>

#define SIZE 512

using namespace std;

__device__ unsigned int device_data[SIZE];

//Test v gpu in nazaj
__global__ void gpucopy(int* src, int* dst)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;;
	dst[i] = src[i];
}

__global__ void radix_sort(int* const sort_tmp, const int num_list, const int num_elements, const int tid, int* const sort_tmp_1) {
	//sort into num_list, lists
	//aplly radixx sort on 32 bits of data
	for (int bit = 0; bit < 32; bit++) {
		const unsigned int bit_mask = (1 << bit);
		unsigned int base_cnt_0 = 0, base_cnt_1 = 0;

		for (int i = 0; i < num_elements; i++) {
			const int element = sort_tmp[i + tid];

			if ((element & bit_mask) > 0) {
				sort_tmp_1[base_cnt_1 + tid] = element;
				base_cnt_1 += num_list;

			}
			else {
				sort_tmp[base_cnt_1 + tid] = element;
				base_cnt_0 += num_elements;
			}
		}
		//copy data back to source from the ones's list
		for (int i = 0; i < base_cnt_1; i++) {
			sort_tmp[base_cnt_0 + i + tid] = sort_tmp_1[i + tid];
		}


	}
}

//__global__ void	ParalelRadixSort()
//{
//	//sprememba
//	__shared__ volatile unsigned int shared_data[SIZE * 2];
//
//	//spremena za commit
//	shared_data[threadIdx.x] = device_data[threadIdx.x];
//		
//	unsigned int bitmask = 1 << 0;
//	unsigned int offset = 0;
//	// -1, -2, -4, -8, -16, -32, -64, -128, -256,...
//	unsigned int thrmask = 0xFFFFFFFFU << threadIdx.x;
//	unsigned int mypos;
//
//	for (int i = 0; i <= 10; i++)
//	{
//		unsigned int mydata = shared_data[((SIZE - 1) - blockIdx.x * blockDim.x + threadIdx.x) + offset];
//		unsigned int mybit = mydata&bitmask;
//		// Get population of ones and zeroes
//		unsigned int ones = __ballot(mybit);
//		unsigned int zeroes = ~ones;
//		// Switch ping-pong buffers
//		offset ^= SIZE;
//
//		// Do zeroes, then ones
//		if (!mybit)
//		{
//			mypos = __popc(zeroes&thrmask);
//		}
//		else {      // Threads with a one bit
//					// Get my position in ping-pong buffer
//			mypos = __popc(zeroes) + __popc(ones&thrmask);
//		}
//
//		// Move to buffer  (or use shfl for cc 3.0)
//		shared_data[mypos - 1 + offset] = mydata;
//		// Repeat for next bit
//		bitmask <<= 1;
//
//		device_data[threadIdx.x] = shared_data[threadIdx.x + offset];
//	}
//		
//
//}

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
		vector.push_back(value * 10);
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
	
	std::vector<int> host_double(9);
	std::vector<int> host_double2(9);

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
	int* device2 = NULL;

	cudaMalloc((void**)&device, size * sizeof(int));
	cudaMalloc((void**)&device2, size * sizeof(int));

	//num_list 32 

	cudaMemcpy(device, host_double.data(), size * sizeof(int), cudaMemcpyHostToDevice);
	
	radix_sort << <(1, size) >> > ();
	//ParalelRadixSort << < 1, size >> >();
	//gpucopy << <(size + size - 1) / size, size >> >(device, device2);
	
	cudaMemcpy(host_double2.data(), device2, size * sizeof(int), cudaMemcpyDeviceToHost);

	cudaDeviceSynchronize();

	cout << endl << "host data" << endl;
	for each (int var in host_double2)
	{
		cout << var << " ";
	}

	cudaFree(device);

	//int size = vector.size();
	////vector.resize(size);

	//std::copy(vector.begin(), vector.end(), hdata);
	////int *a = vector.data();

	//// Copy data from host to device
	//cudaMemcpyToSymbol(device_data, hdata, SIZE * sizeof(unsigned int));

	//// Execution time measurement, that point starts the clock
	//cudaEventRecord(start);
	//ParalelRadixSort << < 1, SIZE >> >();
	//	
	//// Execution time measurement, that point stops the clock
	//cudaEventRecord(stop);

	//// Make kernel function synchronous
	//cudaDeviceSynchronize();

	//cudaEventElapsedTime(&milliseconds, start, stop);
	//totalTime += milliseconds;
	////milliseconds = 0;

	//// Copy data from device to host
	//cudaMemcpyFromSymbol(hdata, device_data, SIZE * sizeof(unsigned int));

	//std::vector<int> v(hdata, hdata + sizeof hdata / sizeof hdata[0]);
	//v.resize(size);

	///*printf("TIME %4.2fs", milliseconds);
	//printf("\n");
	//printf("Effective Bandwidth (GB/s): %fn", SIZE * 4 * 3 / milliseconds / 1e6);
	//printf("\n");*/
	//cout << endl;
	//cout << "Host Data" << endl;
	//for each (int var in v)
	//{
	//	cout << var << " ";
	//}
	///*cout << "Device Data" << endl;
	//for each (int var in device_data)
	//{
	//	cout << var << " ";
	//}
	//*/
	//cout << endl;
	//cout << "TIME %4.2fs " << milliseconds << endl;

    return 0;
}







//
//#define nTPB 512
//
//__global__ void gpucopy(int* src, int* dst)
//{
//	int i = blockIdx.x * blockDim.x + threadIdx.x;;
//	dst[i] = src[i];
//}
//
//int main()
//{
//	const int arraySize = 500; // >= 1025 will fail on my system!
//
//	int* data1 = new int[arraySize];
//	int* data2 = new int[arraySize];
//	// Initialized both data1 and data2
//	// ... 
//	for (int i = 0; i < arraySize; i++) {
//		data1[i] = 2 * i;
//		//cout << data1[i] + " ";
//	}
//		
//
//	int* dev_data1 = NULL;
//	int* dev_data2 = NULL;
//	// Initialized both dev_data1 and dev_data2
//	// ... 
//	cudaMalloc(&dev_data1, arraySize * sizeof(int));
//	cudaMalloc(&dev_data2, arraySize * sizeof(int));
//
//	// copy data1 to device
//	cudaMemcpy(dev_data1, data1, arraySize * sizeof(int), cudaMemcpyHostToDevice);
//
//	//// copy dev_data1 to dev_data2 with gpu
//	//gpucopy << <1, arraySize >> >(dev_data1, dev_data2);
//	gpucopy << <(arraySize + nTPB - 1) / nTPB, nTPB >> >(dev_data1, dev_data2);
//
//	// copy dev_data2 to data
//	cudaMemcpy(data2, dev_data2, arraySize * sizeof(int), cudaMemcpyDeviceToHost);
//
//
//	for (int i = 0; i<arraySize; i++)
//		if (data2[i] != data1[i])
//			cout << "Error: data is different - data2[" << i << "] is " << data2[i] << endl;
//
//	return 0;
//}
//
//
//
