
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include "kernel.h"
#include <iostream>     
#include <fstream>
#include <vector>

#define SIZE 32

using namespace std;

__device__ unsigned int device_data[SIZE];

__global__ void	ParalelRadixSort()
{
	//sprememba
	__shared__ volatile unsigned int shared_data[SIZE * 2];

	shared_data[threadIdx.x] = device_data[threadIdx.x];
		
	unsigned int bitmask = 1 << 0;
	unsigned int offset = 0;
	// -1, -2, -4, -8, -16, -32, -64, -128, -256,...
	unsigned int thrmask = 0xFFFFFFFFU << threadIdx.x;
	unsigned int mypos;

	for (int i = 0; i <= 10; i++)
	{
		unsigned int mydata = shared_data[((SIZE - 1) - threadIdx.x) + offset];
		unsigned int mybit = mydata&bitmask;
		// Get population of ones and zeroes
		unsigned int ones = __ballot(mybit);
		unsigned int zeroes = ~ones;
		// Switch ping-pong buffers
		offset ^= SIZE;

		// Do zeroes, then ones
		if (!mybit)
		{
			mypos = __popc(zeroes&thrmask);
		}
		else {      // Threads with a one bit
					// Get my position in ping-pong buffer
			mypos = __popc(zeroes) + __popc(ones&thrmask);
		}

		// Move to buffer  (or use shfl for cc 3.0)
		shared_data[mypos - 1 + offset] = mydata;
		// Repeat for next bit
		bitmask <<= 1;

		device_data[threadIdx.x] = shared_data[threadIdx.x + offset];
	}
		

}
void WriteFile(vector<int> vector) {
	ofstream myfile;
	myfile.open("out.txt");
	for each (int var in vector)
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
		vector.push_back(value);
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
	

	vector<int> vector;
	

	unsigned int hdata[SIZE];
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
		vector = ReadFile(argv[1]);

		WriteFile(vector);
	}

	int size = vector.size();
	//vector.resize(size);

	std::copy(vector.begin(), vector.end(), hdata);
	//int *a = vector.data();

	// Copy data from host to device
	cudaMemcpyToSymbol(device_data, hdata, SIZE * sizeof(unsigned int));

	// Execution time measurement, that point starts the clock
	cudaEventRecord(start);
	ParalelRadixSort << < 1, SIZE >> >();
		
	// Execution time measurement, that point stops the clock
	cudaEventRecord(stop);

	// Make kernel function synchronous
	cudaDeviceSynchronize();

	cudaEventElapsedTime(&milliseconds, start, stop);
	totalTime += milliseconds;
	//milliseconds = 0;

	// Copy data from device to host
	cudaMemcpyFromSymbol(hdata, device_data, SIZE * sizeof(unsigned int));

	std::vector<int> v(hdata, hdata + sizeof hdata / sizeof hdata[0]);
	v.resize(size);

	/*printf("TIME %4.2fs", milliseconds);
	printf("\n");
	printf("Effective Bandwidth (GB/s): %fn", SIZE * 4 * 3 / milliseconds / 1e6);
	printf("\n");*/
	cout << endl;
	cout << "Host Data" << endl;
	for each (int var in v)
	{
		cout << var << " ";
	}
	/*cout << "Device Data" << endl;
	for each (int var in device_data)
	{
		cout << var << " ";
	}
	*/
	cout << endl;
	cout << "TIME %4.2fs " << milliseconds << endl;

    return 0;
}



