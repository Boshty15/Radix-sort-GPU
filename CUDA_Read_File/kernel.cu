//
//
//
//#pragma once
//#ifdef __INTELLISENSE__
//void __syncthreads();
//#endif


#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cuda_runtime_api.h>
#include <stdio.h>
#include "kernel.h"
#include <iostream>     
#include <fstream>
#include <vector>
#include <bitset>
#include <time.h>
#include <ctime>
#include <chrono>
#include <math.h> 

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/generate.h>
#include <thrust/sort.h>
#include <thrust/copy.h>
#include <thrust/execution_policy.h>
#include <algorithm>
#include <cstdlib>
#include <Windows.h>

using namespace std;
using namespace std::chrono;

size_t getMax2(vector<string> arr, int n) {
	size_t max = arr[0].size();
	for (int i = 1; i < n; i++) {
		if (arr[i].size()>max)
			max = arr[i].size();
	}
	return max;
}

void countSort2(vector<string> a, int size, size_t k) {
	string *b = NULL; int *c = NULL;
	b = new string[size];
	c = new int[257];



	for (int i = 0; i <257; i++) {
		c[i] = 0;
		//cout << c[i] << "\n";
	}
	for (int j = 0; j <size; j++) {
		c[k < a[j].size() ? (int)(unsigned char)a[j][k] + 1 : 0]++;            //a[j] is a string
																			   //cout << c[a[j]] << endl;
	}

	for (int f = 1; f <257; f++) {
		c[f] += c[f - 1];
	}

	for (int r = size - 1; r >= 0; r--) {
		b[c[k < a[r].size() ? (int)(unsigned char)a[r][k] + 1 : 0] - 1] = a[r];
		c[k < a[r].size() ? (int)(unsigned char)a[r][k] + 1 : 0]--;
	}

	for (int l = 0; l < size; l++) {
		a[l] = b[l];
	}

	// avold memory leak
	delete[] b;
	delete[] c;
}


void radixSort2(vector<string> b, int r) {
	size_t max = getMax2(b, r);
	for (size_t digit = max; digit > 0; digit--) { // size_t is unsigned, so avoid using digit >= 0, which is always true
		countSort2(b, r, digit - 1);
	}

}



template <class T>
T getMax(T* arr, int n)
{
	T mx = arr[0];
	for (int i = 1; i < n; i++)
		if (arr[i] > mx)
			mx = arr[i];
	return mx;
}

template <class T>
void countSort(T* arr, int n, int exp)
{
	// output array
	vector<T> ve(n);
	T* output = ve.data();
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

// Radix Sort
template <class T>
void serialRadixSort(T* arr, int n)
{
	// Find the maximum number to know number of digits
	T m = getMax(arr, n);

	// Do counting sort for every digit. Note that instead
	// of passing digit number, exp is passed. exp is 10^i
	// where i is current digit number
	for (int exp = 1; m / exp > 0; exp *= 10)
		countSort(arr, n, exp);
}


int main(int argc, char * argv[]) {

	//host
	//thrust::sort(thrust::host, A, A + N);

	unsigned long long int N;
	cout << "Izberi število elementov:" << endl << endl;
	cout << "(1) 1M elementov:  " << endl;
	cout << "(2) 5M elementov:  " << endl;
	cout << "(3) 10M elementov:  " << endl;
	cout << "(4) 25M elementov:  " << endl;
	cout << "(5) 50M elementov:  " << endl;
	cout << "(6) 75M elementov:  " << endl;
	cout << "(7) 100M elementov:  " << endl;
	cin >> N;

	switch (N)
	{
	case 1:
		N = 1000000;
		break;
	case 2:
		N = 5000000;
		break;
	case 3:
		N = 10000000;
		break;
	case 4:
		N = 25000000;
		break;
	case 5:
		N = 50000000;
		break;
	case 6:
		N = 75000000;
		break;
	case 7:
		N = 100000000;
		break;
	default:
		break;
	}
	cout << "Izbral si: " << N << " elementov" << endl;


	cout << "(1) Integer" << endl;
	cout << "(2) Char" << endl;
	cout << "(3) Long" << endl;
	cout << "(4) Float" << endl;
	cout << "(5) Double" << endl;
	int tip;
	cin >> tip;

	switch (tip)
	{
	case 1: 
	{
		cout << endl;
		cout << endl;
		cout << "Int" << endl << endl;

		

		//cout << "Parralel: HOST Stable Radix sort" << endl;

		//thrust::host_vector<unsigned int> host_int_Parralel_Host(N);
		//std::generate(host_int_Parralel_Host.begin(), host_int_Parralel_Host.end(), rand);

		//cout << "Parralel size: " << host_int_Parralel_Host.size() << endl;

		//thrust::device_vector<unsigned int> d_vec_Host = host_int_Parralel_Host;

		//high_resolution_clock::time_point start_Host = high_resolution_clock::now();
		//// sort data on the device (846M keys per second on GeForce GTX 480)
		//thrust::stable_sort(thrust::host,d_vec_Host.begin(), d_vec_Host.end());
		//high_resolution_clock::time_point stop_Host = high_resolution_clock::now();

		//// transfer data back to host
		//thrust::copy(d_vec_Host.begin(), d_vec_Host.end(), host_int_Parralel_Host.begin());

		//auto durationParralel_Host = duration_cast <milliseconds> (stop_Host - start_Host).count();
		//cout << "Time: " << durationParralel_Host << " milliseconds" << endl;

		//cout << endl;
		//cout << endl;


		cout << "Parralel:  Stable Radix sort" << endl;

		thrust::host_vector<unsigned int> host_int_Parralel_S(N);
		std::generate(host_int_Parralel_S.begin(), host_int_Parralel_S.end(), rand);

		cout << "Parralel size: " << host_int_Parralel_S.size() << endl;

		thrust::device_vector<unsigned int> d_vec_S = host_int_Parralel_S;

		high_resolution_clock::time_point start_S = high_resolution_clock::now();
		// sort data on the device (846M keys per second on GeForce GTX 480)
		thrust::stable_sort(d_vec_S.begin(), d_vec_S.end());
		high_resolution_clock::time_point stop_S = high_resolution_clock::now();

		// transfer data back to host
		thrust::copy(d_vec_S.begin(), d_vec_S.end(), host_int_Parralel_S.begin());

		auto durationParralel_S = duration_cast <milliseconds> (stop_S - start_S).count();
		cout << "Time: " << durationParralel_S << " milliseconds" << endl;

		cout << endl;
		cout << endl;

		//////thrust parralel
		cout << "Parralel: Radix sort" << endl;

		thrust::host_vector<unsigned int> host_int_Parralel(N);
		std::generate(host_int_Parralel.begin(), host_int_Parralel.end(), rand);

		cout << "Parralel size: " << host_int_Parralel.size() << endl;

		thrust::device_vector<unsigned int> d_vec = host_int_Parralel;

		high_resolution_clock::time_point start = high_resolution_clock::now();
		// sort data on the device (846M keys per second on GeForce GTX 480)
		thrust::sort(d_vec.begin(), d_vec.end());
		high_resolution_clock::time_point stop = high_resolution_clock::now();

		// transfer data back to host
		thrust::copy(d_vec.begin(), d_vec.end(), host_int_Parralel.begin());

		auto durationParralel = duration_cast <milliseconds> (stop - start).count();
		cout << "Time: " << durationParralel << " milliseconds" << endl;

		cout << endl;
		cout << endl;

		cout << "Serial Radix sort" << endl;

		thrust::host_vector<unsigned int> host_int_Serial(N);
		std::generate(host_int_Serial.begin(), host_int_Serial.end(), rand);

		cout << "Serial size: " << host_int_Serial.size() << endl;

		//thrust::device_vector<unsigned long int> d_vec_S;
		unsigned int* da = host_int_Serial.data();
		high_resolution_clock::time_point startS = high_resolution_clock::now();
		// sort data on the device (846M keys per second on GeForce GTX 480)
		serialRadixSort<unsigned int>(da, N);
		high_resolution_clock::time_point stopS = high_resolution_clock::now();

		auto durationS = duration_cast <milliseconds> (stopS - startS).count();
		cout << "Time: " << durationS << " milliseconds" << endl;

		cout << endl;
		cout << endl; 
		break;
	}
	case 2:
	{
		cout << endl;
		cout << endl;
		cout << "Char" << endl << endl;

		cout << "Parralel:  Stable Radix sort" << endl;

		thrust::host_vector<char> host_int_Parralel_S(N);
		std::generate(host_int_Parralel_S.begin(), host_int_Parralel_S.end(), rand);

		cout << "Parralel size: " << host_int_Parralel_S.size() << endl;

		thrust::device_vector<char> d_vec_S = host_int_Parralel_S;

		high_resolution_clock::time_point start_S = high_resolution_clock::now();
		// sort data on the device (846M keys per second on GeForce GTX 480)
		thrust::stable_sort(d_vec_S.begin(), d_vec_S.end());
		high_resolution_clock::time_point stop_S = high_resolution_clock::now();

		// transfer data back to host
		thrust::copy(d_vec_S.begin(), d_vec_S.end(), host_int_Parralel_S.begin());

		auto durationParralel_S = duration_cast <milliseconds> (stop_S - start_S).count();
		cout << "Time: " << durationParralel_S << " milliseconds" << endl;

		cout << endl;
		cout << endl;

		//////thrust parralel
		cout << "Parralel: Radix sort" << endl;

		thrust::host_vector<char> host_int_Parralel(N);
		std::generate(host_int_Parralel.begin(), host_int_Parralel.end(), rand);

		cout << "Parralel size: " << host_int_Parralel.size() << endl;
		//cout << host_int_Parralel.data() << endl;

		thrust::device_vector<char> d_vec = host_int_Parralel;

		high_resolution_clock::time_point start = high_resolution_clock::now();
		// sort data on the device (846M keys per second on GeForce GTX 480)
		thrust::sort(d_vec.begin(), d_vec.end());
		high_resolution_clock::time_point stop = high_resolution_clock::now();

		// transfer data back to host
		thrust::copy(d_vec.begin(), d_vec.end(), host_int_Parralel.begin());
		//cout << host_int_Parralel.data() << endl;

		auto durationParralel = duration_cast <milliseconds> (stop - start).count();
		cout << "Time: " << durationParralel << " milliseconds" << endl;

		cout << endl;
		cout << endl;
	
		
		string t(host_int_Parralel.begin(), host_int_Parralel.end());
		vector<string> tmp;
		tmp.push_back(t);

		

		/*countingRadixSort(tmp, N);*/

		//cout << "Serial Radix sort" << endl;

		//thrust::host_vector<char> host_int_Serial(N);
		//std::generate(host_int_Serial.begin(), host_int_Serial.end(), rand);

		//cout << "Serial size: " << host_int_Serial.size() << endl;

		////thrust::device_vector<unsigned long int> d_vec_S;
		//char* da = host_int_Serial.data();
		high_resolution_clock::time_point startS = high_resolution_clock::now();
		//// sort data on the device (846M keys per second on GeForce GTX 480)
		//serialRadixSort<char>(da, N);
		radixSort2(tmp, (int)(sizeof(tmp) / sizeof(tmp[0])));

		high_resolution_clock::time_point stopS = high_resolution_clock::now();

		auto durationS = duration_cast <milliseconds> (stopS - startS).count();
		cout << "Time: " << durationS << " milliseconds" << endl;

		//cout << endl;
		//cout << endl;
		break;
	}
	case 3:
	{
		cout << endl;
		cout << endl;
		cout << "Long" << endl << endl;

		cout << "Parralel:  Stable Radix sort" << endl;

		thrust::host_vector<long> host_int_Parralel_S(N);
		std::generate(host_int_Parralel_S.begin(), host_int_Parralel_S.end(), rand);

		cout << "Parralel size: " << host_int_Parralel_S.size() << endl;

		thrust::device_vector<long> d_vec_S = host_int_Parralel_S;

		high_resolution_clock::time_point start_S = high_resolution_clock::now();
		// sort data on the device (846M keys per second on GeForce GTX 480)
		thrust::stable_sort(d_vec_S.begin(), d_vec_S.end());
		high_resolution_clock::time_point stop_S = high_resolution_clock::now();

		// transfer data back to host
		thrust::copy(d_vec_S.begin(), d_vec_S.end(), host_int_Parralel_S.begin());

		auto durationParralel_S = duration_cast <milliseconds> (stop_S - start_S).count();
		cout << "Time: " << durationParralel_S << " milliseconds" << endl;

		cout << endl;
		cout << endl;

		//////thrust parralel
		cout << "Parralel: Radix sort" << endl;

		thrust::host_vector<long> host_int_Parralel(N);
		std::generate(host_int_Parralel.begin(), host_int_Parralel.end(), rand);

		cout << "Parralel size: " << host_int_Parralel.size() << endl;

		thrust::device_vector<long> d_vec = host_int_Parralel;

		high_resolution_clock::time_point start = high_resolution_clock::now();
		// sort data on the device (846M keys per second on GeForce GTX 480)
		thrust::sort(d_vec.begin(), d_vec.end());
		high_resolution_clock::time_point stop = high_resolution_clock::now();

		// transfer data back to host
		thrust::copy(d_vec.begin(), d_vec.end(), host_int_Parralel.begin());

		auto durationParralel = duration_cast <milliseconds> (stop - start).count();
		cout << "Time: " << durationParralel << " milliseconds" << endl;

		cout << endl;
		cout << endl;

		cout << "Serial Radix sort" << endl;

		thrust::host_vector<long> host_int_Serial(N);
		std::generate(host_int_Serial.begin(), host_int_Serial.end(), rand);

		cout << "Serial size: " << host_int_Serial.size() << endl;

		//thrust::device_vector<unsigned long int> d_vec_S;
		long* da = host_int_Serial.data();
		high_resolution_clock::time_point startS = high_resolution_clock::now();
		// sort data on the device (846M keys per second on GeForce GTX 480)
		serialRadixSort<long>(da, N);
		high_resolution_clock::time_point stopS = high_resolution_clock::now();

		auto durationS = duration_cast <milliseconds> (stopS - startS).count();
		cout << "Time: " << durationS << " milliseconds" << endl;

		cout << endl;
		cout << endl;
		break;
	}
	case 4:
	{
		cout << endl;
		cout << endl;
		cout << "float" << endl << endl;

		cout << "Parralel:  Stable Radix sort" << endl;


		thrust::host_vector<float> host_int_Parralel_S(N);
		std::generate(host_int_Parralel_S.begin(), host_int_Parralel_S.end(), rand);

		cout << "Parralel size: " << host_int_Parralel_S.size() << endl;

		thrust::device_vector<float> d_vec_S = host_int_Parralel_S;

		high_resolution_clock::time_point start_S = high_resolution_clock::now();
		// sort data on the device (846M keys per second on GeForce GTX 480)
		thrust::stable_sort(d_vec_S.begin(), d_vec_S.end());
		high_resolution_clock::time_point stop_S = high_resolution_clock::now();

		// transfer data back to host
		thrust::copy(d_vec_S.begin(), d_vec_S.end(), host_int_Parralel_S.begin());

		auto durationParralel_S = duration_cast <milliseconds> (stop_S - start_S).count();
		cout << "Time: " << durationParralel_S << " milliseconds" << endl;

		cout << endl;
		cout << endl;

		//////thrust parralel
		cout << "Parralel: Radix sort" << endl;

		thrust::host_vector<float> host_int_Parralel(N);
		std::generate(host_int_Parralel.begin(), host_int_Parralel.end(), rand);

		cout << "Parralel size: " << host_int_Parralel.size() << endl;

		thrust::device_vector<float> d_vec = host_int_Parralel;

		high_resolution_clock::time_point start = high_resolution_clock::now();
		// sort data on the device (846M keys per second on GeForce GTX 480)
		thrust::sort(d_vec.begin(), d_vec.end());
		high_resolution_clock::time_point stop = high_resolution_clock::now();

		// transfer data back to host
		thrust::copy(d_vec.begin(), d_vec.end(), host_int_Parralel.begin());

		auto durationParralel = duration_cast <milliseconds> (stop - start).count();
		cout << "Time: " << durationParralel << " milliseconds" << endl;

		cout << endl;
		cout << endl;

		cout << "Serial Radix sort" << endl;

		thrust::host_vector<float> host_int_Serial(N);
		thrust::host_vector<unsigned int> host_int_Serial_Int(0);
		std::generate(host_int_Serial.begin(), host_int_Serial.end(), rand);
		for each (float var in host_int_Serial)
		{
			host_int_Serial_Int.push_back(var);
			//cout << var << endl;
		}
		cout << host_int_Serial_Int.size() << endl;

		cout << "Serial size: " << host_int_Serial.size() << endl;

		//thrust::device_vector<unsigned long int> d_vec_S;
		unsigned int* da = host_int_Serial_Int.data();
		high_resolution_clock::time_point startS = high_resolution_clock::now();
		// sort data on the device (846M keys per second on GeForce GTX 480)
		serialRadixSort<unsigned int>(da, N);
		high_resolution_clock::time_point stopS = high_resolution_clock::now();

		auto durationS = duration_cast <milliseconds> (stopS - startS).count();
		cout << "Time: " << durationS << " milliseconds" << endl;

		cout << endl;
		cout << endl;
		break;
	}
	case 5:
	{
		cout << endl;
		cout << endl;
		cout << "double" << endl << endl;

		cout << "Parralel:  Stable Radix sort" << endl;

		thrust::host_vector<double> host_int_Parralel_S(N);
		std::generate(host_int_Parralel_S.begin(), host_int_Parralel_S.end(), rand);

		cout << "Parralel size: " << host_int_Parralel_S.size() << endl;

		thrust::device_vector<double> d_vec_S = host_int_Parralel_S;

		high_resolution_clock::time_point start_S = high_resolution_clock::now();
		// sort data on the device (846M keys per second on GeForce GTX 480)
		thrust::stable_sort(d_vec_S.begin(), d_vec_S.end());
		high_resolution_clock::time_point stop_S = high_resolution_clock::now();

		// transfer data back to host
		thrust::copy(d_vec_S.begin(), d_vec_S.end(), host_int_Parralel_S.begin());

		auto durationParralel_S = duration_cast <milliseconds> (stop_S - start_S).count();
		cout << "Time: " << durationParralel_S << " milliseconds" << endl;

		cout << endl;
		cout << endl;

		//////thrust parralel
		cout << "Parralel: Radix sort" << endl;

		thrust::host_vector<double> host_int_Parralel(N);
		std::generate(host_int_Parralel.begin(), host_int_Parralel.end(), rand);

		cout << "Parralel size: " << host_int_Parralel.size() << endl;

		thrust::device_vector<double> d_vec = host_int_Parralel;

		high_resolution_clock::time_point start = high_resolution_clock::now();
		// sort data on the device (846M keys per second on GeForce GTX 480)
		thrust::sort(d_vec.begin(), d_vec.end());
		high_resolution_clock::time_point stop = high_resolution_clock::now();

		// transfer data back to host
		thrust::copy(d_vec.begin(), d_vec.end(), host_int_Parralel.begin());

		auto durationParralel = duration_cast <milliseconds> (stop - start).count();
		cout << "Time: " << durationParralel << " milliseconds" << endl;

		cout << endl;
		cout << endl;

		cout << "Serial Radix sort" << endl;

		thrust::host_vector<double> host_int_Serial(N);
		thrust::host_vector<unsigned int> host_int_Serial_IntD(0);
		std::generate(host_int_Serial.begin(), host_int_Serial.end(), rand);

		cout << "Serial size: " << host_int_Serial.size() << endl;
		for each (double var in host_int_Serial)
		{
			host_int_Serial_IntD.push_back(var);
			//cout << var << endl;
		}

		//thrust::device_vector<unsigned long int> d_vec_S;
		unsigned int* da = host_int_Serial_IntD.data();
		high_resolution_clock::time_point startS = high_resolution_clock::now();
		// sort data on the device (846M keys per second on GeForce GTX 480)
		serialRadixSort<unsigned int>(da, N);
		high_resolution_clock::time_point stopS = high_resolution_clock::now();

		auto durationS = duration_cast <milliseconds> (stopS - startS).count();
		cout << "Time: " << durationS << " milliseconds" << endl;

		cout << endl;
		cout << endl;
		break;
	}
	default:
		
		break;
	}

	
	
	 //sort data on the device (846M keys per second on GeForce GTX 480)
	//thrust::sort(thrust::host, d_vec_S.begin(), d_vec_S.end());
	
	
	 //transfer data back to host
	//thrust::copy(d_vec_S.begin(), d_vec_S.end(), host_int_Serial.begin());
//	auto durationSerial = duration_cast <milliseconds> (stopS - startS).count();
	/*for (int i = 0; i < N; i++) {
		if (da[i] > da[i + 1]) {
			printf("sort error at, hdata[%d] = %d, hdata[%d] = %d\n", i, da[i], i + 1, da[i + 1]);
			return 1;
		}
	}*/



	

	 
	

}





//
//#define BIT 32
//#define SIZE (33 * 1024)
//
///*Serial radix sort*/



//__device__ void radixSortParralel(unsigned int* sort_tmp, const int num_lists, const int num_element, const int tid, unsigned int* sort_tmp_1) {
//	//if (sort_tmp == NULL) {
//	//	printf("\nsort_tmp is null");
//	//}
//
//	for (int bit = 0; bit < BIT; bit++) {
//		const int bit_mask = (1 << bit);
//		int base_cnt_0 = 0, base_cnt_1 = 0;
//		//printf("\ntest1 + %d",tid );
//		for (int i = 0; i < num_element; i += num_lists) {
//			int elem = sort_tmp[i + tid];
//
//			//printf("\ntest %d + %d bit %d",tid, elem, bit_mask);
//			if ((elem & bit_mask) > 0) {
//				sort_tmp_1[base_cnt_1 + tid] = elem;
//				base_cnt_1 += num_lists;
//				//printf("\n element 1 %d = %d",i, elem);
//			}
//			else {
//				sort_tmp[base_cnt_0 + tid] = elem;
//				base_cnt_0 += num_lists;
//				//printf("\n element 0 %d = %d",i, elem);
//			}
//		}
//		//copy data back to source - first the zero list
//		for (int i = 0; i < base_cnt_1; i += num_lists) {
//			sort_tmp[base_cnt_0 + i + tid] = sort_tmp_1[i + tid];
//			//printf("\nsort tmp: %d", sort_tmp_0[i + tid]);
//		}
//	}
//	__syncthreads();
//
//	//for (int i = 0; i < num_element; i++) {
//	//	printf("\n%d", sort_tmp[i]);
//	//}
//}
//int find_min(unsigned int * const src_array,
//	int * const list_indexes,
//	const int num_lists,
//	const int num_elements_per_list)
//{
//	int min_val = 0xFFFFFFFF;
//	int min_idx = 0;
//	// Iterate over each of the lists
//	for (int i = 0; i < num_lists; i++)
//	{
//		// If the current list has already been emptied
//		// then ignore it
//		if (list_indexes[i] < num_elements_per_list)
//		{
//			const int src_idx = i + (list_indexes[i] * num_lists);
//			const int data = src_array[src_idx];
//			if (data <= min_val)
//			{
//				min_val = data;
//				min_idx = i;
//			}
//		}
//	}
//	list_indexes[min_idx]++;
//	return min_val;
//}
//__device__ void copy_data_to_shared(unsigned int * data,
//	unsigned int * sort_tmp,
//	const int num_lists,
//	const int num_elements,
//	const int tid)
//{
//	// Copy data into temp store
//	for (int i = 0; i < num_elements; i += num_lists)
//	{
//		sort_tmp[i + tid] = data[i + tid];
//	}
//	__syncthreads();
//}
//__device__ void merge_array6(unsigned int * src_array,
//	unsigned int * dest_array,
//	const int num_lists,
//	const int num_elements,
//	const int tid)
//{
//	const int num_elements_per_list = (num_elements / num_lists);
//	
//	__shared__ int list_indexes[32];
//
//	list_indexes[tid] = 0;
//	// Wait for list_indexes[tid] to be cleared
//	__syncthreads();
//	// Iterate over all elements
//	for (int i = 0; i < num_elements; i++)
//	{
//		// Create a value shared with the other threads
//		__shared__ int min_val;
//		__shared__ int min_tid;
//		// Use a temp register for work purposes
//		int data;
//		// If the current list has not already been
//		// emptied then read from it, else ignore it
//		if (list_indexes[tid] < num_elements_per_list)
//		{
//			// Work out from the list_index, the index into
//			// the linear array
//			const int src_idx = tid + (list_indexes[tid] * num_lists);
//			// Read the data from the list for the given
//			// thread
//			data = src_array[src_idx];
//		}
//		else
//		{
//			data = 0xFFFFFFFF;
//		}
//		// Have thread zero clear the min values
//		if (tid == 0)
//		{
//			// Write a very large value so the first
//			// thread thread wins the min
//			min_val = 0xFFFFFFFF;
//			min_tid = 0xFFFFFFFF;
//		}
//		// Wait for all threads
//		__syncthreads();
//		// Have every thread try to store it’s value into
//		// min_val. Only the thread with the lowest value
//		// will win
//
//		atomicMin(&min_val, data);
//
//		// Make sure all threads have taken their turn.
//		__syncthreads();
//		// If this thread was the one with the minimum
//		if (min_val == data)
//		{
//			// Check for equal values
//			// Lowest tid wins and does the write
//
//			atomicMin(&min_tid, tid);
//		}
//		// Make sure all threads have taken their turn.
//		__syncthreads();
//		// If this thread has the lowest tid
//		if (tid == min_tid)
//		{
//			// Incremene the list pointer for this thread
//			list_indexes[tid]++;
//			// Store the winning value
//			dest_array[i] = data;
//		}
//	}
//}
//
//__global__ void gpu_sort_array_array(unsigned int * data, const int num_lists, const int num_elements)
//{
//	int tid = (blockIdx.x * blockDim.x) + threadIdx.x;
//
//	__shared__ unsigned int sort_tmp[64];
//	__shared__ unsigned int sort_tmp_1[64];
//
//	copy_data_to_shared(data, sort_tmp, num_lists,
//		num_elements, tid);
//	radixSortParralel(sort_tmp, num_lists, num_elements,
//		tid, sort_tmp_1);
//	merge_array6(sort_tmp, data, num_lists,
//		num_elements, tid);
//
//}
//
//void print(unsigned int* arr, int n)
//{
//	for (int i = 0; i < n; i++)
//		cout << (double)arr[i] / 100 << endl;
//}
//
//void menu() {
//	cout << "Izberi: " << endl;
//	cout << " 0 Exit " << endl;
//	cout << " 1 Serial sort" << endl;
//	cout << " 2 Parrallel sort" << endl;
//	cout << " 3 Thrust Parrallel sort" << endl;
//}
//
////void WriteFile(vector<double> vector) {
////	ofstream myfile;
////	myfile.open("out.txt");
////	for each (double var in vector)
////	{
////		myfile << var << " ";
////	}
////	myfile.close();
////}
//
//vector<unsigned int> ReadFile(string filepath) {
//	vector<unsigned int> vector;
//	ifstream inputFile(filepath);
//	double value;
//
//	//read the elements in the file into a vector  
//	while (inputFile >> value) {
//		vector.push_back(value * 100);
//	}
//	inputFile.close();
//	return vector;
//}
//
//int main(int argc, char * argv[]) {
//
//	int blockSize = 32;
//	int numBlocks = (SIZE + blockSize - 1) / blockSize;
//
//	cudaEvent_t start, stop;
//	cudaEventCreate(&start);
//	cudaEventCreate(&stop);
//	float totalTime = 0;
//
//	std::ifstream fin(argv[1]);
//	//std::vector<unsigned int> host_double(1);
//
//	unsigned int host_data[SIZE];
//
//	// generate 32M random numbers serially
//	
//
//	if (!fin) {
//		cerr << "Datoteka ne obstaja " << argv[1] << endl;
//	}
//	else {
//		cout << "Datoteka obstaja " << argv[1] << endl;
//		//host_double = ReadFile(argv[1]);
//
//		//WriteFile(host_double);
//	}
//	
//
//	while (true) {
//		thrust::host_vector<unsigned int> host_double(1 << 20);
//		std::generate(host_double.begin(), host_double.end(), rand);
//		int size = host_double.size();
//
//		int izbiraAlg;
//		menu();
//		cin >> izbiraAlg;
//		if (izbiraAlg == 0) {
//			return 0;
//		}
//		else if (izbiraAlg == 1) {
//
//			//Serial radix sort
//
//			unsigned int* da = host_double.data();
//			int n = sizeof(host_double) / sizeof(host_double[0]);
//			//auto t1 = high_resolution_clock::now();
//			high_resolution_clock::time_point start = high_resolution_clock::now();
//			serialRadixSort(da, size);
//			high_resolution_clock::time_point stop = high_resolution_clock::now();
//			auto duration = duration_cast <milliseconds> (stop - start).count();
//			cout << endl << "Serial radix sort " << endl;
//			//print(da, size);
//			for (int i = 0; i < size - 1; i++) {
//				if (host_double[i] > host_double[i + 1]) {
//					printf("sort error at, hdata[%d] = %d, hdata[%d] = %d\n", i, host_double[i], i + 1, host_double[i + 1]);
//					return 1;
//				}
//			}
//			cout << "Total time: " << duration << "ms" << endl;
//			cout << "Success" << endl;
//
//		}
//		else if (izbiraAlg == 2) {
//
//			//Parrallel radix sort
//			cout << endl;
//			cout << "Parrallel sort " << endl;
//			std::copy(host_double.begin(), host_double.end(), host_data);
//			/*cout << "Not sorted!" << endl;
//			print(host_data, size);*/
//
//			unsigned int* ddata;
//			unsigned int* ddata_tmp;
//			//
//
//			//for (int lcount = 0; lcount < LOOPS; lcount++) {
//			//	unsigned int range = 1U << UPPER_BIT;
//			//	//for (int i = 0; i < SIZE; i++) host_data[i] = rand() % range;
//				//cudaMemcpyToSymbol(device_data, host_data, SIZE * sizeof(unsigned int));
//			cudaMalloc((void**)&ddata_tmp, SIZE * sizeof(unsigned int));
//			cudaMalloc((void**)&ddata, SIZE * sizeof(unsigned int));
//			cudaMemcpy(ddata, host_data, SIZE * sizeof(unsigned int), cudaMemcpyHostToDevice);
//
//			high_resolution_clock::time_point start = high_resolution_clock::now();
//			//radixSortParralel << <1, 10 >> >(ddata, ddata_tmp);
//			//radixSort << <(numBlocks, blockSize) >> > (ddata, ddata_tmp);
//			gpu_sort_array_array << < 1, 2 >> > (ddata, 32, 64);
//			high_resolution_clock::time_point stop = high_resolution_clock::now();
//			host_data[SIZE];
//
//			cudaMemcpy(host_data, ddata, SIZE * sizeof(unsigned int), cudaMemcpyDeviceToHost);
//
//			//cudaMemcpyFromSymbol(host_data, device_data, SIZE * sizeof(unsigned int));
//			auto durationParralel = duration_cast <milliseconds> (stop - start).count();
//			//	/*for (int i = 0; i < SIZE - 1; i++) {
//			//	if (host_data[i] > host_data[i + 1]) {
//			//	printf("sort error at loop %d, hdata[%d] = %d, hdata[%d] = %d\n", lcount, i, host_data[i], i + 1, host_data[i + 1]);
//			//	return 1;
//			//	}
//			//	}*/
//			//	cout << " Sorted data: " << endl;
//			//	for (int i = 0; i < size; i++) {
//			//		cout << i << "     " <<(double)host_data[i] / 100 << endl;
//			//	}
//			//	cout << " Time parralel sort: " << durationParralel << endl;
//			//}
//			cout << endl << "Sorted" << endl;
//			print(host_data, size);
//
//		}
//		else {
//			//thrust
//
//			cout << "Parralel Radix sort Thrust" << endl;
//
//			// transfer data to the device
//			thrust::device_vector<int> d_vec = host_double;
//
//			high_resolution_clock::time_point start = high_resolution_clock::now();
//			// sort data on the device (846M keys per second on GeForce GTX 480)
//			thrust::sort(d_vec.begin(), d_vec.end());
//			high_resolution_clock::time_point stop = high_resolution_clock::now();
//
//			// transfer data back to host
//			thrust::copy(d_vec.begin(), d_vec.end(), host_double.begin());
//
//			auto durationParralel = duration_cast <milliseconds> (stop - start).count();
//			cout << "Time: " << durationParralel << "ms" << endl;
//			/*for each (unsigned int var in host_double)
//			{
//				double tmp = (double)var / 100;
//				cout << tmp << endl;
//			}*/
//		}
//	}
//
//
//	return 0;
//}


