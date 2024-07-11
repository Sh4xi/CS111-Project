//PROGRAMMING PROJECT IN CS111
//PROGRAM USED: C LANGUAGE
//GROUP MEMBERS:
// *ERJAY CLEOFE
// *IGNACIO TABUG III
// *IAN BANARES
// *ARAGORN MUNCAL

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define MAXRANGE 10000000	//MAX RANGE OF RANDOMLY GENERATED NUMBERS			
#define MAXSIZE 10000000	//MAX SIZE OF AN ARRAY ELEMENT

//Prototype function declarations
void generate_random(int N, int array[]);
void generate_sorted(int N, int X, int array[]);

double insertionSort(int array[], int N);
void bubbleSort(int arr[], int n);
void selectionSort(int arr[], int n);

void mergeSort(int array[], int low, int high);
void merge(int arr[], int low, int mid, int high);

void quickSort(int array[], int low, int high);
int partition(int array[], int low, int high);

void heapSort(int array[], int N);
void heapify(int arr[], int n, int i);

void swap(int *a, int *b);
void printArr(int N, int array[]);
void fprintArr(int N, int array[], FILE *fptr);


int main(void) { 				//MAIN FUNCTION
	clock_t start, end;			//variable declarations
	double cpu_time_used;
	int N, X, i;
	
	int *arrRand = malloc(MAXSIZE*sizeof(int));			//initialize the array on the Heap
	int *arrSort = malloc(MAXSIZE*sizeof(int));     	//int arrRand[MAXSIZE], arrSort[MAXSIZE];			
														//arrRand stores array of random inputs; arrSORT sorted in increasing order
	
	printf("\nPLEASE ENTER THE ARRAY SIZE: ");	//prompt user input N
	do{										//loop makes sure only positive input
		scanf("%d", &N);					//number of integers as input used in sorting
		if(N < 1)
			printf("INVALID INPUT! POSITIVE INTERGER ONLY: ");
	}while(N < 1);		
	generate_random(N, arrRand);
	printf("\tGENERATED ARRAY OF RANDOM NUMBERS: \n");
		//printArr(N, arrRand);
	
	
	printf("\nENTER A POSITIVE INTEGER FOR X (Interval for sorted array): ");						//prompt user input X
	do{										//loop makes sure only positive input
		scanf("%d", &X);					//number to be used as increment for sorted array
		if(X < 1)
			printf("INVALID INPUT! POSITIVE INTERGER ONLY:  ");
	}while(X < 1);
	generate_sorted(N, X, arrSort);
	printf("\tGENERATED ARRAY OF SORTED NUMBERS: \n");
		//printArr(N, arrSort);

	printf("\n\nPROGRAM RESULTS:\n");

	//make 6 copies of the contents of arrRand, that will be sorted in respective arrays
	int *iSort = malloc(MAXSIZE*sizeof(int));
	int *bSort = malloc(MAXSIZE*sizeof(int));
	int *sSort = malloc(MAXSIZE*sizeof(int));
	int *mSort = malloc(MAXSIZE*sizeof(int));
	int *qSort = malloc(MAXSIZE*sizeof(int));
	int *hSort = malloc(MAXSIZE*sizeof(int));
	for(i=0; i<N; i++){
		iSort[i] = arrRand[i];
		bSort[i] = arrRand[i];
		sSort[i] = arrRand[i];
		mSort[i] = arrRand[i];
		qSort[i] = arrRand[i];
		hSort[i] = arrRand[i];		
	}
	
	double insertionS_time;
	start = clock();
		insertionSort(iSort, N);
	end = clock();
	insertionS_time = ((double)(end-start)) / CLOCKS_PER_SEC;
	printf("\nArray sorted by INSERTION SORT:\n");
	//printArr(N, iSort);
	printf("\tTime Insertion Sort took for N =%d, Random Input : %2.10lf seconds \n", N, insertionS_time);
	
	double bubbleS_time;
	start = clock();
		bubbleSort(bSort, N);
	end = clock();
	bubbleS_time = ((double)(end-start)) / CLOCKS_PER_SEC;
	printf("\nArray sorted by BUBBLE SORT:\n");
	//printArr(N, bSort);
	printf("\tTime Bubble Sort took for N =%d, Random Input : %2.10lf seconds \n", N, bubbleS_time);

	double selectionS_time;
	start = clock();
		selectionSort(sSort, N);
	end = clock();
	selectionS_time = ((double)(end-start)) / CLOCKS_PER_SEC;
	printf("\nArray sorted by SELECTION SORT:\n");
	//printArr(N, sSort);
	printf("\tTime Selection Sort took for N =%d, Random Input : %2.10lf seconds \n", N, selectionS_time);

	double mergeS_time;
	start = clock();
		mergeSort(mSort, 0, N-1);
	end = clock();
	mergeS_time = ((double)(end-start)) / CLOCKS_PER_SEC;
	printf("\nArray sorted by MERGE SORT:\n");
	//printArr(N, mSort);
	printf("\tTime Merge Sort took for N =%d, Random Input     : %2.10lf seconds \n", N, mergeS_time);

	double quickS_time;
	start = clock();
		quickSort(qSort, 0, N-1);
	end = clock();
	quickS_time = ((double)(end-start)) / CLOCKS_PER_SEC;
	printf("\nArray sorted by QUICK SORT:\n");
	//printArr(N, qSort);
	printf("\tTime Quick Sort took for N=%d, Random Input: %2.10lf seconds \n", N, quickS_time);

	double heapS_time;
	start = clock();
		heapSort(hSort, N);
	end = clock();
	heapS_time = ((double)(end-start)) / CLOCKS_PER_SEC;
	printf("\nArray sorted by HEAP SORT:\n");
	//printArr(N, hSort);
	printf("\tTime Heap  Sort took for N =%d, Random Input     : %2.10lf seconds \n", N, heapS_time);
	printf("\n\n");
	
	//make 6 copies of the contents of arrSort, that will be sorted in respective arrays
	int *iSorted = (int*)malloc(MAXSIZE*sizeof(int));
	int *bSorted = (int*)malloc(MAXSIZE*sizeof(int));
	int *sSorted = (int*)malloc(MAXSIZE*sizeof(int));
	int *mSorted = (int*)malloc(MAXSIZE*sizeof(int));
	int *qSorted = (int*)malloc(MAXSIZE*sizeof(int));
	int *hSorted = (int*)malloc(MAXSIZE*sizeof(int));
	for(i=0; i<N; i++){
		iSorted[i] = arrSort[i];
		bSorted[i] = arrSort[i];
		sSorted[i] = arrSort[i];
		mSorted[i] = arrSort[i];
		qSorted[i] = arrSort[i];
		hSorted[i] = arrSort[i];		
	}
	//Find running time of array that is already sorted:
	double insertionS_time2;
	insertionS_time2 = insertionSort(iSorted, N);
	printf("\tTime Insertion Sort took for N =%d, Sorted Input : %2.10lf seconds \n", N, insertionS_time2);
	
	double bubbleS_time2;
	bubbleS_time2 = insertionSort(bSorted, N);
	printf("\tTime Bubble Sort took for N =%d, Sorted Input    : %2.10lf seconds \n", N, bubbleS_time2);

	double selectionS_time2;
	selectionS_time2 = insertionSort(sSorted, N);
	printf("\tTime Selection Sort took for N =%d, Sorted Input : %2.10lf seconds \n", N, selectionS_time2);

	double mergeS_time2;
	start = clock();
	mergeSort(mSorted, 0, N-1);
	end = clock();
	mergeS_time2 = ((double)(end-start)) / CLOCKS_PER_SEC;
	printf("\tTime Merge Sort took for N =%d, Sorted Input     : %2.10lf seconds \n", N, mergeS_time2);
	
	double quickS_time2;
	start = clock();
	quickSort(qSorted, 0, N-1);
	end = clock();
	quickS_time2 = ((double)(end-start)) / CLOCKS_PER_SEC;
	printf("\tTime Quick Sort took for N =%d, Sorted Input     : %2.10lf seconds \n", N, quickS_time2);
	
	double heapS_time2;
	start = clock();
	heapSort(hSorted, N);
	end = clock();
	heapS_time2 = ((double)(end-start)) / CLOCKS_PER_SEC;
	printf("\tTime Heap  Sort took for N =%d, Sorted Input     : %2.10lf seconds \n", N, heapS_time2);

	
	//FILE PRINTING
	FILE *fptr;
	fptr = fopen("DAA-project-Array_Values.txt","w");
	if(fptr == NULL){
		printf("ERROR! no output File was produced.");
		return 1;
	}
	
	fprintf(fptr, "Time Insertion Sort took for N = %6d,     Random Input:\t%3.10lf seconds \n", N, insertionS_time);
	fprintf(fptr, "Time Bubble Sort took for    N = %6d,     Random Input:\t%3.10lf seconds \n", N, bubbleS_time);
	fprintf(fptr, "Time Selection Sort took for N = %6d,     Random Input:\t%3.10lf seconds \n", N, selectionS_time);
	fprintf(fptr, "Time Merge Sort took for     N = %6d,     Random Input:\t%3.10lf seconds \n", N, mergeS_time);
	fprintf(fptr, "Time Quick Sort took for     N = %6d,     Random Input:\t%3.10lf seconds \n", N, quickS_time);
	fprintf(fptr, "Time Heap  Sort took for     N = %6d,     Random Input:\t%3.10lf seconds \n", N, heapS_time);
	fprintf(fptr,"\n\n");
	fprintf(fptr, "Time Insertion Sort took for N = %6d,     Sorted Input:\t%3.10lf seconds \n", N, insertionS_time2);
	fprintf(fptr, "Time Bubble Sort took for    N = %6d,     Sorted Input:\t%3.10lf seconds \n", N, bubbleS_time2);
	fprintf(fptr, "Time Selection Sort took for N = %6d,     Sorted Input:\t%3.10lf seconds \n", N, selectionS_time2);
	fprintf(fptr, "Time Merge Sort took for     N = %6d,     Sorted Input:\t%3.10lf seconds \n", N, mergeS_time2);
	fprintf(fptr, "Time Quick Sort took for     N = %6d,     Sorted Input:\t%3.10lf seconds \n", N, quickS_time2);
	fprintf(fptr, "Time Heap  Sort took for     N = %6d,     Sorted Input:\t%3.10lf seconds \n", N, heapS_time2);

	fprintf(fptr,"\n\nORIGINAL ARRAY WITH RANDOMY GENERATE VALUES:\n");
	fprintArr(N,arrRand,fptr);
	fprintf(fptr,"\n\nSORTED ARRAY WITH RANDOMLY GENERATED VALUES:\n");
	fprintArr(N,iSort,fptr);
	fprintf(fptr,"\n\nArray with increasing sequence based on X:\n");
	fprintArr(N,arrSort,fptr);
	fprintf(fptr,"\n\nArray Sorted by [BUBBLE SORT]:\n");
	fprintArr(N,bSorted,fptr);
	fprintf(fptr,"\n\nArray Sorted by [SELECTION SORT]:\n");
	fprintArr(N,sSorted,fptr);
	fprintf(fptr,"\n\nArray Sorted by [MERGE SORT]:\n");
	fprintArr(N,mSorted,fptr);
	fprintf(fptr,"\n\nArray Sorted by [QUICK SORT]:\n");
	fprintArr(N,qSorted,fptr);
	fprintf(fptr,"\n\nArray Sorted by [HEAP SORT]:\n");
	fprintArr(N,hSorted,fptr);
	
	fclose(fptr);
	
	
	//free all the allocated memory used earlier
	free(arrRand);
	free(iSort);free(bSort);free(sSort);free(mSort);free(qSort);free(hSort);
	free(arrSort);
	free(iSorted);free(bSorted);free(sSorted);free(mSorted);free(qSorted);free(hSorted);
	
	return 0;
}	//end of main



void generate_random(int N, int array[]){				//generates array of random nums
	int i;
	srand(time(NULL));		//sets the seed of rand, by predefined function time (seconds elapsed since Jan. 1, 1970)
	unsigned long x;

	for(i=0; i<N; i++){
		x = rand();				
		x <<= 15; 				//shift by 15 bits, to increase max value
		x ^= rand();			//takes 2 random numbers and performs XOR on every bit. result is 1 if 2 bits are different.
		x %= MAXRANGE+1;			//perform modulo by 1M+1, so number can only be between 0 to 1M.
				
		array[i] = x ;
	}
	
	/*	NOTES	
	rand can give at least 15 bits, and at max the value 32767 on windows
	shift by 15bits, XOR the results, then modulo by 1000000
	*/
	return;
}
void generate_sorted(int N, int X, int array[]){			//generates an array already sorted in increasing order
	int i;
	for (i=0; i<N; i++){		
		array[i] = N+(1+i)*X ;
	}
	return;
}

double insertionSort(int array[], int N){			//function for insertion sort
	clock_t start, end;
	double cpu_time_used;
	
	start = clock();
	int i, j, key;
	for(j=1; j<N; j++){								//
		key = array[j];								//
		i = j-1;									//
		while ((i>=0) && (array[i]>key)){			//
			array[i+1] = array[i];					//
			i = i-1;								//
		}											//
		array[i+1] = key;							//
	}
	end = clock();
	
	cpu_time_used = ((double)(end-start)) / CLOCKS_PER_SEC;
	return cpu_time_used;
}

// Selection sort algorithm
void selectionSort(int arr[], int n) {
    int minIndex, temp;
    for (int i = 0; i < n - 1; i++) {
        minIndex = i;
        for (int j = i + 1; j < n; j++) {
            if (arr[j] < arr[minIndex]) {
                minIndex = j;
            }
        }
        temp = arr[minIndex];
        arr[minIndex] = arr[i];
        arr[i] = temp;
    }
}

// Bubble sort algorithm
void bubbleSort(int arr[], int n) {
    int temp;
    for (int i = 0; i < n - 1; i++) {
        for (int j = 0; j < n - i - 1; j++) {
            if (arr[j] > arr[j+1]) {
                temp = arr[j];
                arr[j] = arr[j+1];
                arr[j+1] = temp;
            }
        }
    }
}

void mergeSort(int array[], int low, int high){
	if (low < high){
		int mid = (low+high)/2;
		
		mergeSort(array, low, mid);
		mergeSort(array, mid+1, high);
		
		merge(array, low, mid, high);		
	}
	return;
}
void merge(int arr[], int low, int mid, int high){	//function to merge subarrays for MergeSort
	
	//find length of temp arrays to be used
	int n1 = mid - low+1;
	int n2 = high - mid;
	
	//create temp array
//	int Left[n1], Right[n2];
	int *Left = malloc(n1*sizeof(int));
	int *Right = malloc(n2*sizeof(int));
	
	
	//copy data to temp array
	int i, j, k;
	for(i=0; i<n1; i++)
		Left[i] = arr[low+i];
	for(j=0; j<n2;j++)
		Right[j] = arr[mid+1+j];
		
	//merge the temp arrays (with leftovers)
	i=0, j=0, k=low;
	while((i<n1)&&(j<n2)){
		if (Left[i] <= Right[j]){
			arr[k] = Left[i];
			i++;
		}
		else {
			arr[k] = Right[j];
			j++;
		}
		k++;
	}
	
	//copy the remaining elements 
	while (i<n1){
		arr[k] = Left[i];
		i++; k++;
 	}
 	while (j<n2){
 		arr[k] = Right[j];
 		j++; k++;
	 }
	 
	 
	free(Left);
	free(Right);
}

void quickSort(int arr[], int low, int high) {
    if (low < high) {
        int pivot = partition(arr, low, high);
        quickSort(arr, low, pivot - 1);
        quickSort(arr, pivot + 1, high);
    }
}

int partition(int arr[], int low, int high) {
    int pivot = arr[high];
    int i = (low - 1);

    for (int j = low; j <= high - 1; j++) {
        if (arr[j] < pivot) {
            i++;
            swap(&arr[i], &arr[j]);
        }
    }
    swap(&arr[i + 1], &arr[high]);
    return (i + 1);
}

void heapSort(int array[], int N){
	int i;
	for(i=N/2-1; i>=0; i--)				//build maxheap
		heapify(array, N, i);
		
	//Heap sort
	for(i=N-1; i>=0; i--){
		swap(&array[0], &array[i]);
		heapify(array, i, 0);			//heapify root element to get highest element at root again
	}
	return;
}
void heapify(int arr[], int n, int i){
	//Find largest in root, left and right child
	int largest = i;
	int left = 2*i +1;
	int right = 2*i +2;
	
	if((left<n) && (arr[left] > arr[largest]))
		largest = left;
		
	if((right<n) && (arr[right] > arr[largest]))
		largest = right;
		
	if(largest != i){					//swap and continue heapifying if root is not largest
		swap(&arr[i], &arr[largest]);
		heapify(arr, n, largest);
	}
	return;
}


void swap(int *a, int *b){								//swaps two elements
	int t = *a;
	*a = *b;
	*b = t;
}

void printArr(int N, int array[]){						//Prints the array on screen
	int i, j, k;
	printf("\t\t");
		for(i=0; i<17; i++)
			printf("------");							//prints top border
	int rows = (N-1)/10;

	i=0;
		for(j=0; j<=rows; j++){
			printf("\n\t\t|");
				for (k=0; (k<10)&&(i<N); i++,k++){
					printf(" %8d ",array[i]);
				}
				if (k != 10){							//makes sure proper closing border.
					for(;k<10;k++)
						printf("          ");
				}
			printf("|");
		}
	printf("\n\t\t");
		for(i=0; i<17; i++)
			printf("------");							//prints bottom border
		printf("\n");
	
	return;
}
void fprintArr(int N, int array[], FILE *fptr){						//prints array on a file
	int i, j, k;
	fprintf(fptr,"\t\t");
		for(i=0; i<17; i++)
			fprintf(fptr,"------");							//prints top border
	int rows = (N-1)/10;

	i=0;
		for(j=0; j<=rows; j++){
			fprintf(fptr,"\n\t\t|");
				for (k=0; (k<10)&&(i<N); i++,k++){
					fprintf(fptr," %8d ",array[i]);
				}
				if (k != 10){							//makes sure proper closing border.
					for(;k<10;k++)
						fprintf(fptr,"          ");
				}
			fprintf(fptr,"|");
		}
	fprintf(fptr,"\n\t\t");
		for(i=0; i<17; i++)
			fprintf(fptr,"------");							//prints bottom border
		fprintf(fptr,"\n");
	
	return;	
}
