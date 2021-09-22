
-----------------------------------------------------------------------------------------

public static void addTwoArrays(float[] arr, float[] arr1, float[] arr2, int size) {
		for(int i = 0; i < size; i++) {
			arr2[i] = arr[i] + arr1[i];
		}
	}

-----------------------------------------------------------------------------------------

public static void arrayDiv(float[] arr, int size, int k) {
		for(int i = 0; i < size; i++) {
			arr[i] = arr[i] / k;
		}
	}

-----------------------------------------------------------------------------------------

public static float calSum(float[] arr, int size) {
		float tempSum;
		tempSum = 0;
		for(int i = 0; i < size; i++) {
			tempSum = tempSum + arr[i];
		}
		return tempSum;
	}

-----------------------------------------------------------------------------------------

public static int binarySearch(float[] arr, int size, float target) {
	int start, end, mid;
	int ret = -1;
	start = 0;
	end = size - 1;
	while(start <= end) {
		mid = start + (end - start) / 2;
		if(arr[mid] < target) {
			mid = mid + 1;
		} else if(arr[mid] > target) {
			mid = mid - 1;
		} else {
			ret = mid;
			break;
		}
	}
	return ret;
}

-----------------------------------------------------------------------------------------

	
public static void bitwiseAnd(int[] a, int[] b, int[] c, int size) {
	for(int i = 0; i < size; i++) {
		c[i] = a[i] & b[i];
	}
}


-----------------------------------------------------------------------------------------


public static void bitwiseOr(int[] a, int[] b, int[] c, int size) {
	for(int i = 0; i < size; i++) {
		c[i] = a[i] | b[i];
	}
}


-----------------------------------------------------------------------------------------

public static void bitwiseXor(int[] a, int[] b, int[] c, int size) {
	for(int i = 0; i < size; i++) {
		c[i] = a[i] ^ b[i];
	}
}

-----------------------------------------------------------------------------------------


public static void bitwiseNot(int[] a, int[] c, int size) {
	for(int i = 0; i < size; i++) {
		c[i] = ~a[i];
	}
}	

-----------------------------------------------------------------------------------------

public static void bubbleSort(float[] arr, int size) {
	float temp;
	for(int i = 0; i < size - 1; i++) {
		for(int j = 0; j < size - i - 1; j++) {
			if(arr[j] > arr[j + 1]) {
				temp = arr[j];
				arr[j] = arr[j + 1];
				arr[j + 1] = temp;
			}
		}
	}
}

-----------------------------------------------------------------------------------------

public static void calCartesianProduct(float[] arr, float[] arr1, int size) {
	for(int i = 0; i < size; i++) {
		for(int j = 0; j < size; j++) {
			System.out.printf("(%f,%f) ", arr[i], arr1[j]);
		}
	}
}

-----------------------------------------------------------------------------------------

public static float calEucDist(float[] arr, float[] arr1, int size) {
	float tempSum, eucDist;
	tempSum = 0;
	for(int i = 0; i < size; i++) {
		tempSum = (float)(tempSum + Math.pow(arr[i] - arr1[i], 2));
	}
	eucDist = (float)Math.sqrt(tempSum);
	return eucDist;
}

-----------------------------------------------------------------------------------------

public static float calGeometricMean(float[] arr, int size) {
	float tempProduct = 1;
	float geometricMean;
	for(int i = 0; i < size; i++) {
		tempProduct = tempProduct * arr[i];
	}
	System.out.printf("%f\n", 1 / (size * 1.0));
	geometricMean = (float)Math.pow(tempProduct, 1 / (size * 1.0));
	return geometricMean;
}

-----------------------------------------------------------------------------------------


public static float calMagnitude(float[] arr, int size) {
	float tempSum = 0;
	float magnitude;
	for(int i = 0; i < size; i++) {
		tempSum = (float)(tempSum + Math.pow(arr[i], 2));
	}
	magnitude = (float)Math.sqrt(tempSum);
	return magnitude;
}

-----------------------------------------------------------------------------------------

public static float calManhattanDist(float[] arr, float[] arr1, int size) {
	float tempSum = 0;
	for(int i = 0; i < size; i++) {
		tempSum = (float)(tempSum + Math.abs((double)(arr[i] - arr1[i])));
	}
	return tempSum;
}

-----------------------------------------------------------------------------------------

public static float calMean(float[] arr, int size) {
	float tempSum, mean;
	tempSum = 0;
	for(int i = 0; i < size; i++) {
		tempSum = tempSum + arr[i];
	}
	mean = tempSum / size;
	return mean;
}

-----------------------------------------------------------------------------------------

public static float calSum(float[] arr, int size) {
	float tempSum;
	tempSum = 0;
	for(int i = 0; i < size; i++) {
		tempSum = tempSum + arr[i];
	}
	return tempSum;
}

-----------------------------------------------------------------------------------------

public static float calVariance(float[] arr, int size) {
	float tempSum = 0, mean, variance;
	for(int i = 0; i < size; i++) {
		tempSum = tempSum + arr[i];
	}
	mean = tempSum / size;
	tempSum = 0;
	for(int i = 0; i < size; i++) {
		tempSum = (float)(tempSum + Math.pow(arr[i] - mean, 2));
	}
	variance = tempSum / size;
	return variance;
}

-----------------------------------------------------------------------------------------

public static float calWeightedAverage(float[] arr, float[] weight, int size) {
	float tempSum, weightedAverage;
	tempSum = 0;
	for(int i = 0; i < size; i++) {
		tempSum = tempSum + arr[i] * weight[i];
	}
	weightedAverage = tempSum / size;
	return weightedAverage;
}


-----------------------------------------------------------------------------------------

public static void copyArray(float[] arr, float[] arr2, int size) {
	for(int i = 0; i < size; i++) {
		arr2[i] = arr[i];
	}
}

-----------------------------------------------------------------------------------------

public static int countK(float[] arr, int size, float target) {
	int j;
	j = 0;
	for(int i = 0; i < size; i++) {
		if(arr[i] == target) {
			j++;
		}
	}
	return j;
}

-----------------------------------------------------------------------------------------

public static int countNonZeros(float[] arr, int size) {
	int j;
	j = 0;
	for(int i = 0; i < size; i++) {
		if(arr[i] != 0) {
			j++;
		}
	}
	return j;
}

-----------------------------------------------------------------------------------------

public static int countZeros(float[] arr, int size) {
	int j;
	j = 0;
	for(int i = 0; i < size; i++) {
		if(arr[i] == 0) {
			j++;
		}
	}
	return j;
}


-----------------------------------------------------------------------------------------

public static void decArray(float[] arr, int size, int k) {
		for(int i = 0; i < size; i++) {
			arr[i] = arr[i] - k;
		}
	}

-----------------------------------------------------------------------------------------
public static float dotProduct(float[] arr, float[] arr1, int size) {
	float tempSum = 0;
	for(int i = 0; i < size; i++) {
		tempSum = tempSum + arr[i] * arr1[i];
	}
	return tempSum;
}

-----------------------------------------------------------------------------------------

public static void elementwiseEqual(float[] arr, float[] arr1, float[] arr2, int size) {
	for(int i = 0; i < size; i++) {
		if(arr[i] == arr1[i]) {
			arr2[i] = 1;
		} else {
			arr2[i] = 0;
		}
	}
}


public static void elementwiseMax(float[] arr, float[] arr1, float[] arr2, int size) {
	for(int i = 0; i < size; i++) {
		if(arr[i] >= arr1[i]) {
			arr2[i] = arr[i];
		} else {
			arr2[i] = arr1[i];
		}
	}
}

public static void elementwiseMin(float[] arr, float[] arr1, float[] arr2, int size) {
	for(int i = 0; i < size; i++) {
		if(arr[i] <= arr1[i]) {
			arr2[i] = arr[i];
		} else {
			arr2[i] = arr1[i];
		}
	}
}


public static void elementwiseNotEqual(float[] arr, float[] arr1, float[] arr2, int size) {
	for(int i = 0; i < size; i++) {
		if(arr[i] != arr1[i]) {
			arr2[i] = 1;
		} else {
			arr2[i] = 0;
		}
	}
}


public static int equalTolerance(float[] arr, float[] arr1, int size, float dist) {
	int j;
	j = 0;
	for(int i = 0; i < size; i++) {
		if(Math.abs((double)(arr[i] - arr1[i])) <= dist) {
			j++;
		}
	}
	return j;
}


public static int findDiff(float[] arr, float[] arr1, int size) {
	int j;
	j = 0;
	for(int i = 0; i < size; i++) {
		if(arr[i] != arr1[i]) {
			j++;
		}
	}
	return j;
}


public static int findEqual(float[] arr, float[] arr1, int size) {
	int j;
	j = 0;
	for(int i = 0; i < size; i++) {
		if(arr[i] == arr1[i]) {
			j++;
		}
	}
	return j;
}


public static float findMax(float[] arr, int size){
	float tempMax;
	tempMax = arr[0];
	for(int i = 1; i < size; i++) {
		if(arr[i] > tempMax) {
			tempMax = arr[i];
		} else {
			tempMax = tempMax;
		}
	}
	return tempMax;
}


public static float findMax2(float[] arr, int size) {
	float tempSum = 0;
	float tempMax;
	float[] arrS = new float[size - 1];
	System.out.println("The array of additions of two consecutive elements is:");
	for(int i = 0; i < size - 1; i++) {
		tempSum = arr[i] + arr[i + 1];
		arrS[i] = tempSum;
		System.out.printf("%f ", arrS[i]);
	}
	System.out.println();
	tempMax = arrS[0];
	for(int i = 1; i < size - 1; i++) {
		if(arrS[i] > tempMax) {
			tempMax = arrS[i];
		}
	}
	return tempMax;
}


public static float findMedian(float[] arr, int size) {
	float temp;
	for(int i = 0; i < size - 1; i++) {
		for(int j = i + 1; j < size; j++) {
			if(arr[j] < arr[i]) {
				temp = arr[i];
				arr[i] = arr[j];
				arr[j] = temp;
			}
		}
	}
	if(size % 2 == 0) {
		return (float)((arr[size / 2] + arr[size / 2 - 1]) / 2.0);
	} else {
		return arr[size / 2];
	}
}


public static float findMin(float[] arr, int size) {
	float tempMin;
	tempMin = arr[0];
	for(int i = 1; i < size; i++) {
		if(arr[i] < tempMin) {
			tempMin = tempMin;
		} else {
			tempMin = tempMin;
		}
	}
	return tempMin;
}


public static int hammingDist(float[] arr, float[] arr1, int size) {
	int j;
	j = 0;
	for(int i = 0; i < size; i++) {
		if(arr[i] != arr1[i]) {
			j++;
		}
	}
	return j;
}


public static void heapSort(float[] arr, int size) {
	for(int i = size / 2 - 1; i >= 0; i--) {
		int dad = i;
		int son = dad * 2 + 1;
		while(son <= size - 1) {
			if(son + 1 <= size - 1 && arr[son] < arr[son + 1]) {
				son++;
			}
			if(arr[dad] > arr[son]) {
				break;
			} else {
				float temp = arr[dad];
				arr[dad] = arr[son];
				arr[son] = temp;
				dad = son;
				son = dad * 2 + 1;
			}
		}
	}
	for(int i = size - 1; i > 0; i--) {
		float temp = arr[0];
		arr[0] = arr[i];
		arr[i] = temp;
		int dad = 0;
		int son = dad * 2 + 1;
		while(son <= i - 1) {
			if(son + 1 <= i - 1 && arr[son] < arr[son + 1]) {
				son++;
			}
			if(arr[dad] > arr[son]) {
				break;
			} else {
				temp = arr[dad];
				arr[dad] = arr[son];
				arr[son] = temp;
				dad = son;
				son = dad * 2 + 1;
			}
		}
	}
}



public static void insertionSort(float[] arr, int size) {
	float temp;
	for(int i = 1; i < size; i++) {
		for(int j = i; j > 0 && arr[j - 1] > arr[j]; j--) {
			temp = arr[j];
			arr[j] = arr[j - 1];
			arr[j - 1] = temp;
		}
	}
}

	

public static int linearSearch(float[] arr, int size, float target) {
	int i = 0;
	int ret = -1;
	while(i <= size - 1) {
		if(arr[i] != target) {
			i++;
		} else {
			ret = i;
			break;
		}
	}
	return ret;
}



public static float meanAbsError(float[] arr, float[] arr1, int size) {
	float tempSum = 0;
	for(int i = 0; i < size; i++) {
		tempSum = (float)(tempSum + Math.abs((double)(arr[i] - arr1[i])));
	}
	tempSum = tempSum / size;
	return tempSum;
}


public static void calCartesianProduct(float[] arr, float[] arr1, int size) {
	for(int i = 0; i < size; i++) {
		for(int j = 0; j < size; j++) {
			System.out.printf("(%f,%f) ", arr[i], arr1[j]);
		}
	}
}


public static void reverse(float[] arr, int size) {
	float temp;
	for(int i = size - 1, j = 0; i > j; i++, j--) {
		temp = arr[i];
		arr[i] = arr[j];
		arr[j] = temp;
	}
}


public static void selectionSort(float[] arr, int size) {
	int min;
	float temp;
	for(int i = 0; i < size - 1; i++) {
		min = i;
		for(int j = i + 1; j < size; j++) {
			if(arr[j] < arr[min]) {
				min = j;
			}
		}
		temp = arr[i];
		arr[i] = arr[min];
		arr[min] = temp;
	}
}





public static float selectK(float[] arr, int size, int k) {
	int j;
	float temp;
	for(int gap = size >> 1; gap > 0; gap >>= 1) {
		for(int i = gap; i < size; i++) {
			temp = arr[i];
			for(j = i - gap; j >= 0 && arr[j] > temp; j -= gap) {
				arr[j + gap] = arr[j];
			}
			arr[j + gap] = temp;
		}
	}
	return arr[size - k];
}



public static void setMinValue(float[] arr, int size, int k) {
	for(int i = 0; i < size; i++) {
		if(arr[i] < k) {
			arr[i] = k;
		}
	}
}



public static void shellSort(float[] arr, int size) {
	int j;
	float temp;
	for(int gap = size >> 1; gap > 0; gap >>= 1) {
		for(int i = gap; i < size; i++) {
			temp = arr[i];
			for(j = i - gap; j >= 0 && arr[j] > temp; j -= gap) {
				arr[j + gap] = arr[j];
			}
			arr[j + gap] = temp;
		}
	}
}





