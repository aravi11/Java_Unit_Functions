static void displayUnionAndIntersection(int[] arrayOne,int[] arrayTwo){
        Set<Integer> obj = new HashSet<>();
        int i,j;
        for(i=0; i<arrayOne.length; i++){
            obj.add(arrayOne[i]);
        }
        for(j=0; j<arrayTwo.length; j++){
            obj.add(arrayTwo[j]);
        }
        System.out.println("The union of both the arrays is");
        for(Integer I: obj){
            System.out.print(I + " ");
        }
        System.out.println();
        obj.clear();
        System.out.println("The intersection of both the arrays is");
        for(i=0; i<arrayOne.length; i++){
            obj.add(arrayOne[i]);
        }
        for(j=0; j<arrayTwo.length; j++){
            if(obj.contains(arrayTwo[j]))
                System.out.print(arrayTwo[j] + " ");
        }
    }


public static int findMin(int[] list)
	{	assert list != null && list.length > 0 : "failed precondition";

		int indexOfMin = 0;
		for(int i = 1; i < list.length; i++)
		{	if(list[i] < list[indexOfMin])
			{	indexOfMin = i;
			}
		}

		return indexOfMin;
	}

public static void badResize(int[] list, int newSize)
	{	assert list != null && newSize >= 0 : "failed precondition";

		int[] temp = new int[newSize];
		int limit = Math.min(list.length, newSize);

		for(int i = 0; i < limit; i++)
		{	temp[i] = list[i];
		}

		// uh oh!! Changing pointer, not pointee. This breaks the
		// relationship between the parameter and argument
		list = temp;
	}

public static void findAndPrintPairs(int[] list, int target)
	{	assert list != null : "failed precondition";

		for(int i = 0; i < list.length; i++)
		{	for(int j = i + 1; j < list.length; j++)
			{	if(list[i] + list[j] == target)
				{	System.out.println("The two elements at indices " + i + " and " + j
						+ " are " + list[i] + " and " + list[j] + " add up to " + target);
				}
			}
		}


	}

public static boolean isAscending( int[] list )
	{	boolean ascending = true;
		int index = 1;
		while( ascending && index < list.length )
		{	assert index >= 0 && index < list.length;

			ascending = (list[index - 1] <= list[index]);
			index++;
		}

		return ascending;
	}

public static void addBonus(double[] array, double bonus)
{
    for (int k = 0; k < array.length; k++)
   {
        array[k] = array[k] + bonus;
   }
}



public static string GetSortedName(String[] names) 
    {
        int n;
        String temp;
        
        for (int i = 0; i < n; i++) 
        {
            for (int j = i + 1; j < n; j++) 
            {
                if (names[i].compareTo(names[j])>0) 
                {
                    temp = names[i];
                    names[i] = names[j];
                    names[j] = temp;
                }


            }
        }
        System.out.print("Names in Sorted Order:");
        for (int i = 0; i < n - 1; i++) 
        {
            System.out.print(names[i] + ",");
        }
        System.out.print(names[n - 1]);
    }
}


 public static void SeperateArrayWithPosition(int a[], int loc) 
    {
        int n, x, flag = 1, loc = 0, k = 0,j = 0;
     
    
        int b[] = new int[n];
        int c[] = new int[n];
      
        for(int i = 0; i < loc; i++)
        {
            b[k] = a[i];
            k++;
        }
        for(int i = loc; i < n; i++)
        {
            c[j] = a[i];
            j++;
        }
        System.out.print("First array:");
        for(int i = 0;i < k; i++)
        {
            System.out.print(b[i]+" ");
        }
        System.out.println("");
        System.out.print("Second array:");
        for(int i = 0; i < j; i++)
        {
            System.out.print(c[i]+" ");
        }
    }
}


 static void segregate0and1(int array[], int size)
    {
        int left = 0, right = size-1;
    	while (left < right)
    	{
            /* Increment left index while we see 0 at left */
            while (array[left] == 0 && left < right)
            left++;
            /* Decrement right index while we see 1 at right */
            while (array[right] == 1 && left < right)
            right--;
            /* If left is smaller than right then there is a 1 at left and a 0 at right.  Exchange it */
            if (left < right)
            {
                array[left] = 0;
                array[right] = 1;
                left++;
                right--;
            }
        }



	for(int i = 0 ; i < 6; i++)
        {
            System.out.print(array[i]+"\t");
        }
    }


public static void CountOccurrence(int a[], int x) 
    {
        int n, x, count = 0, i = 0;
        Scanner s = new Scanner(System.in);
        System.out.print("Enter no. of elements you want in array:");
        n = s.nextInt();
        int a[] = new int[n];
      
        
        for(i = 0; i < n; i++)
        {
            if(a[i] == x)
            {
                count++;
            }
        }
        System.out.println("Number of Occurrence of the Element:"+count);
    }



    static int[] fisherYatesShuffling(int []arr, int n)
    {
        int []a = new int[n];
        int []ind = new int[n];
        for(int i=0; i<n; i++)
            ind[i] = 0;
        int index;
        Random rand = new Random();
        for(int i=0; i<n; i++)
        {
            do
            {
                index = rand.nextInt(n);
            } while(ind[index] != 0);
 
            ind[index] = 1;
            a[i] = arr[index];
        }
        return a;
    }



public void FermatFactor(long N)
    {
        long a = (long) Math.ceil(Math.sqrt(N));
        long b2 = a * a - N;
        while (!isSquare(b2))
        {
            a++;
            b2 = a * a - N;
        }
        long r1 = a - (long)Math.sqrt(b2);
        long r2 = N / r1;
        display(r1, r2);
    }
    /** func


public int[] RotatedArray(int[] A, int K) {
    int [] rotatedA = new int[A.length];
    
    for(int i=0; i<A.length; i++) {
      //rotated index needs to "wrap" around end of array
      int rotatedIndex = (i + K) % A.length;

      rotatedA[rotatedIndex] = A[i];
    }
    return rotatedA;
  }

public int MissingElement(int[] A) {
    int max = A.length + 1;
    
    //load array elements into array so would be quick to check if elements exist
    Set incompleteSet = new HashSet();
    for(int i=0; i<A.length; i++) {
      incompleteSet.add(A[i]);
    }

    for(int i=1; i<max+1; i++) {
      if(!incompleteSet.contains(i)) {
        return (i);
      }
    }
    throw new RuntimeException("shouldn't reach here");
  }

public int MinSumOfTwoSubArray(int[] A) {
    long sumAllElements = 0;
    for(int i=0; i<A.length; i++) {
      sumAllElements += A[i];
    }
    
    int minDifference = Integer.MAX_VALUE;
    int currentDifference = Integer.MAX_VALUE;
    long sumFirstPart = 0;
    long sumSecondPart = 0;

    for(int p=0; p<A.length-1; p++) {
      sumFirstPart += A[p];
      sumSecondPart = sumAllElements - sumFirstPart;
      currentDifference = (int) Math.abs(sumFirstPart - sumSecondPart);
      minDifference = Math.min(currentDifference, minDifference);
    }
    return minDifference;
  }


public int[] MaxCounters(int N, int[] A) {
    int [] counters = new int[N];
    
    int maxCounter = 0; //for the next re-set will know what high score to set all values
    int lastResetCounter = 0; //for setting values that were never explicitly set - at the end
    
    for(int i=0; i<A.length; i++) {
      if(A[i] <= N) {
        if(counters[A[i]-1] < lastResetCounter) { counters[A[i]-1] = lastResetCounter; //bring it up to last reset value } counters[A[i]-1]++; //store max counter in case need to set all counters to this value on next reset if(counters[A[i]-1] > maxCounter) {
          maxCounter = counters[A[i]-1]; 
        }
        
      }
      else {
        //keep track of last reset counter
        lastResetCounter = maxCounter;
      }
    }
    //set all values to last reset value that was never explicitly changed after last reset
    for(int i=0; i<counters.length; i++) {
      if(counters[i] < lastResetCounter) {
        counters[i]  = lastResetCounter;
      }
    }

    return counters;
  }


public int PassingCars(int[] A) {
    int zeros = 0;
    int carPasses = 0;
    
    for(int i=0; i<A.length; i++) { if(A[i] == 0) { zeros++; } else if(A[i] == 1) 
	{ 
	//for every 1 - there will be an extra car pass for ALL the 0's that came before it carPasses += zeros; if(carPasses > 1000000000) {
          return -1;
        }
      }
      else throw new RuntimeException("shouldn't reach here");
    }
    return carPasses;
  }


public int[] GenomeRangeQuery(String S, int[] P, int[] Q) {
    int [] answers = new int[P.length];
    int stringLength = S.length();
    
    int [][] occurrences = new int [stringLength][4];
    
    //step 1 - for each row, count occurrences of each nucleotide (can only have 1 occurrence / row)
    //e.g. if S=CAGCCTA array will be
    //{ 
    //  {0,1,0,0}  C 
    //  {1,0,0,0}  A
    //  {0,0,1,0}  G
    //  {0,1,0,0}  C
    //  {0,1,0,0}  C
    //  {0,0,0,1}  T
      //  {1,0,0,0}  A
    // }
    for(int i=0; i<occurrences.length; i++) {
      char c = S.charAt(i);
      if(c == 'A')      occurrences[i][0] = 1;
      else if(c == 'C') occurrences[i][1] = 1;
      else if(c == 'G') occurrences[i][2] = 1;
      else if(c == 'T') occurrences[i][3] = 1;
    }

    //step 2 - for each row (starting from 2nd row), add up occurrences of each nucleotide by adding
    //occurrences from previous row to current row
    //now have running sum of each nucleotide for any row
    //e.g. if S=CAGCCTA array will be
    //{ 
    //  {0,1,0,0}  C 
    //  {1,1,0,0}  A
    //  {1,1,1,0}  G
    //  {1,2,1,0}  C
    //  {1,3,1,0}  C
    //  {1,3,1,1}  T
      //  {2,3,1,1}  A
    // }
    for(int i=1; i<stringLength; i++) {
      for(int j=0; j<4; j++) {
        occurrences[i][j] += occurrences[i-1][j];
      }
    }

    //check each slice for min by simple subtraction
    for(int i=0; i<P.length; i++) {
      int index1 = P[i];
      int index2 = Q[i];

      for(int j=0; j<4; j++) { int lowerIndexCount = 0; //when index1 not at beginning of String - need to get occurrences from just before //beginning of slice - to see if that nucleotide occurred within slice //e.g. if slice is (2, 4), need to check for occurrences of A, C, G, T from index 1 to 4 if(index1-1 >= 0) { 
          lowerIndexCount = occurrences[index1-1][j];
        }
        
        if(occurrences[index2][j] - lowerIndexCount > 0) {
          answers[i] = j+1; //nucleotide value is 1 more than loop value (A=1, C=2, G=3, T=4)
          //no need to keep checking since always checking from smallest impact factor
          //as soon as find occurrence, that has to be minimum, cause subsequent nucleotides have
          //larger impact factor
          break; 
        }
      }
    }
    return answers;
  }


public int MinAverageTwoSlice(int[] A) {

    //main idea: will find min average by checking only 2 and 3 contiguous elements at a time
    int sum1, sum2 = 0;
    double minAverage = Double.MAX_VALUE;
    double currentAverage1 = Double.MAX_VALUE;
    double currentAverage2 = Double.MAX_VALUE;
    int minAverageSliceIndex = 0; //for size 2 arrays, this will always be true

    //if array is > 2 elements
    for(int i=0; i<A.length-2; i++) {
      sum1 = A[i] + A[i+1];
      currentAverage1 = sum1 / 2.0d;
      if(currentAverage1 < minAverage) {
        minAverage = currentAverage1;
        minAverageSliceIndex = i;
      }

      sum2 = sum1 + A[i+2];
      currentAverage2 = sum2 / 3.0d;
      if(currentAverage2 < minAverage) {
        minAverage = currentAverage2;
        minAverageSliceIndex = i;
      }
    }

    //check last 2 contiguous elements from the end - they won't otherwise be checked because 
    //when checking 2 and 3 contiguous elements at a time, will stop checking 3 elements from the end 
    currentAverage1 = (A[A.length-2] + A[A.length-1]) / 2.0d;
    if(currentAverage1 < minAverage) {
      minAverage = currentAverage1;
      minAverageSliceIndex = A.length-2;
    }
    
    return minAverageSliceIndex;
  }




public int IsTrianglePossible(int[] A) {
      if(A.length < 3) return 0;
      
    List aList = new ArrayList();
    for(int i=0; i<A.length; i++) { aList.add(A[i]); } Collections.sort(aList); //made long array because each int element can be as high as Integer.MAX_VALUE so when add them //can overflow int long [] aOrdered = new long[A.length]; int index = 0; for(Integer i : aList) { aOrdered[index++] = i; } //start from the end (largest) //if previous 2 elements have sum > current element, found a triangle 
    for(int i=aOrdered.length-1; i>=2; i--) {
      if(aOrdered[i-1] + aOrdered[i-2] > aOrdered[i]) {
        return 1;
      }
    }
      return 0;
    }

public int MaximumProfitOfStockTicker(int[] A) {
    if(A.length < 2) return 0; //for empty array or 1 element array, no profit can be realized
    
    //convert profit table to delta table so can use max slice technique
    int [] deltaA = new int[A.length];
    deltaA[0] = 0;
    for(int i=1; i<A.length; i++) {
      deltaA[i] = A[i] - A[i-1];
    }
    
    int absoluteMax = deltaA[0];
    int localMax = deltaA[0];
    int nextSum = 0;
    int currentElement = 0;
    
    for (int i = 1; i < deltaA.length; i++) { 
	currentElement = deltaA[i]; 
	nextSum = localMax + currentElement; 
	localMax = Math.max(deltaA[i], nextSum); 
	
	absoluteMax = Math.max(absoluteMax, localMax); 
   } 
	if(absoluteMax > 0) return absoluteMax;
    
    return 0; 
  }



public int MaxDoubleSliceSum(int[] A) {
    int[] slice1LocalMax = new int[A.length];
    int[] slice2LocalMax = new int[A.length];
    
    //start from i=1 because slice can't start at index 0
    for(int i = 1; i < A.length-1; i++) { slice1LocalMax[i] = Math.max(slice1LocalMax[i-1] + A[i], 0); } //start from i=A.length-2 because slice can't end at index A.length-1 for(int i = A.length-2; i > 0; i--){
      slice2LocalMax[i] = Math.max(slice2LocalMax[i+1]+A[i], 0);
    }
    
    int maxDoubleSliceSum = 0;
    
    //compute sums of all slices to find absolute max
    for(int i = 1; i < A.length-1; i++) {
      maxDoubleSliceSum = Math.max(maxDoubleSliceSum, slice1LocalMax[i-1] + slice2LocalMax[i+1]);
    }
    
    return maxDoubleSliceSum;
  }



public int[] ArePrime(int N) {
    //make size N+1 so will have direct mapping from array index
    boolean [] arePrimes = new boolean[N+1];

    arePrimes[0] = false; //0 is never prime
    arePrimes[1] = false; //1 is never prime
    for(int i=2; i<arePrimes.length; i++) {
      arePrimes[i] = true;
    }

    int nSquareRoot = (int) Math.sqrt(N);
    for(int i=2; i<=nSquareRoot; i++) {
      if(arePrimes[i]) {
        //start checking from i^2 because lower numbers will have already been checked
        //keep checking very multiple of i
        for(int j=i*i; j<=N; j+=i) {
          arePrimes[j] = false;
        }
      }
    }

    List primesList = new ArrayList();

    for(int i=2; i<arePrimes.length; i++) {
      if(arePrimes[i]) {
        primesList.add(i);
      }
    }
    //https://stackoverflow.com/questions/960431/how-to-convert-listinteger-to-int-in-java
    int[] primes = primesList.stream().mapToInt(i->i).toArray();
    return primes;
  }


public static int lonelyinteger(int[] a) {

		int result = 0;

		for (int i = 0; i < a.length; i++) {

			result ^= a[i];

		}

		return result;

}


static int minimumLoss(long[] price) {

		Map<Long,Integer> map=new HashMap<Long,Integer>();

		for(int i=0;i<price.length;i++) {

			map.put(price[i], i);

		}

		Arrays.sort(price);

		long min=Long.MAX_VALUE;

		for(int i=price.length-1;i>0;i--) {

			long diff=price[i]-price[i-1];

			if( diff<min && map.get(price[i])<map.get(price[i-1])) {

				min=Math.abs(diff);

			}

		}

		return (int) min;

    }


public static int ToyContainerRequired(int[] w) {

		Arrays.sort(w);

		int containerCouter = 1, weightLimit = w[0];

		for (int i = 0; i < w.length; i++) {

			if (w[i] > weightLimit + 4) {

				weightLimit = w[i];

				containerCouter++;

			}

		}

		return containerCouter;

	}



public static int shiftCount(int[] a) {

		int count = 0;

		for (int i = 0; i < a.length; i++) {

			for (int j = i; j > 0; j--) {

				if (a[j] < a[j - 1]) {

					int temp = a[j - 1];

					a[j - 1] = a[j];

					a[j] = temp;

					count++;

				}

			}

		}



		return count;

	}




static int birthdayCakeCandles(int[] ar) {

		int maxCandleHeight = Integer.MIN_VALUE;

		int maxCandleFreqCount = 0;



		for (int i = 0; i < ar.length; i++) {



			if (ar[i] == maxCandleHeight) {

				maxCandleFreqCount++;

			}



			if (ar[i] > maxCandleHeight) {

				maxCandleHeight = ar[i];

				maxCandleFreqCount = 1;

			}



		}

		return maxCandleFreqCount;



	}


static int[] matchingStrings(String[] strings, String[] queries) {

		Map<String, Integer> map = new HashMap<>();

		int result[] = new int[queries.length];



		for (int i = 0; i < strings.length; i++) {

			String inputString = strings[i];

			if (map.containsKey(inputString)) {

				map.put(inputString, map.get(inputString) + 1);

			} else {

				map.put(inputString, 1);

			}

		}



		for (int i = 0; i < queries.length; i++) {

			String queryString = queries[i];

			if (map.containsKey(queryString)) {

				result[i] = map.get(queryString);

			}

		}



	static int divisibleSumPairs(int n, int k, int[] ar) {

		int[] bucket = new int[k];

		int pairCounter = 0;



		for (int element : ar) {



			int remainder = element % k;

			int complement = (k-remainder)%k;

			pairCounter += ar[complement];

			bucket[remainder]++;



		}

		return pairCounter;

	}



   public static int[] replaceMiddle(int[] arr) {

        if (arr.length % 2 == 0) {

            int[] sol = new int[arr.length - 1];

            for (int i = 0, j = 0; i < arr.length; j++) {

                if (i == arr.length / 2 - 1) {

                    sol[j] = arr[arr.length / 2] + arr[arr.length / 2 - 1];

                    i += 2;

                } else {

                    sol[j] = arr[i];

                    i++;

                }

            }

            return sol;



        }

        return arr;



   public static int largestDistance(int[] a) {



        int[] mx = new int[] { a[0], a[1] };

        int[] mn = new int[] { a[a.length - 2], a[a.length - 1] };

        for (int i = 0; i < a.length; i++) {

            int k = i % 2;

            if (a[i] > mx[k]) {

                mx[k] = a[i];

            } else if (a[i] < mn[k]) {

                mn[k] = a[i];

            }

        }

        return Math.max(mx[0] - mn[0], mx[1] - mn[1]);

    }


  private static double findMedian(int [] array, int start, int end) {

        if ((end - start) % 2 == 0) { // odd number of elements

            return (array[(end + start) / 2]);

        } else { // even number of elements

            int value1 = array[(end + start) / 2];

            int value2 = array[(end + start) / 2 + 1];

            return (value1 + value2) / 2.0;

        }

    }


  static int[] reverseArray(int[] a) {

        int[] reverseArr = new int[a.length];



        for (int i = 0; i < a.length; i++) {

            reverseArr[i] = a[a.length-1-i];

        }



        return reverseArr;

    }

    }

	return result;

}

public int rob(int[] nums) {

        if(nums.length == 0) {

            return 0;

        }



        if(nums.length == 1) {

            return nums[0];

        }

        

        int[] dp = new int[nums.length];

        

        dp[0] = nums[0];

        dp[1] = nums[0] > nums[1] ? nums[0] : nums[1];



        for(int i = 2; i < nums.length; i++) {

            dp[i] = Math.max(dp[i - 2] + nums[i], dp[i - 1]);

        }

        

        return dp[dp.length - 1];

    }


public boolean wordBreak(String s, Set<String> wordDict) {

        boolean[] dp = new boolean[s.length() + 1];

        

        dp[0] = true;

        

        for(int i = 1; i <= s.length(); i++) {

            for(int j = 0; j < i; j++) {

                if(dp[j] && wordDict.contains(s.substring(j, i))) {

                    dp[i] = true;

                    break;

                }

            }

        }
        return dp[s.length()];
    }




public int[] productExceptSelf(int[] nums) {

        int n = nums.length;

        int[] result = new int[n];

        int left = 1;

        for(int i = 0; i < nums.length; i++) {

            if(i > 0) {

                left *= nums[i - 1];
            }

            result[i] = left;
        }

        int right = 1;

        for(int i = n - 1; i >= 0; i--) {

            if(i < n - 1) {

                right *= nums[i + 1];

            }
            result[i] *= right;
        }
        return result;
    }
}





    public int TrappingRainWater(int[] height) {

        int water = 0;

        

        int leftIndex = 0;

        int rightIndex = height.length - 1;       

        int leftMax = 0;

        int rightMax = 0;

        while(leftIndex <= rightIndex) {

            leftMax = Math.max(leftMax, height[leftIndex]);

            rightMax = Math.max(rightMax, height[rightIndex]);
  
            if(leftMax < rightMax) {

                water += leftMax - height[leftIndex];

                leftIndex++;

            } else {

                water += rightMax - height[rightIndex];
                rightIndex--;
            }
        }
        return water;

    }



public int minCostClimbingStairs(int[] cost) {

        if(cost == null || cost.length == 0) {

            return 0;

        }

        if(cost.length == 1) {

            return cost[0];

        }

        if(cost.length == 2) {

            return Math.min(cost[0], cost[1]);

        }
        int[] dp = new int[cost.length];

        dp[0] = cost[0];

        dp[1] = cost[1];

        for(int i = 2; i < cost.length; i++) {

            dp[i] = Math.min(dp[i - 1] + cost[i], dp[i - 2] + cost[i]);

        }
        return Math.min(dp[cost.length - 1], dp[cost.length -2]);

    }




 public int combinationSum4(int[] nums, int target) {

        int[] dp = new int[target + 1];

        dp[0] = 1;

        for(int i = 1; i < dp.length; i++) {

            for(int j = 0; j < nums.length; j++) {                

                if(i - nums[j] >= 0) {

                    dp[i] += dp[i - nums[j]];
                }
               
            }
        }

       return dp[target];     
    }



 public void sortColors(int[] nums) {

        int wall = 0;

        

        for(int i = 0; i < nums.length; i++) {

            if(nums[i] < 1) {

                int temp = nums[i];

                nums[i] = nums[wall];

                nums[wall] = temp;

                wall++;
            }
        }

        for(int i = 0; i < nums.length; i++) {

            if(nums[i] == 1) {

                int temp = nums[i];

                nums[i] = nums[wall];

                nums[wall] = temp;

                wall++;

            }

        }

    }



    public int minSubArrayLen(int s, int[] nums) {

        if(nums == null || nums.length == 0) {

            return 0;

        }

        

        int i = 0;

        int j = 0;

        int result = Integer.MAX_VALUE;

        int total = 0;

        while(i < nums.length) {

            total += nums[i++];

            while(total >= s) {

                result = Math.min(result, i - j);

                total -= nums[j++];
            }
        }

        return result == Integer.MAX_VALUE ? 0 : result;

    }



public void moveZeroes(int[] nums) {

        if(nums == null || nums.length == 0) {

            return;

        }

        

        int index = 0;

        for(int num : nums) {

            if(num != 0) {

                nums[index] = num;

                index++;

            }
        }
        while(index < nums.length) {

            nums[index] = 0;

            index++;
        }

    }