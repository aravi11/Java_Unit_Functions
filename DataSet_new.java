public static Set<String> getPermutation(String str) {
    // create a set to avoid duplicate permutation
    Set<String> permutations = new HashSet<String>();

    // check if string is null
    if (str == null) {
      return null;
    } else if (str.length() == 0) {
      // terminating condition for recursion
      permutations.add("");
      return permutations;
    }

    // get the first character
    char first = str.charAt(0);

    // get the remaining substring
    String sub = str.substring(1);

    // make recursive call to getPermutation()
    Set<String> words = getPermutation(sub);

    // access each element from words
    for (String strNew : words) {
      for (int i = 0;i<=strNew.length();i++){

        // insert the permutation to the set
        permutations.add(strNew.substring(0, i) + first + strNew.substring(i));
      }
    }
    return permutations;
  }
  
  

public static int reverseNumber(int num) {

        int reversed = 0;

        while(num != 0) {
            int digit = num % 10;
            reversed = reversed * 10 + digit;
            num /= 10;
        }

        System.out.println("Reversed Number: " + reversed);
    }
	
	
public static bool isPalindrom(int num) {

        int reversedInteger = 0, remainder, originalInteger;

        originalInteger = num;

        // reversed integer is stored in variable 
        while( num != 0 )
        {
            remainder = num % 10;
            reversedInteger = reversedInteger * 10 + remainder;
            num  /= 10;
        }

        // palindrome if orignalInteger and reversedInteger are equal
        if (originalInteger == reversedInteger)
            return true;
        else
            return false;
    }
	
	

public static bool isArmstrongNumber(int number) {

        int originalNumber, remainder, result = 0;

        originalNumber = number;

        while (originalNumber != 0)
        {
            remainder = originalNumber % 10;
            result += Math.pow(remainder, 3);
            originalNumber /= 10;
        }

        if(result == number)
            return true;
        else
            return false;
    }
	

public static void printFactors(int number) {
        System.out.print("Factors of " + number + " are: ");
        for(int i = 1; i <= number; ++i) {
            if (number % i == 0) {
                System.out.print(i + " ");
            }
        }
    }
	
	

public static int addNumbers(int num) {
        if (num != 0)
            return num + addNumbers(num - 1);
        else
            return num;
    }
	
	

public static long multiplyNumbers(int num)
    {
        if (num >= 1)
            return num * multiplyNumbers(num - 1);
        else
            return 1;
    }
	


 public  static boolean checkPrime(int num) {
        boolean isPrime = true;

        for (int i = 2; i <= num / 2; ++i) {
            if (num % i == 0) {
                isPrime = false;
                break;
            }
        }

        return isPrime;
    }
	

public static double calculateSD(double numArray[])
    {
        double sum = 0.0, standardDeviation = 0.0;
        int length = numArray.length;

        for(double num : numArray) {
            sum += num;
        }

        double mean = sum/length;

        for(double num: numArray) {
            standardDeviation += Math.pow(num - mean, 2);
        }

        return Math.sqrt(standardDeviation/length);
    }
	
	
public static int[][] addTwoMatrix(int[][] firstMatrix, int[][] secondMatrix) {
        int rows = 2, columns = 3;


        // Adding Two matrices
        int[][] sum = new int[rows][columns];
        for(int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                sum[i][j] = firstMatrix[i][j] + secondMatrix[i][j];
            }
        }

        // Displaying the result
        System.out.println("Sum of two matrices is: ");
        for(int[] row : sum) {
            for (int column : row) {
                System.out.print(column + "    ");
            }
            System.out.println();
        }
    }
	

public static int characterFrequency(String str, char ch ) {
        int frequency = 0;

        for(int i = 0; i < str.length(); i++) {
            if(ch == str.charAt(i)) {
                ++frequency;
            }
        }

        return frequency;
    }
	
	
public static void countVowelAndConsonant(String line) {

	int vowels = 0, consonants = 0, digits = 0, spaces = 0;

	line = line.toLowerCase();
	for(int i = 0; i < line.length(); ++i)
	{
		char ch = line.charAt(i);
		if(ch == 'a' || ch == 'e' || ch == 'i'
			|| ch == 'o' || ch == 'u') {
			++vowels;
		}
		else if((ch >= 'a'&& ch <= 'z')) {
			++consonants;
		}
		else if( ch >= '0' && ch <= '9')
		{
			++digits;
		}
		else if (ch ==' ')
		{
			++spaces;
		}
	}

	System.out.println("Vowels: " + vowels);
	System.out.println("Consonants: " + consonants);
	System.out.println("Digits: " + digits);
	System.out.println("White spaces: " + spaces);
}


public static void sortWords(String[] words) {

	for(int i = 0; i < 3; ++i) {
		for (int j = i + 1; j < 4; ++j) {
			if (words[i].compareTo(words[j]) > 0) {

				// swap words[i] with words[j[
				String temp = words[i];
				words[i] = words[j];
				words[j] = temp;
			}
		}
	}

	System.out.println("In lexicographical order:");
	for(int i = 0; i < 4; i++) {
		System.out.println(words[i]);
	}
}


public static Time difference(Time start, Time stop)
{
	Time diff = new Time(0, 0, 0);

	if(start.seconds > stop.seconds){
		--stop.minutes;
		stop.seconds += 60;
	}

	diff.seconds = stop.seconds - start.seconds;

	if(start.minutes > stop.minutes){
		--stop.hours;
		stop.minutes += 60;
	}

	diff.minutes = stop.minutes - start.minutes;
	diff.hours = stop.hours - start.hours;

	// return the difference time
	return(diff);
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


public static int[] arrayOfMultiples(int num, int length) {
	int [] arr = new int [length];
	for (int i = 0; i < length; i++) {
		arr[i] = num * (i+1);
	}
return arr;
}


public static int[][] squarePatch(int n) {
		int[][] finalArray = new int[n][n];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				finalArray[i][j] = n;
			}
		}
		return finalArray;
	}
	
public static int solutions(int a, int b, int c) {
	int discriminant = (b * b) - (4 * a * c);
	if(discriminant < 0 ) return 0;
	return discriminant == 0 ? 1 : 2;
}


public static boolean isPandigital(long num) {
	Set<Long> newSet = new HashSet<>();
	while(num != 0){
		newSet.add(num % 10);
		num /= 10;
	}
	return newSet.size() == 10;
}


public static boolean rectangleInCircle(int w, int h, int radius) {
		int num=h*h+w*w;
		int dia=2*radius;
		int squ=dia*dia;
		if(num<=squ)
			return true;
		return false;
  }
  

public static String moveToNextString(String word) {
	  	String result = "";
	  
	  	for (int i = 0; i < word.length(); i++) {
		  	char newLetter = word.charAt(i);
		  	newLetter += 1;
		  	result = result + newLetter;
		}
	  
	  	return result;
	}


public static boolean isAutomorphic(int n) {
		int x = n;
		int i = 0;
		while(x!=0){
			x/=10;
			i++;
		}
		return (Math.pow(n,2)%Math.pow(10,i)==n);
  }
  
public static int minSwaps(String s1, String s2) {
    int diff =0;
		for (int i=0;i<s1.length();i++) {
				if (s1.charAt(i) != s2.charAt(i)) diff++;
		}
		return diff/2;				
  }
  

public <T extends Comparable<T>> int ternarySearch(T[] arr, T key, int start, int end) {
        if (start > end) {
            return -1;
        }
        /* First boundary: add 1/3 of length to start */
        int mid1 = start + (end - start) / 3;
        /* Second boundary: add 2/3 of length to start */
        int mid2 = start + 2 * (end - start) / 3;

        if (key.compareTo(arr[mid1]) == 0) {
            return mid1;
        } else if (key.compareTo(arr[mid2]) == 0) {
            return mid2;
        }

        /* Search the first (1/3) rd part of the array.*/

        else if (key.compareTo(arr[mid1]) < 0) {
            return ternarySearch(arr, key, start, --mid1);
        }
        /* Search 3rd (1/3)rd part of the array */

        else if (key.compareTo(arr[mid2]) > 0) {
            return ternarySearch(arr, key, ++mid2, end);
        }
        /* Search middle (1/3)rd part of the array */

        else {
            return ternarySearch(arr, key, mid1, mid2);
        }
    }
	
	
public static int[] find(int arr[][], int row, int col, int key) {

	//array to store the answer row and column
	int ans[] = {-1, -1};
	if (row < 0 || col >= arr[row].length) {
		return ans;
	}
	if (arr[row][col] == key) {
		ans[0] = row;
		ans[1] = col;
		return ans;
	}
	//if the current element is greater than the given element then we move up
	else if (arr[row][col] > key) {
		return find(arr, row - 1, col, key);
	}
	//else we move right
	return find(arr, row, col + 1, key);
}

public int perfectBinarySearch(int[] arr, int target) 
    {
        int low = 0 ;
        int high = arr.length - 1 ;

        while(low <= high) {
            int mid =(low + high) / 2;

            if(arr[mid] == target) {
                return mid;
            }
            else if(arr[mid] > target) {
                high = mid - 1;
            }
            else {
                low = mid + 1;
            }

        }
        return -1;
    }
	
public static boolean validForBase(String n, int base) {
	char[] validDigits = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E',
			'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V',
			'W', 'X', 'Y', 'Z'};
	// digitsForBase contains all the valid digits for the base given
	char[] digitsForBase = Arrays.copyOfRange(validDigits, 0, base);

	// Convert character array into set for convenience of contains() method
	HashSet<Character> digitsList = new HashSet<>();
	for (int i = 0; i < digitsForBase.length; i++)
		digitsList.add(digitsForBase[i]);

	// Check that every digit in n is within the list of valid digits for that base.
	for (char c : n.toCharArray())
		if (!digitsList.contains(c))
			return false;

	return true;
}


public static String base2base(String n, int b1, int b2) {
        int decimalValue = 0, charB2;
        char charB1;
        String output = "";
        // Go through every character of n
        for (int i = 0; i < n.length(); i++) {
            // store the character in charB1
            charB1 = n.charAt(i);
            // if it is a non-number, convert it to a decimal value >9 and store it in charB2
            if (charB1 >= 'A' && charB1 <= 'Z')
                charB2 = 10 + (charB1 - 'A');
                // Else, store the integer value in charB2
            else
                charB2 = charB1 - '0';
            decimalValue = decimalValue * b1 + charB2;
        }

        while (decimalValue != 0) {
            // If the remainder is a digit < 10, simply add it to
            // the left side of the new number.
            if (decimalValue % b2 < 10)
                output = Integer.toString(decimalValue % b2) + output;
                // If the remainder is >= 10, add a character with the
                // corresponding value to the new number. (A = 10, B = 11, C = 12, ...)
            else
                output = (char) ((decimalValue % b2) + 55) + output;
            // Divide by the new base again
            decimalValue /= b2;
        }
        return output;
    }
	
public static int convertToDecimal(String s, int radix) {
        int num = 0;
        int pow = 1;

        for (int i = s.length() - 1; i >= 0; i--) {
            int digit = valOfChar(s.charAt(i));
            if (digit >= radix) {
                throw new NumberFormatException("For input string " + s);
            }
            num += valOfChar(s.charAt(i)) * pow;
            pow *= radix;
        }
        return num;
    }
	
static String binToHex(int binary) {
        //hm to store hexadecimal codes for binary numbers within the range: 0000 to 1111 i.e. for decimal numbers 0 to 15
        HashMap<Integer, String> hm = new HashMap<>();
        //String to store hexadecimal code
        String hex = "";
        int i;
        for (i = 0; i < 10; i++) {
            hm.put(i, String.valueOf(i));
        }
        for (i = 10; i < 16; i++) hm.put(i, String.valueOf((char) ('A' + i - 10)));
        int currbit;
        while (binary != 0) {
            int code4 = 0;    //to store decimal equivalent of number formed by 4 decimal digits
            for (i = 0; i < 4; i++) {
                currbit = binary % 10;
                binary = binary / 10;
                code4 += currbit * Math.pow(2, i);
            }
            hex = hm.get(code4) + hex;
        }
        return hex;
    }
	
public static String convertBinaryToOctal(int binary) {
        String octal = "";
        int currBit = 0, j = 1;
        while (binary != 0) {
            int code3 = 0;
            for (int i = 0; i < 3; i++) {
                currBit = binary % 10;
                binary = binary / 10;
                code3 += currBit * j;
                j *= 2;
            }
            octal = code3 + octal;
            j = 1;
        }
        return octal;
    }
	
public static String convertToAnyBase(int inp, int base) {
	ArrayList<Character> charArr = new ArrayList<>();

	while (inp > 0) {
		charArr.add(reVal(inp%base));
		inp /= base;
	}

	StringBuilder str = new StringBuilder(charArr.size());

	for(Character ch: charArr)
	{
		str.append(ch);
	}

	return str.reverse().toString();
}

public static String decToHex(int dec) {
	StringBuilder hexBuilder = new StringBuilder(sizeOfIntInHalfBytes);
	hexBuilder.setLength(sizeOfIntInHalfBytes);
	for (int i = sizeOfIntInHalfBytes - 1; i >= 0; --i) {
		int j = dec & halfByte;
		hexBuilder.setCharAt(i, hexDigits[j]);
		dec >>= numberOfBitsInAHalfByte;
	}
	return hexBuilder.toString().toLowerCase();
}

public static int decimal2octal(int q) {
	int now;
	int i = 1;
	int octnum = 0;
	while (q > 0) {
		now = q % 8;
		octnum = (now * (int) (Math.pow(10, i))) + octnum;
		q /= 8;
		i++;
	}
	octnum /= 10;
	return octnum;
}

public static String integerToRoman(int num) {
	int[] allArabianRomanNumbers = new int[]{1000, 900, 500, 400, 100, 90, 50, 40, 10, 9, 5, 4, 1};
	String[] allRomanNumbers = new String[]{"M", "CM", "D", "CD", "C", "XC", "L", "XL", "X", "IX", "V", "IV", "I"};
	if (num <= 0) {
		return "";
	}

	StringBuilder builder = new StringBuilder();

	for (int a = 0; a < allArabianRomanNumbers.length; a++) {
		int times = num / allArabianRomanNumbers[a];
		for (int b = 0; b < times; b++) {
			builder.append(allRomanNumbers[a]);
		}

		num -= times * allArabianRomanNumbers[a];
	}

	return builder.toString();
}


public static int romanToInt(String A) {

	A = A.toUpperCase();
	char prev = ' ';

	int sum = 0;

	int newPrev = 0;
	for (int i = A.length() - 1; i >= 0; i--) {
		char c = A.charAt(i);

		if (prev != ' ') {
			// checking current Number greater then previous or not
			newPrev = map.get(prev) > newPrev ? map.get(prev) : newPrev;
		}

		int currentNum = map.get(c);

		// if current number greater then prev max previous then add
		if (currentNum >= newPrev) {
			sum += currentNum;
		} else {
			// subtract upcoming number until upcoming number not greater then prev max
			sum -= currentNum;
		}

		prev = c;
	}

	return sum;
}

public static int bpR(int start,int end){
	if(start==end) {
		return 1;
	}
	else if(start>end)
		return 0;
	int count=0;
	for(int dice=1;dice<=6;dice++) {
		count+=bpR(start+dice,end);
	}
	return count;
}

public static int bpRS(int curr,int end,int strg[]){
	if(curr==end) {
		return 1;
	}
	else if(curr>end)
		return 0;
	if(strg[curr]!=0)
		return strg[curr];
	int count=0;
	for(int dice=1;dice<=6;dice++) {
		count+=bpRS(curr+dice,end,strg);
	}
	strg[curr]=count;
	return count;
}

public static int bpIS(int curr,int end,int[]strg){
	strg[end]=1;
	for(int i=end-1;i>=0;i--) {
		int count=0;
		for(int dice=1;dice<=6&&dice+i<strg.length;dice++) {
			count+=strg[i+dice];
		}
		strg[i]=count;
	}
	return strg[0];
}


 public static int change(int[] coins, int amount) {

	int[] combinations = new int[amount + 1];
	combinations[0] = 1;

	for (int coin : coins) {
		for (int i = coin; i < amount + 1; i++) {
			combinations[i] += combinations[i - coin];
		}
		// Uncomment the below line to see the state of combinations for each coin
		// printAmount(combinations);
	}

	return combinations[amount];
}

public static int minimumCoins(int[] coins, int amount) {
	//minimumCoins[i] will store the minimum coins needed for amount i
	int[] minimumCoins = new int[amount + 1];

	minimumCoins[0] = 0;

	for (int i = 1; i <= amount; i++) {
		minimumCoins[i] = Integer.MAX_VALUE;
	}
	for (int i = 1; i <= amount; i++) {
		for (int coin : coins) {
			if (coin <= i) {
				int sub_res = minimumCoins[i - coin];
				if (sub_res != Integer.MAX_VALUE && sub_res + 1 < minimumCoins[i])
					minimumCoins[i] = sub_res + 1;
			}
		}
	}
	// Uncomment the below line to see the state of combinations for each coin
	//printAmount(minimumCoins);
	return minimumCoins[amount];
}


public static int numStrIS(int n) {
	int[] zeros=new int[n];
	int []ones=new int[n];
	//seed
	zeros[0]=1;
	ones[0]=1;
	for(int i=1;i<n;i++) {
		zeros[i]=zeros[i-1]+ones[i-1];
		ones[i]=zeros[i-1];
	}
	int ans=zeros[n-1]+ones[n-1];
	return ans;
}

public static int countStrings(int n, int lastDigit)
{
	if (n == 0) {
		return 0;
	}

	// if only one digit is left
	if (n == 1) {
		return (lastDigit == 1) ? 1: 2;
	}

	// if last digit is 0, we can have both 0 and 1 at current pos
	if (lastDigit == 0) {
		return countStrings(n - 1, 0) + countStrings(n - 1, 1);
	}
	// if last digit is 1, we can have only 0 at current position
	else {
		return countStrings(n - 1, 0);
	}
}

public static int minDistance(String word1, String word2) {
	int len1 = word1.length();
	int len2 = word2.length();
	// len1+1, len2+1, because finally return dp[len1][len2]
	int[][] dp = new int[len1 + 1][len2 + 1];
	/* If second string is empty, the only option is to
  insert all characters of first string into second*/
	for (int i = 0; i <= len1; i++) {
		dp[i][0] = i;
	}
	/* If first string is empty, the only option is to
  insert all characters of second string into first*/
	for (int j = 0; j <= len2; j++) {
		dp[0][j] = j;
	}
	//iterate though, and check last char
	for (int i = 0; i < len1; i++) {
		char c1 = word1.charAt(i);
		for (int j = 0; j < len2; j++) {
			char c2 = word2.charAt(j);
			//if last two chars equal
			if (c1 == c2) {
				//update dp value for +1 length
				dp[i + 1][j + 1] = dp[i][j];
			} else {
		/* if two characters are different ,
		then take the minimum of the various operations(i.e insertion,removal,substitution)*/
				int replace = dp[i][j] + 1;
				int insert = dp[i][j + 1] + 1;
				int delete = dp[i + 1][j] + 1;

				int min = replace > insert ? insert : replace;
				min = delete > min ? min : delete;
				dp[i + 1][j + 1] = min;
			}
		}
	}
	/* return the final answer , after traversing through both the strings*/
	return dp[len1][len2];
}
	
	
public static int minTrials(int n, int m) {

	int[][] eggFloor = new int[n + 1][m + 1];
	int result, x;

	for (int i = 1; i <= n; i++) {
		eggFloor[i][0] = 0;   // Zero trial for zero floor.
		eggFloor[i][1] = 1;   // One trial for one floor 
	}

	// j trials for only 1 egg

	for (int j = 1; j <= m; j++)
		eggFloor[1][j] = j;

	// Using bottom-up approach in DP

	for (int i = 2; i <= n; i++) {
		for (int j = 2; j <= m; j++) {
			eggFloor[i][j] = Integer.MAX_VALUE;
			for (x = 1; x <= j; x++) {
				result = 1 + Math.max(eggFloor[i - 1][x - 1], eggFloor[i][j - x]);

				// choose min of all values for particular x
				if (result < eggFloor[i][j])
					eggFloor[i][j] = result;
			}
		}
	}

	return eggFloor[n][m];
}


public static int networkFlow(int source, int sink) {
	flow = new int[V][V];
	int totalFlow = 0;
	while (true) {
		Vector<Integer> parent = new Vector<>(V);
		for (int i = 0; i < V; i++)
			parent.add(-1);
		Queue<Integer> q = new LinkedList<>();
		parent.set(source, source);
		q.add(source);
		while (!q.isEmpty() && parent.get(sink) == -1) {
			int here = q.peek();
			q.poll();
			for (int there = 0; there < V; ++there)
				if (capacity[here][there] - flow[here][there] > 0 && parent.get(there) == -1) {
					q.add(there);
					parent.set(there, here);
				}
		}
		if (parent.get(sink) == -1)
			break;

		int amount = INF;
		String printer = "path : ";
		StringBuilder sb = new StringBuilder();
		for (int p = sink; p != source; p = parent.get(p)) {
			amount = Math.min(capacity[parent.get(p)][p] - flow[parent.get(p)][p], amount);
			sb.append(p + "-");
		}
		sb.append(source);
		for (int p = sink; p != source; p = parent.get(p)) {
			flow[parent.get(p)][p] += amount;
			flow[p][parent.get(p)] -= amount;
		}
		totalFlow += amount;
		printer += sb.reverse() + " / max flow : " + totalFlow;
		System.out.println(printer);
	}

	return totalFlow;
}


public static int largestContiguousSum(int arr[]) {
	int i, len = arr.length, cursum = 0, maxsum = Integer.MIN_VALUE;
	if (len == 0)    //empty array
		return 0;
	for (i = 0; i < len; i++) {
		cursum += arr[i];
		if (cursum > maxsum) {
			maxsum = cursum;
		}
		if (cursum <= 0) {
			cursum = 0;
		}
	}
	return maxsum;
}


public static int knapSack(int W, int wt[], int val[], int n) throws IllegalArgumentException {
	if(wt == null || val == null)
		throw new IllegalArgumentException();
	int i, w;
	int rv[][] = new int[n + 1][W + 1];    //rv means return value

	// Build table rv[][] in bottom up manner
	for (i = 0; i <= n; i++) {
		for (w = 0; w <= W; w++) {
			if (i == 0 || w == 0)
				rv[i][w] = 0;
			else if (wt[i - 1] <= w)
				rv[i][w] = Math.max(val[i - 1] + rv[i - 1][w - wt[i - 1]], rv[i - 1][w]);
			else
				rv[i][w] = rv[i - 1][w];
		}
	}

	return rv[n][W];
}


public static int calculate_distance(String a, String b) {
	int len_a = a.length() + 1;
	int len_b = b.length() + 1;
	int[][] distance_mat = new int[len_a][len_b];
	for (int i = 0; i < len_a; i++) {
		distance_mat[i][0] = i;
	}
	for (int j = 0; j < len_b; j++) {
		distance_mat[0][j] = j;
	}
	for (int i = 0; i < len_a; i++) {
		for (int j = 0; j < len_b; j++) {
			int cost;
			if (a.charAt(i) == b.charAt(j)) {
				cost = 0;
			} else {
				cost = 1;
			}
			distance_mat[i][j] = minimum(distance_mat[i - 1][j], distance_mat[i - 1][j - 1], distance_mat[i][j - 1]) + cost;


		}

	}
	return distance_mat[len_a - 1][len_b - 1];

}

public static String lcsString(String str1, String str2, int[][] lcsMatrix) {
	StringBuilder lcs = new StringBuilder();
	int i = str1.length(),
			j = str2.length();
	while (i > 0 && j > 0) {
		if (str1.charAt(i - 1) == str2.charAt(j - 1)) {
			lcs.append(str1.charAt(i - 1));
			i--;
			j--;
		} else if (lcsMatrix[i - 1][j] > lcsMatrix[i][j - 1]) {
			i--;
		} else {
			j--;
		}
	}
	return lcs.reverse().toString();
}


public static String getLCS(String str1, String str2) {

	//At least one string is null
	if (str1 == null || str2 == null)
		return null;

	//At least one string is empty
	if (str1.length() == 0 || str2.length() == 0)
		return "";

	String[] arr1 = str1.split("");
	String[] arr2 = str2.split("");

	//lcsMatrix[i][j]  = LCS of first i elements of arr1 and first j characters of arr2
	int[][] lcsMatrix = new int[arr1.length + 1][arr2.length + 1];

	for (int i = 0; i < arr1.length + 1; i++)
		lcsMatrix[i][0] = 0;
	for (int j = 1; j < arr2.length + 1; j++)
		lcsMatrix[0][j] = 0;
	for (int i = 1; i < arr1.length + 1; i++) {
		for (int j = 1; j < arr2.length + 1; j++) {
			if (arr1[i - 1].equals(arr2[j - 1])) {
				lcsMatrix[i][j] = lcsMatrix[i - 1][j - 1] + 1;
			} else {
				lcsMatrix[i][j] = lcsMatrix[i - 1][j] > lcsMatrix[i][j - 1] ? lcsMatrix[i - 1][j] : lcsMatrix[i][j - 1];
			}
		}
	}
	return lcsString(str1, str2, lcsMatrix);
}


public static int getLongestValidParentheses(String s) {
	if (s == null || s.length() < 2) {
		return 0;
	}
	char[] chars = s.toCharArray();
	int n = chars.length;
	int[] res = new int[n];
	res[0] = 0;
	res[1] = chars[1] == ')' && chars[0] == '(' ? 2 : 0;

	int max = res[1];

	for (int i = 2; i < n; ++i) {
		if (chars[i] == ')') {
			if (chars[i - 1] == '(') {
				res[i] = res[i - 2] + 2;
			} else {
				int index = i - res[i - 1] - 1;
				if (index >= 0 && chars[index] == '(') {
					// ()(())
					res[i] = res[i - 1] + 2 + (index - 1 >= 0 ? res[index - 1] : 0);
				}
			}
		}
		max = Math.max(max, res[i]);
	}

	return max;

}


public static void matrixChainOrder() {
	for (int i = 1; i < size + 1; i++) {
		m[i][i] = 0;
	}

	for (int l = 2; l < size + 1; l++) {
		for (int i = 1; i < size - l + 2; i++) {
			int j = i + l - 1;
			m[i][j] = Integer.MAX_VALUE;

			for (int k = i; k < j; k++) {
				int q = m[i][k] + m[k + 1][j] + p[i - 1] * p[k] * p[j];
				if (q < m[i][j]) {
					m[i][j] = q;
					s[i][j] = k;
				}
			}
		}
	}
}


public static int[] subset(int arr[],int sum)
 {
	int n = arr.length;
		boolean dp[][] = new boolean[n+1][sum+1];
		for(int i = 0; i <= n; i++)
			dp[i][0] = true;
		for(int i = 1; i <= sum; i++)
			dp[0][i] = false;
		// subset sum concept
		for(int i = 1; i <= n; i++)
		{
			for(int j = 1; j <= sum; j++)
			{
					if(arr[i-1] <= j)
						dp[i][j] = dp[i-1][j-arr[i-1]] || dp[i-1][j];
					else
						dp[i][j] = dp[i-1][j];
			}
		}
		//storing last dp column whose value is true till sum/2
		int index[] = new int[sum];
		int p = 0;
		for(int i = 0 ; i <= sum / 2; i++)
		{
			if(dp[n][i] == true)
					index[p++] = i;
		}
		return index;
 }
 
 
 public static int cutRod(int[] price, int n) {
	int val[] = new int[n + 1];
	val[0] = 0;

	for (int i = 1; i <= n; i++) {
		int max_val = Integer.MIN_VALUE;
		for (int j = 0; j < i; j++)
			max_val = Math.max(max_val, price[j] + val[i - j - 1]);

		val[i] = max_val;
	}

	return val[n];
}


public static int absMax(int[] numbers) {
	int absMaxValue = numbers[0];
	for (int i = 1, length = numbers.length; i < length; ++i) {
		if (Math.abs(numbers[i]) > Math.abs(absMaxValue)) {
			absMaxValue = numbers[i];
		}
	}
	return absMaxValue;
}

public static int absMin(int[] numbers) {
	int absMinValue = numbers[0];
	for (int i = 1, length = numbers.length; i < length; ++i) {
		if (Math.abs(numbers[i]) < Math.abs(absMinValue)) {
			absMinValue = numbers[i];
		}
	}
	return absMinValue;
}


public static int absVal(int value) {
	return value < 0 ? -value : value;
}

public static int aliquotSum(int number) {
	int sum = 0;
	for (int i = 1, limit = number / 2; i <= limit; ++i) {
		if (number % i == 0) {
			sum += i;
		}
	}
	return sum;
}

public static int recursiveCalcOfDividerSum(int number, int div) {

	if (div == 1) {
		return 0;
	} else if (number % --div == 0) {
		return recursiveCalcOfDividerSum(number, div) + div;
	} else {
		return recursiveCalcOfDividerSum(number, div);
	}
}


public static void findAllInRange(int startValue, int stopValue) {

	/* the 2 for loops are to avoid to double check tuple. For example (200,100) and (100,200) is the same calculation
	 * also to avoid is to check the number with it self. a number with itself is always a AmicableNumber
	 * */
	StringBuilder res = new StringBuilder();
	int countofRes = 0;

	for (int i = startValue; i < stopValue; i++) {
		for (int j = i + 1; j <= stopValue; j++) {
			if (isAmicableNumber(i, j)) {
				countofRes++;
				res.append("" + countofRes + ": = ( " + i + "," + j + ")" + "\t");
			}
		}
	}
	res.insert(0, "Int Range of " + startValue + " till " + stopValue + " there are " + countofRes + " Amicable_numbers.These are \n ");
	System.out.println(res.toString());
}

public static int average(int[] array) {
	long sum = 0;
	for (int i = 0 ; i < array.length; ++i) {
		sum += array[i];
	}
	return (int)(sum / array.length);
}

public static double ceil(double number) {
	if (number - (int) number == 0) {
		return number;
	} else if (number - (int) number > 0) {
		return (int) (number + 1);
	} else {
		return (int) number;
	}
}

public static int maxRecursion(int[] array, int low, int high) {
	if (low == high) {
		return array[low]; //or array[high]
	}

	int mid = (low + high) >>> 1;

	int leftMax = max(array, low, mid); //get max in [low, mid]
	int rightMax = max(array, mid + 1, high); //get max in [mid+1, high]

	return Math.max(leftMax, rightMax);
}

public static int minRecursion(int[] array, int low, int high) {
	if (low == high) {
		return array[low]; //or array[high]
	}

	int mid = (low + high) >>> 1;

	int leftMin = min(array, low, mid); //get min in [low, mid]
	int rightMin = min(array, mid + 1, high); //get min in [mid+1, high]

	return Math.min(leftMin, rightMin);
}

public static int[] mode(int[] numbers) {
	
	if(numbers.length == 0) return null;
	
	HashMap<Integer, Integer> count = new HashMap<>();
	
	for(int num : numbers) {
		if(count.containsKey(num)) {
			
			count.put(num, count.get(num) + 1);
			
		} else {
			
			count.put(num, 1);
			
		}
		
	}
	
	int max = Collections.max(count.values());
	ArrayList<Integer> modes = new ArrayList<>();
	
	for(int num : count.keySet()) {
		if(count.get(num) == max) {
			modes.add(num);
		}
	}
	return modes.stream().mapToInt(n -> n).toArray();
}


public static String splitIntoDigits(int num, int num2) {

	StringBuilder res = new StringBuilder();

	ArrayList<Integer> digits = new ArrayList<>();
	while (num > 0) {
		digits.add(num%10);
		num /= 10;
	}
	while (num2 > 0) {
		digits.add(num2%10);
		num2/= 10;
	}
	Collections.sort(digits);
	for (int i : digits) {
		res.append(i);
	}


	return res.toString();
}

public static boolean isPythagTriple(int a, int b, int c) {
	if(a <= 0 || b <= 0 || c <= 0) {
		return false;
	} else {
		return (a * a) + (b * b) == (c * c);
	}
}

public static boolean isAbecedarian(String s) {
	int index = s.length() - 1;

	for (int i = 0; i < index; i++) {

		if (s.charAt(i) <= s.charAt(i + 1)) {
		} //Need to check if each letter for the whole word is less than the one before it

		else {
			return false;
		}
	}
	return true;
}

public static int findWorstFit(int[] blockSizes, int processSize) {
	int max = -1;
	int index = -1;
	for(int i=0 ; i < blockSizes.length ; i++) { // Find the index of the biggest memory block available.
		if(blockSizes[i] > max) {
			max = blockSizes[i];
			index = i;
		}
	}
	// If the biggest memory block cannot fit the process, return -255 as the result
	if(processSize > blockSizes[index]) {
		return NO_ALLOCATION;
	}
	return index;
}


public static Test[] joinTestArray(Test[] one, Test[] two) {
	Test[] output = new Test[one.length + two.length];
	int i = 0;
	for (i = 0; i < one.length; i++) {
		output[i] = one[i];
	}         
	int j = 0;
	for (j = 0; j < two.length; j++) {
		output[i+j] = two[j];
	}
	return output;    
}

public int getIndexFromPile (Card toFind, List<Card> pile) {
   int index = -1;
   boolean found = false;
   for (int i = 0; i < pile.size() && !found; i++) {
	  //System.out.println("Finding: " + toFind + " Found: " + pile.get(i));
	  if (pile.get(i).equals(toFind)) {
		 found = true;
		 index = i;
	  }
   }
   return index;
}


public static boolean isSafe(int v, boolean[][] arr, int[] col, int c,int V) {

	// TODO Auto-generated method stub

	for(int i=0;i<V;i++)

	{

		if(arr[v][i]&&c==col[i])

		{

			return false;

		}

	}

	return true;

}


public static boolean colorpossible(boolean[][] arr, int m, int[] col, int v,int V) {
		if(v==V)
		{
			return true;
		}
		for(int i=1;i<=m;i++)
		{
			if(isSafe(v,arr,col,i,V))

			{

				col[v]=i;

				if(colorpossible(arr, m, col, v+1, V))

				{
					return true;
				}
				col[v]=0;
			}
		}
		return false;
}



public static void createGraph(Scanner in) {

	int v = in.nextInt();
	adj = new LinkedList[v];

	for (int i = 0; i < v; i++) {
		adj[i] = new LinkedList<>();
	}

	int e = in.nextInt();

	for (int i = 0; i < e; i++) {
		int v1 = in.nextInt();
		int v2 = in.nextInt();
		adj[v1].add(v2);
	}
}


public static void bfs(LinkedList<Integer>[] adj, int s) {
        boolean marked[] = new boolean[adj.length];
        int[] dist = new int[adj.length];
        marked[s] = true;
        dist[s] = 0;
        Queue<Integer> queue = new LinkedList<>();
        queue.add(s);
        while (!queue.isEmpty()) {
            int v = queue.poll();
            for (Integer i : adj[v]) {
                if (!marked[i]) {
                    marked[i] = true;
                    dist[i] = dist[v] + 1;
                    queue.add(i);
                }
            }
        }
        System.out.println(Arrays.toString(dist));
    }
	
public static void dfs(int i) {
	for (int v : adj[i]) {
		if (!marked[v]) {
			marked[v] = true;
			dist[v] = dist[i] + 1;
			dfs(v);
		}
	}
}


public static int findMinDiff(int[] arr, int m) {
	if (arr == null || arr.length == 0) {

		return -1;
	}

	if (m > arr.length) {

		return -1;
	}

	Arrays.sort(arr);
	int minDiff = Integer.MAX_VALUE;

	for (int i = 0; i + m - 1 < arr.length; i++) {

		int diff = arr[i + m - 1] - arr[i];

		if (diff < minDiff) {

			minDiff = diff;
		}
	}

	return minDiff;
}


public static int tappingWater(int[] heights) {

	int n = heights.length;
	int left[] = new int[n];
	int right[] = new int[n];
	left[0] = heights[0];

	for (int i = 1; i < n; i++) {
		left[i] = Math.max(left[i - 1], heights[i]);
	}

	right[n - 1] = heights[n - 1];
	for (int i = n - 2; i >= 0; i--) {
		right[i] = Math.max(right[i + 1], heights[i]);
	}

	int water = 0;
	for (int i = 0; i < n; i++) {
		water += Math.min(left[i], right[i]) - heights[i];
	}
	return water;

}


public void getAllWordsStartingWith(TrieNode node, StringBuilder prefix, Set<String> words) {

	if (node.isEnd) {
		words.add(prefix.toString());
	}

	for (int i = 0; i < 26; i++) {
		if (node.links[i] != null) {
			prefix.append((char) ('a' + i));
			getAllWordsStartingWith(node.links[i], prefix, words);
			prefix.deleteCharAt(prefix.length() - 1);
		}
	}
}


public static int countWays(int steps, int[] memo, int count) {
        if (steps == 1 || steps == 2) {
            return steps;
        }
        if (memo[steps] != 0) {
            return memo[steps];
        }
        count = countWays(steps - 1, memo, count) + countWays(steps - 2, memo, count);
        memo[steps] = count;
        return count;
    }
	
	
public boolean isLexicographicallyValid(String leftWord, String rightWord, int[] idxOrder) {
	int i = 0, j = 0;
	for (; i < leftWord.length() && j < rightWord.length(); i++, j++) {
		char leftChar = leftWord.charAt(i);
		char rightChar = rightWord.charAt(j);

		if (idxOrder[leftChar - 'a'] > idxOrder[rightChar - 'a']) {
			return false;
		} else if (idxOrder[leftChar - 'a'] < idxOrder[rightChar - 'a']) {
			return true;
		}
	}

	if (i < leftWord.length()) {
		return false;
	}
	return true;
}



public int[] convertToIdxArray(String order) {
	int[] idxOrder = new int[26];
	for (int i = 0; i < 26; i++) {
		idxOrder[order.charAt(i) - 'a'] = i;
	}
	return idxOrder;
}


public static int numberOfArithmeticSlices(int[] A) {
	int count = 0;
	for (int i = 0; i + 2 < A.length; i++) {
		int diff = A[i + 1] - A[i];
		int j = i + 2;
		while (j < A.length && A[j] - A[j - 1] == diff) {
			count++;
			j++;
		}
	}
	return count;
}


public static List<String> possibleMinutes(int num) {
	List<String> allMinutes = new ArrayList<>();
	for (int i = 0; i < 59; i++) {
		if (countOne(Integer.toBinaryString(i)) == num) {
			allMinutes.add(String.format("%02d", i));
		}
	}
	return allMinutes;
}


public static int maxKilledEnemies(char[][] grid) {
	int m = grid.length, n = grid[0].length;
	int maxHits = 0, rowHits = 0;
	int[] colHits = new int[n];

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			if (j == 0 || grid[i][j - 1] == 'W') {
				rowHits = 0;
				for (int k = j; k < n && grid[i][k] != 'W'; k++) {
					rowHits += grid[i][k] == 'E' ? 1 : 0;
				}
			}

			if (i == 0 || grid[i - 1][j] == 'W') {
				colHits[j] = 0;
				for (int k = i; k < m && grid[k][j] != 'W'; k++) {

					colHits[j] += grid[k][j] == 'E' ? 1 : 0;
				}
			}

			if (grid[i][j] == 0) {
				maxHits = Math.max(maxHits, rowHits + colHits[j])
			}
		}
	}
	return maxHits;

}


public static int candy(int[] ratings) {
	if (ratings == null || ratings.length == 0) {
		return 0;
	}
	int n = ratings.length;
	int[] candies = new int[n];
	Arrays.fill(candies, 1);

	for (int i = 1; i < n; i++) {
		if (ratings[i] > ratings[i - 1]) {
			candies[i] = candies[i - 1] + 1;
		}
	}

	for (int i = n - 2; i >= 0; i--) {
		if (ratings[i] > ratings[i + 1]) {
			candies[i] = Math.max(candies[i], candies[i + 1] + 1);
		}
	}
	return Arrays.stream(candies).sum();

}


public static String convertToHex(int n) {

	String[] hex = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "a", "b", "c", "d", "e", "f"};
	if (n < 0) {
		return "";
	}
	String output = "";
	while (n > 0) {
		output = hex[n % 16] + output;
		n = n / 16;
	}
	return output;
}


public static int firstUniqChar(String s) {
	int[] counts = new int[26]; // char counts
	LinkedHashMap<Character, Integer> firstIndex = new LinkedHashMap<>(26); // first indexes
	for (int i = 0; i < s.length(); i++) {
		if (counts[s.charAt(i) - 'a']++ == 0) {
			firstIndex.put(s.charAt(i), i);
		}
	}

	for (Map.Entry<Character, Integer> entry : firstIndex.entrySet()) {
		int i = entry.getKey() - 'a';
		if (counts[i] == 1) {
			return entry.getValue();
		}
	}
	return -1;
}


public static int integerReplacement(int n) {
	int count = 0;
	while (n > 1) {
		if (n % 2 == 0) {
			n = n / 2;
		} else if (n == 3) {
			n--;
		} else if (Integer.bitCount(n - 1) < Integer.bitCount(n + 1)) {
			n--;
		} else {
			n++;
		}
		count++;
	}
	return count;
}


public static List<Integer> largestDivisibleSubset(int[] nums) {

	if (nums == null || nums.length == 0) {
		return new ArrayList<>();
	}
	Arrays.sort(nums);
	int max = 0;
	int maxInd = -1;
	int[] count = new int[nums.length];
	Arrays.fill(count, 1);
	int[] parent = new int[nums.length];

	for (int i = nums.length - 1; i >= 0; i--) {
		for (int j = i + 1; j < nums.length; j++) {
			if (nums[j] % nums[i] == 0
					&& count[j] + 1 > count[i]) {
				count[i] = count[j] + 1;
				parent[i] = j;
			}
		}
		if (count[i] > max) {
			max = count[i];
			maxInd = i;
		}
	}

	List<Integer> subset = new ArrayList<>();
	for (int i = 0; i < max; i++) {
		subset.add(nums[maxInd]);
		maxInd = parent[maxInd];
	}

	return subset;
}


public static List<Integer> lexicalOrder(int n) {
	List<Integer> result = new ArrayList<>(n);
	int val = 1;
	for (int i = 1; i <= n; i++) {
		result.add(val);
		if (val * 10 <= n) {
			val *= 10;
		} else if (val % 10 != 9 && val + 1 <= n) {
			val++;
		} else {
			while ((val / 10) % 10 == 9) {
				val /= 10;
			}
			val = val / 10 + 1;
		}

	}
 
	return result;
}

public static int getDepth(String subpath) {

	int count = 0, i = 0;
	while (i < subpath.length() - 1) {
		if (subpath.charAt(i) == '\\' && subpath.charAt(i + 1) == 't') {
			count++;
		}
		i++;
	}
	return count;
}


public static Map<Integer, Integer> createIndexes(char[] chars, char c) {

	Map<Integer, Integer> indexes = new HashMap<>();
	int start = -1;

	for (int j = 0; j < chars.length; j++) {
		if (chars[j] == c) {
			if (start == -1) {
				start = j;
			}
			indexes.put(start, j);
		} else 
			start = -1;
		}
	}
	return indexes;

}

public static int[] charCount(String s, int start, int end) {

	int[] charCount = new int[26];
	for (int i = start; i <= end; i++) {
		charCount[s.charAt(i) - 'a']++;
	}
	return charCount;
}

public List<Integer> generatePowers(int number, int bound) {

	if (number == 0 || number == 1) {
		ArrayList<Integer> list = new ArrayList<>();
		list.add(number);
		return list;
	}
	List<Integer> powers = new ArrayList<>();
	int x = 0;

	while (true) {
		int value = (int) Math.pow(number, x++);
		powers.add(value);
		if (value > bound) {
			break;
		}
	}
	return powers;
}

public static int[] getModifiedArray(int length, int[][] updates) {

	if (length <= 0) {
		return null;
	}

	int res[] = new int[length];
	for (int i = 0; i < updates.length; i++) {
		res[updates[i][0]] += updates[i][2];
		if (updates[i][1] + 1 < length) {
			res[updates[i][1] + 1] -= updates[i][2];
		}
	}

	for (int i = 1; i < length; i++) {
		res[i] += res[i - 1];
	}

	return res;

}

public boolean isBoardValid(int n, int[] board) {

	for (int row = 0; row < n; row++) {
		if (board[n] == board[row]) {
			return false; // same column, different rows
		}
		if (Math.abs(board[n] - board[row]) == (n - row)) {
			return false;
		}
	}
	return true;
}

public static boolean canConstruct(String ransomNote, String magazine) {

	int[] counts = new int[26];

	for (int i = 0; i < magazine.length(); i++) {
		counts[magazine.charAt(i) - 'a']++;
	}

	for (int i = 0; i < ransomNote.length(); i++) {
		counts[ransomNote.charAt(i) - 'a']--;
	}

	for (int i = 0; i < 26; i++) {
		if (counts[i] < 0) {
			return false;
		}
	}
	return true;
}

private void appendAlphabet(StringBuilder result, char alphabet, int times) {

	for (int i = 0; i < times; i++) {
		result.append(alphabet);
	}
}


public int repeatedNTimes(int[] A) {
	if (A == null || A.length == 0) {
		return -1;
	}

	Set<Integer> uniqueNumbers = new HashSet<>();
	for (int num : A) {
		if (uniqueNumbers.contains(num)) {
			return num;
		}
		uniqueNumbers.add(num);
	}
	return -1;
}

public static int maxRotateFunction(int[] A) {
	int n = A.length;
	int max = 0;

	for (int i = 0; i < n; i++) {
		int currMax = 0;
		for (int j = 0; j < n; j++) {
			currMax += A[(j + i) % n] * j;
		}
		max = Math.max(max, currMax);
	}
	return max;

}

public static int fitOnScreen(String sentence, int rows, int cols) {

String[] words = sentence.split(" ");
int i = 0, j = 0, wi = 0, count = 0;

	while (i < rows) {
		if (j == cols - 1 || j + words[wi].length() - 1 >= cols) {
			i++;
			j = 0;
		} else {
			j += words[wi].length();
			j++; // for spac

			if (wi == words.length - 1) {
				count++;
				wi = 0;
			} else {
				wi++;
			}
		}
	}
	return count;
}

public int indexOfSmallestNegative(int[] sortedArray) {
	int idx = 0;
	while (idx < sortedArray.length && sortedArray[idx] < 0) {
		idx++;
	}
	return idx - 1;
}


public static int countFalses(boolean... values) {
	int count = 0;
	for (boolean b : values) {
		if (!b) {
			count++;
		}
	}
	return count;

}


public static int tripletCount(String str) {

	int triplets = 0;
	int i = 0;

	while (i + 2 < str.length()) {
		char ch = str.charAt(i);
		if (ch == str.charAt(i + 1)
			&& ch == str.charAt(i + 2)) {
			triplets++;
			i += 3;
		} else {
			i++;
		}
	}
	return triplets;
}


public static boolean hasLowerCase(String str) {

	for (int i = 0; i < str.length(); i++) {
		if (str.charAt(i) >= 'a' && str.charAt(i) <= 'z') {
			return true;
		}
	}
	return false;
}



public static boolean hasUpperCase(String str) {

	for (int i = 0; i < str.length(); i++) {
		if (str.charAt(i) >= 'A' && str.charAt(i) <= 'Z') {
			return true;
		}
	}
	return false;
}



public static boolean hasDigit(String str) {
	for (int i = 0; i <= 9; i++) {
		if (str.contains(Integer.toString(i))) {
			return true;
		}
	}
	return false;
}