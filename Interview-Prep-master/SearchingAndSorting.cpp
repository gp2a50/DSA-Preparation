#include <iostream>
#include <vector>
#include <unordered_map>
#include <list>
#include <set>
#include <map>
#include <algorithm>
#include <math.h>

using namespace std;

// 169. Majority Element (Boyer-Moore Voting Algorithm)
/*
In case the majority element always exists, then sort the array and return arr[n/2]
Otherwise
Method 1- Brute Force
Method 2- Hashmap to keep freq of elements, and for each element check if their freq is >n/2: O(n),O(n)
Method 3- Use self balancing BST, similar to ordered map: O(nlogn), O(n), for normal BST O(n^2) time
Method 4- Sort the array, then keep track of count of elements: O(nlogn), O(1)
Method 5- Moore's Voting Algorithm: O(n), O(1)
*/
// Boyer-Moore Voting Algorithm
int majorityElement(vector<int> &nums)
{
    //find candidate
    int me = nums[0], count = 1;
    for (int i = 1; i < nums.size(); i++)
    {
        if (nums[i] == me)
            count++;
        else
            count--;

        if (count == 0)
        {
            me = nums[i];
            count = 1;
        }
    }

    //verify
    count = 0;
    for (int i = 0; i < nums.size(); i++)
        if (nums[i] == me)
            count++;

    if (count > nums.size() / 2)
        return me;

    return -1;
}
//Using Hashmap or BST
int majorityElement(vector<int> &nums)
{
    int n = nums.size();
    unordered_map<int, int> freq;
    for (int i = 0; i < n; i++)
    {
        freq[nums[i]]++;
        if (freq[nums[i]] > n / 2)
            return nums[i];
    }

    return -1;
}
//Incase it always exists
int majorityElement(vector<int> &nums)
{
    sort(nums.begin(), nums.end());
    return nums[nums.size() / 2];
}

// 229. Majority Element II
//Approach 1- Hashmap
vector<int> majorityElement2(vector<int> &nums)
{
    vector<int> res;
    int n = nums.size();
    unordered_map<int, int> freq;
    for (int i = 0; i < n; i++)
    {
        if (freq[nums[i]] != -1)
        {
            freq[nums[i]]++;
            if (freq[nums[i]] > n / 3)
            {
                res.push_back(nums[i]);
                freq[nums[i]] = -1;
            }
        }
    }

    return res;
}
//Approach 2- Modified Boyer Moore Voting Algorithm
vector<int> majorityElement2(vector<int> &nums)
{
    int n = nums.size();
    vector<int> res;
    int c1 = 0, c2 = 0, count1 = 0, count2 = 0;

    for (int i = 0; i < n; i++)
    {
        if (nums[i] == c1)
            count1++;
        else if (nums[i] == c2)
            count2++;
        else if (count1 == 0)
        {
            c1 = nums[i];
            count1 = 1;
        }
        else if (count2 == 0)
        {
            c2 = nums[i];
            count2 = 1;
        }
        else
        {
            count1--;
            count2--;
        }
    }

    count1 = 0;
    count2 = 0;
    for (int i = 0; i < n; i++)
    {
        if (nums[i] == c1)
            count1++;
        else if (nums[i] == c2)
            count2++;
    }

    if (count1 > n / 3)
        res.push_back(c1);
    if (count2 > n / 3)
        res.push_back(c2);

    return res;
}

// Searching in an array where adjacent differ by at most k
/*
A Simple Approach is to traverse the given array one by one and compare every element with given element ‘x’.
If matches, then return index.
The above solution can be Optimized using the fact that difference between all adjacent elements is at most k.
The idea is to start comparing from the leftmost element and find the difference between current array element and x.
Let this difference be ‘diff’. From the given property of array,
we always know that x must be at-least ‘diff/k’ away, so instead of searching one by one, we jump ‘diff/k’.
*/
int searchElement(vector<int> &nums, int k, int tar)
{
    int i = 0;
    while (i < nums.size())
    {
        // If x is found at index i
        if (nums[i] == tar)
            return i;

        // Jump the difference between current
        // array element and tar divided by k
        // We use max here to make sure that i
        // moves atleast one step ahead.
        i = i + max(1, abs((nums[i] - tar) / k));
    }

    //not found
    return -1;
}

// Find Missing And Repeating
/*
7 Methods-
https://www.geeksforgeeks.org/find-a-repeating-and-a-missing-number/

My method-
for all numbers mark arr[abs(arr[i])-1] -ve
if while marking -ve the target is already -ve, then that is repeating number
then for missing, iterate again, the index that has +ve element gives missing number
*/
int *findTwoElement(int *arr, int n)
{
    int missing, repeating;
    int *res = new int(2);
    for (int i = 0; i < n; i++)
    {
        if (arr[abs(arr[i]) - 1] < 0)
            repeating = abs(arr[i]);
        else
            arr[abs(arr[i]) - 1] *= -1;
    }

    for (int i = 0; i < n; i++)
    {
        if (arr[i] > 0)
        {
            missing = i + 1;
            break;
        }
    }

    res[0] = repeating;
    res[1] = missing;
    return res;
}

// Find four elements that sum to a given value
/*
Approach 1-Recursion(Subsequence method) O(2^n)

Approach 2- (Best) O(NlogN), (N=n^2)
find all pair sums (n^2)
sort all pair sums, and use two pointer approach to find pairSum: O(NlogN)

Approach 3- O(n^2)
Make a hashmap of all pair sums
for elements in map, check if Sum-X exists in map, if it exists, then check if all 4 indexes are unique
to avoid duplicate quadruples, we use a set.
Now for an ordered_Set, complexity will be logn for insert, so overall complexity is not O(n^2)
*/
void fourElementSum(vector<int> &nums, int tar)
{
    int n = nums.size();
    unordered_map<int, vector<pair<int, int>>> mp; //{pair sum: list of indexes}
    set<vector<int>> res;

    //group all index pairs having same pair sum
    for (int i = 0; i < n - 1; i++)
    {
        for (int j = i + 1; j < n; j++)
            mp[nums[i] + nums[j]].push_back({i, j});
    }

    for (int i = 0; i < n - 1; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            //for each pair check give tar-pairSum exists in map
            int pairSum = nums[i] + nums[j];
            if (mp.find(tar - pairSum) != mp.end())
            {
                //if it exists, then pair it with all the indexes with that sum
                for (pair<int, int> ele : mp[tar - pairSum])
                {
                    int a = i;
                    int b = j;
                    int c = ele.first;
                    int d = ele.second;
                    //check that all indexes are unique
                    if (a != c && a != d && b != c && b != d)
                    {
                        //add them all in a vector and sort the vector
                        vector<int> smallAns;
                        smallAns.push_back(nums[a]);
                        smallAns.push_back(nums[b]);
                        smallAns.push_back(nums[c]);
                        smallAns.push_back(nums[d]);
                        sort(smallAns.begin(), smallAns.end());
                        //add it to a set so that duplicate ans are removed
                        res.insert(smallAns);
                    }
                }
            }
        }
    }

    if (res.size() == 0)
        cout << -1;

    for (auto ans : res)
    {
        for (auto ele : ans)
        {
            cout << ele << " ";
        }
        cout << "$";
    }
}

// 18. 4Sum
vector<vector<int>> fourSum(vector<int> &nums, int tar)
{
    int n = nums.size();
    unordered_map<int, vector<pair<int, int>>> mp; //{pair sum: list of indexes}
    set<vector<int>> res;

    //group all index pairs having same pair sum
    for (int i = 0; i < n - 1; i++)
    {
        for (int j = i + 1; j < n; j++)
            mp[nums[i] + nums[j]].push_back({i, j});
    }

    for (int i = 0; i < n - 1; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            //for each pair check give tar-pairSum exists in map
            int pairSum = nums[i] + nums[j];
            if (mp.find(tar - pairSum) != mp.end())
            {
                //if it exists, then pair it with all the indexes with that sum
                for (pair<int, int> ele : mp[tar - pairSum])
                {
                    int a = i;
                    int b = j;
                    int c = ele.first;
                    int d = ele.second;
                    //check that all indexes are unique
                    if (a != c && a != d && b != c && b != d)
                    {
                        //add them all in a vector and sort the vector
                        vector<int> smallAns;
                        smallAns.push_back(nums[a]);
                        smallAns.push_back(nums[b]);
                        smallAns.push_back(nums[c]);
                        smallAns.push_back(nums[d]);
                        sort(smallAns.begin(), smallAns.end());
                        //add it to a set so that duplicate ans are removed
                        res.insert(smallAns);
                    }
                }
            }
        }
    }

    vector<vector<int>> ans;
    for (auto ele : res)
    {
        ans.push_back(ele);
    }

    return ans;
}

// Maximum sum such that no two elements are adjacent
/*
Use the recursion and apply dp
*/
long long maxSubsequenceSum(vector<int> &arr, int idx, vector<int> &dp)
{
    if (idx >= arr.size())
        return dp[idx] = 0;

    if (dp[idx] != -1)
        return dp[idx];

    //including the current element
    long long including = maxSubsequenceSum(arr, idx + 2, dp) + arr[idx];
    //excluding the current element
    long long excluding = maxSubsequenceSum(arr, idx + 1, dp);

    //return the max of inlcuding and excluding
    return dp[idx] = max(including, excluding);
}
long long maxSubsequenceSum_DP(vector<int> &arr)
{
    int n = arr.size();
    vector<int> dp(n + 2, -1);

    for (int idx = n + 1; idx >= 0; idx--)
    {
        if (idx >= arr.size())
        {
            dp[idx] = 0;
            continue;
        }

        //including the current element
        long long including = dp[idx + 2] + arr[idx];
        //excluding the current element
        long long excluding = dp[idx + 1];

        //return the max of inlcuding and excluding
        dp[idx] = max(including, excluding);
    }

    return dp[0];
}
long long maxSubsequenceSum(vector<int> &arr)
{
    int n = arr.size();
    vector<int> dp(n + 2, -1);
    return maxSubsequenceSum(arr, 0, dp);
}

// 349. Intersection of Two Arrays
/*
Approach 1- O(n+m) time, O(n) space
Use Hashmap to get count of elements in first array
then iterate over second array and check if it is present in map. 
If it is present, add it to result and remove it from map to avoid adding it again

Approach 2-O(1) space O(nlogn + mlogm) time
sort both arrays, and use two pointers to find the intersection
*/
vector<int> intersection(vector<int> &nums1, vector<int> &nums2)
{
    int n = nums1.size();
    int m = nums2.size();
    vector<int> res;
    sort(nums1.begin(), nums1.end());
    sort(nums2.begin(), nums2.end());

    int i = 0, j = 0;
    while (i < n && j < m)
    {
        //if it is equal , then add to result
        if (nums1[i] == nums2[j])
        {
            res.push_back(nums1[i]);
            //skip all other equal elements as it has been added already
            while (i < n - 1 && nums1[i] == nums1[i + 1])
                i++;
            while (j < m - 1 && nums2[j] == nums2[j + 1])
                j++;
            i++;
            j++;
        }
        else if (nums1[i] < nums2[j])
            i++;
        else
            j++;
    }

    return res;
}

// 350. Intersection of Two Arrays II
/*
Similar approaches as previous, just add the duplicates as well
*/
//Approach 1-Sorting (O(1) space)
vector<int> intersect(vector<int> &nums1, vector<int> &nums2)
{
    int n = nums1.size();
    int m = nums2.size();
    vector<int> res;

    sort(nums1.begin(), nums1.end());
    sort(nums2.begin(), nums2.end());

    int i = 0, j = 0;
    while (i < n && j < m)
    {
        if (nums1[i] == nums2[j])
        {
            res.push_back(nums1[i]);
            i++;
            j++;
        }
        else if (nums1[i] < nums2[j])
            i++;
        else
            j++;
    }

    return res;
}
//Approach 2- Hashmap (O(n) space)
vector<int> intersect(vector<int> &nums1, vector<int> &nums2)
{
    vector<int> res;
    unordered_map<int, int> freq;
    //count frequency of each element in nums1
    for (int ele : nums1)
        freq[ele]++;

    //for each element in nums2, if its count in map >0, add it to result, and decrease count by 1
    for (int ele : nums2)
    {
        if (freq[ele] > 0)
        {
            res.push_back(ele);
            freq[ele]--;
        }
    }

    return res;
}

// Find common elements in three sorted arrays (G4g test cases are wrong)
/*
Approach 1- Hashmap
Find intersection of first two, then find intersection of res and third array
Approach 2- Sorting
Using 3 pointers, similar to 2 pointers in 2 arrays
*/
vector<int> commonElements(int A[], int B[], int C[], int n1, int n2, int n3)
{
    vector<int> res;

    sort(A, A + n1);
    sort(B, B + n2);
    sort(C, C + n3);

    int i = 0, j = 0, k = 0;
    while (i < n1 && j < n2 && k < n3)
    {
        if (A[i] == B[j] && B[j] == C[k])
        {
            res.push_back(A[i]);
            while (i < n1 - 1 && A[i] == A[i + 1])
                i++;
            while (j < n2 - 1 && B[j] == B[j + 1])
                j++;
            while (k < n3 - 1 && C[k] == C[k + 1])
                k++;
            i++;
            j++;
            k++;
        }
        else if (A[i] < B[j])
            i++;
        else if (B[j] < C[k])
            j++;
        else
            k++;
    }

    return res;
}

// Count triplets with sum smaller than a given value
/*
Approach 1- Brute Force O(n^3)

Approach 2- O(nlogn + n^2)
Sort the array
Then use on for loop , for the first element
inside the for loop find the other two using 2 pointers approach
*/
long long countTriplets(long long arr[], int n, long long sum)
{
    sort(arr, arr + n);

    long long count = 0;

    for (int k = 0; k < n; k++)
    {
        int ele1 = arr[k];

        //use meet in the middle approach
        int i = k + 1, j = n - 1;
        while (i < j)
        {
            //sum is less so it is part of the ans
            //if for a fixed k, (i,j) gives sum(i,j,k)<target, then every element before j, will also give <sum with i,k
            //as in sorted array everything before j is smaller, so sum will keep decreasing for fixed i,k.
            if (arr[i] + arr[j] + arr[k] < sum)
            {
                //include all indexes from i+1 to in the ans, as they can all pair up with i,k
                count += j - i;
                i++;
            }
            //mySum>=target, so decrease j
            else
                j--;
        }
    }

    return count;
}

// Print all subarrays with 0 sum (count the subarrays)
int zeroSumSubarray(vector<int> &arr)
{
    unordered_map<long long, int> mp; //prefixSum:count
    long long currSum = 0;
    int count = 0;

    for (int i = 0; i < arr.size(); i++)
    {
        currSum += arr[i];
        if (currSum == 0)
            count++;
        if (mp.find(currSum) != mp.end())
            count += mp[currSum];
        mp[currSum]++;
    }

    return count;
}

// 462. Minimum Moves to Equal Array Elements II (Make all array elements equal with minimum cost)
/*
The quickest way to make them all equal is to make them all equal to the median
Median is the element in the middle of the sorted array
*/
int minMoves2(vector<int> &nums)
{
    int n = nums.size();

    //sort and get median
    sort(nums.begin(), nums.end());
    int median = nums[n / 2], count = 0;

    //find the count
    for (int i = 0; i < n; i++)
        count += abs(median - nums[i]);

    return count;
}

// Check if reversing a sub array make the array sorted (no submit option)
/*
for the array to be sorted every arr[i]>arr[i-1]
so we keep going until array is increasing
when it starts decreasing, we keep going,
and then for the last element of the decreasing subarray, we check if it is > last element of increasing subarray
and the first element of decreasing subarray is < first of next increasing subarray
and everything after that is also sorted

For eg- {1,2,3,4,20,9,16,17} 
till 4 it is first part
then 20,9 is second part
but 20>16 so false

Similarly like, {1,2,3,4,10,3,16,17}
now after rotation 3<4 so false
*/
bool checkSubarray(vector<int> &arr)
{
    int n = arr.size();

    int i = 0;
    //check for the first increasing subarray
    while (i < n - 1)
    {
        if (arr[i] > arr[i + 1])
            break;
        i++;
    }

    int lastSortedEle = arr[i - 1]; //last element in the first part
    int firstUnsorted = arr[i];

    //check for decreasing subarray, which can be rotated
    while (i < n - 1)
    {
        if (arr[i] < arr[i + 1])
            break;
        i++;
    }

    //check if the last element of second part is greater than the last of first part
    //or first element of second is > first element of third part
    //i.e. after rotation, it will be sorted
    if (arr[i] < lastSortedEle || (i < n && firstUnsorted > arr[i + 1]))
        return false;

    //check for the increasing subarray after the rotated subarray
    while (i < n - 1)
    {
        if (arr[i] > arr[i + 1])
            break;
        i++;
    }

    return i == n;
}

// A Product Array Puzzle
/*
Given an array arr[] of n integers, construct a Product Array prod[] (of same size) 
such that prod[i] is equal to the product of all the elements of arr[] except arr[i].

Approach 1-
Get the total product of array. For each element res[i]=totalProduct/arr[i]

Without using Division-
Approach 1 - Prefix and suffix array
Make a prefix and suffix array
then for each res[i]=prefixProduct[i-1] * suffixProduct[i+1]

Approach 2 - Instead of taking O(2n) space for prefix and suffix, use a temp variable to store the prefix and suffix product
Do same as prefix and suffix, But build the prefix and suffix array in the res array itself,
just keep a temp variable that stores the prefix product for that index
then use the temp variable to keep the suffix product
*/
void productArray(vector<int> &arr)
{
    int n = arr.size();
    long long temp = 1;
    vector<long long> res(n);

    //get the prefix product of 0...i-1
    for (int i = 0; i < n; i++)
    {
        res[i] = temp;
        temp *= arr[i];
    }

    //multiply the prefix products with suffix product of i+1....n-1
    temp = 1;
    for (int i = n - 1; i >= 0; i--)
    {
        res[i] *= temp;
        temp *= arr[i];
    }

    //display the result
    for (int i = 0; i < n; i++)
        cout << res[i] << " ";
    cout << endl;
}

// 4. Median of Two Sorted Arrays
int getMax(vector<int> &nums, int idx) //find max of left region of partition
{
    if (idx == 0)
        return -1e8;
    return nums[idx - 1];
}
int getMin(vector<int> &nums, int idx) //find min of right region of partition
{
    if (idx == nums.size())
        return 1e8;
    return nums[idx];
}
double findMedianSortedArrays(vector<int> &nums1, vector<int> &nums2)
{
    if (nums1.size() > nums2.size())
        return findMedianSortedArrays(nums2, nums1);

    int lo = 0;
    int hi = nums1.size();
    int totalLength = nums1.size() + nums2.size();

    while (lo <= hi)
    {
        //get the partition for both arrays
        int partX = (lo + hi) / 2;
        int partY = (totalLength + 1) / 2 - partX; //+1 to account for odd and even, as in odd (7+1/2)=4 and in even (8+1)/2=4 too

        //get the leftMax and rightMin for both
        int leftX = getMax(nums1, partX);
        int leftY = getMax(nums2, partY);

        int rightX = getMin(nums1, partX);
        int rightY = getMin(nums2, partY);

        //if partition is correct
        if (leftX <= rightY && leftY <= rightX)
        {
            //avg in case of even length
            if (totalLength % 2 == 0)
                return (max(leftX, leftY) + min(rightX, rightY)) / 2.0;
            //leftMax in case of odd length
            return max(leftX, leftY);
        }
        //lies in left region
        else if (leftX > rightY)
            hi = partX - 1;
        //lies in right region
        else
            lo = partX + 1;
    }

    return -1;
}

// Merge Without Extra Space
/*
Approach 1 - Merge using extra space O(n+m) time, O(n+m) space

Approach 2 - Insertion sort - O(n*m) time O(1) space
iterate through first array, if arr1[i]>arr2[j] then swap them
and then put the new element in the correct position in arr2.
For eg-
1 7 9 11
3 4 6 8

so 1<3 move on
7>3 -> swap(7,3) -> arr2= 7 4 6 8, now put 7 in its correct position using linear traversal
for that you will have to move all other elements <7 to left, so it O(m) time for this step

Approach 3 - Shell Sort kind of algo to merge in O(nlogn), O(1) time,space
Find the gap
swap elements according to gap in 1st array,
swap elements according to gap in 1st and 2nd when gap includes element from both,
swap in second array
update gap
*/
void merge(int arr1[], int arr2[], int n, int m)
{
    // int gap = (n + m + 1) / 2; //gives ceil in case of odd, and normal number in even  case

    for (int gap = (n + m + 1) / 2; gap > 0; gap = (gap == 1) ? 0 : (gap + 1) / 2)
    {
        int i;
        for (i = 0; i + gap < n; i++)
        {
            if (arr1[i] > arr1[i + gap])
                swap(arr1[i], arr1[i + gap]);
        }

        if (gap >= n)
            i = 0;

        for (int idx1 = i, idx2 = (i + gap) % n; idx1 < n && idx2 < m; idx1++, idx2++)
        {
            if (arr1[idx1] > arr2[idx2])
                swap(arr1[idx1], arr2[idx2]);
        }

        for (int j = 0; j + gap < m; j++)
        {
            if (arr2[j] > arr2[j + gap])
                swap(arr2[j], arr2[j + gap]);
        }
    }
}

// Sort an array according to count of set bits
/*
Approach 1- O(nlogn),O(n)
Make another array of setbit count of each element, then sort it normally according to setbit count

Approach 2- Count Sort O(n), O(1)
As at max int can have 32 set bits so take a vector<vector<int>> of 32 size
at each index we have all element with setbit count=i
then just traverse in reverse and add them to result, to get the decreasing order of elements
*/
//Approach 2-
int setBitCount(int n)
{
    int count = 0;

    while (n != 0)
    {
        int rsbmask = (n & -n);
        n -= rsbmask;
        count++;
    }

    return count;
}
void setBitSort(vector<int> &arr)
{
    vector<vector<int>> setBits(32);
    for (int i = 0; i < arr.size(); i++)
    {
        int count = setBitCount(arr[i]);
        setBits[count].push_back(arr[i]);
    }

    for (int i = 31; i >= 0; i--)
    {
        if (setBits[i].size() > 0)
        {
            for (int ele : setBits[i])
                cout << ele << " ";
        }
    }

    cout << endl;
}

// Permute two arrays such that sum of every pair is greater or equal to K
int findAns(vector<long> &nums1, vector<long> &nums2, long k)
{
    sort(nums1.begin(), nums1.end());
    sort(nums2.begin(), nums2.end(), greater<long>());

    for (int i = 0; i < nums1.size(); i++)
        if (nums1[i] + nums2[i] < k)
            return 0;

    return 1;
}

// Minimum number of swaps required to sort an array
//Approach 1- Using Graphs (Not done)
int minSwaps(int arr[], int N)
{
    return 0;
}
//Approach 2- O(nlogn)
/*
copy the array to a temp array, and sort the temp array
then make a hashmap of ele:correct index
then for each index i keep swapping until correct element is at ith index
*/
int minSwaps(int arr[], int N)
{
    int temp[N];
    unordered_map<int, int> mp;
    for (int i = 0; i < N; i++)
    {
        temp[i] = arr[i];
    }

    sort(temp, temp + N);

    for (int i = 0; i < N; i++)
    {
        mp[temp[i]] = i;
    }

    int count = 0;
    for (int i = 0; i < N; i++)
    {
        while (arr[i] != temp[i])
        {
            swap(arr[i], arr[mp[arr[i]]]);
            count++;
        }
    }

    return count;
}

// 1363. Largest Multiple of Three
/*
Obviously, trying combinations of numbers won't work as we can have up to 10,000 numbers. Luckily, there is a handy divisibility test:
A number is divisible by 3 if the sum of all its digits is divisible by 3.
Observation 1: since the order does not matter, the largest number can be formed by adding digits from largest (9) to smallest (0), e.g. 9999966330000.
Therefore, we can just count the occurrences of each digit, and then generate the string.

Observation 2: we need to use all digits to form the maximum number.
If we sum all digits, and the modulo of 3 is not zero, we need to remove 1 (preferably) or 2 smallest digits.
If modulo 3 of the sum is 1, for example, we will try to remove 1, 4, or 7, if exists, or two of 2, 5, or 8.

More examples:
9965341 % 3 == 1; we remove 1 to get the largest number.
9952000 % 3 == 1; now we need to remove two digits, 2 and 5, as there is no 1, 4, or 7.
These observations yield the following algorithm.

Significance of sequence in m1[] = {1, 4, 7, 2, 5, 8}, m2[] = {2, 5, 8, 1, 4, 7}
If sum%3 == 1 then priority should be removing 1,4,7 before removing 2,5,8 as mentioned in m1.
If sum%3 == 2 then priority should be removing 2,5,8 before removing 1,4,7 as mentioned in m2.
*/
string largestMultipleOfThree(vector<int> &digits)
{
    if (digits.size() == 0)
        return "0";

    int m1[] = {1, 4, 7, 2, 5, 8}, m2[] = {2, 5, 8, 1, 4, 7};

    //make a freq array of all digits, and also get the total Sum
    int sum = 0;
    vector<int> freq(10, 0);
    for (auto d : digits)
    {
        freq[d]++;
        sum += d;
    }

    //remove the numbers until sum is divisible be 3
    while (sum % 3 != 0)
    {
        //if remainder is 1 -> remove numbers from m1
        if (sum % 3 == 1)
        {
            for (int num : m1)
            {
                if (freq[num] > 0)
                {
                    freq[num]--;
                    sum -= num;
                    break;
                }
            }
        }
        //if remainder is 2 -> remove numbers from m2
        else
        {
            for (int num : m2)
            {
                if (freq[num] > 0)
                {
                    freq[num]--;
                    sum -= num;
                    break;
                }
            }
        }
    }

    string res = "";
    for (int i = 9; i >= 0; i--)
        while (freq[i]-- > 0)
            res += i + '0';

    if (res[0] == '0')
        return "0";

    return res;
}

// Find pair with greatest product in array
/*
make a freq map of all elements in array and sort the array
start from end of array, for element at i, check from 0 till arr[j]<=sqrt(arr[i])
product=arr[i], num1=arr[j], num2=arr[i]/arr[j]
if mp[num2]>0 then we found the ans

Case:
1. Product, num1,num2 are all different, so mp[num2]>0 gives ans
2. num1==num2, then mp[num1]>1 gives ans
3. num1==1 , then num2==product, so mp[product]>1 gives ans

Eg. 1 93 1 5
here no ans

but in 1 93 93 5
here ans is 93, as there are two 93

*/
int largestProduct(vector<int> &arr)
{
    int n = arr.size();
    unordered_map<int, int> mp;

    for (int ele : arr)
        mp[ele]++;

    sort(arr.begin(), arr.end());

    for (int i = n - 1; i > 1; i--)
    {
        for (int j = 0; j < i && arr[j] <= sqrt(arr[i]); j++)
        {

            //if not divisible, skip
            if (arr[i] % arr[j] != 0)
                continue;

            int num1 = arr[j];
            int num2 = arr[i] / arr[j];

            //cases for num1==1
            if (num1 == 1 && mp[num2] > 1)
                return arr[i];
            else if (num1 == 1 && mp[num2] <= 1)
                continue;
            
            //cases for num2==1 && num2!=num1
            else if (num2 != num1 && mp[num2] > 0)
                return arr[i];
            else if (num2 == num1 && mp[num2] > 1)
                return arr[i];
        }
    }

    return -1;
}

void solve()
{
}
int main()
{
    solve();
    return 0;
}