#include <iostream>
#include <vector>
#include <list>
#include <algorithm>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <stack>

using namespace std;

// 287. Find the Duplicate Number
int findDuplicate(vector<int> &nums)
{
    int slow = nums[0];
    int fast = nums[0];

    do
    {
        slow = nums[slow];
        fast = nums[nums[fast]];
    } while (slow != fast);

    slow = nums[0];
    while (slow != fast)
    {
        slow = nums[slow];
        fast = nums[fast];
    }

    return slow;
}

//KCON (Codechef)
long long maxSubarraySum(int arr[], int n)
{
    long long maxSoFar = arr[0];
    long long currMax = 0;

    for (int i = 0; i < n; i++)
    {
        currMax += arr[i];
        maxSoFar = max(maxSoFar, currMax);
        if (currMax < 0)
            currMax = 0;
    }

    return maxSoFar;
}
void KCON()
{
    int t;
    cin >> t;
    while (t-- > 0)
    {
        int n, k;
        cin >> n >> k;
        int a[n], b[n * 2];
        for (long long i = 0; i < n; i++)
        {
            cin >> a[i];
            b[i] = b[i + n] = a[i];
        }

        //if k is 1 then simply print max subarray sum
        if (k == 1)
        {
            cout << maxSubarraySum(a, n) << endl;
            continue;
        }

        //calculate maxSubarray sum for A*2 i.e A put 2 times
        long long maxSum = maxSubarraySum(b, n * 2);
        long long maxPref = -1e9, maxSuff = -1e9, currPref = 0, currSuff = 0, totalSum = 0;

        //calculate prefix, suffix, total sum
        for (int i = 0; i < n; i++)
        {
            currPref += a[i];
            currSuff += a[n - i - 1];
            totalSum += a[i];
            maxPref = max(maxPref, currPref);
            maxSuff = max(maxSuff, currSuff);
        }

        //if totaSum>0, then only we should add all repetitions of A, otherwise it will just get smaller with every repetition
        if (totalSum > 0)
            maxSum = maxPref + totalSum * (k - 2) + maxSuff;

        cout << maxSum << endl;
    }
}

// Max Circular Subarray Sum
long long maxCircularSum(vector<int> &arr)
{
    long long minSum = arr[0], maxSum = arr[0], minSoFar = 0, maxSoFar = 0, totalSum = 0;

    //Kadannes algo to calculate minSubarray sum and maxSubarray sum
    for (int i = 0; i < arr.size(); i++)
    {
        totalSum += arr[i];
        minSoFar += arr[i];
        maxSoFar += arr[i];
        minSum = min(minSoFar, minSum);
        maxSum = max(maxSoFar, maxSum);
        if (minSoFar > 0)
            minSoFar = 0;
        if (maxSoFar < 0)
            maxSoFar = 0;
    }

    //if all -ve , then return the maxSubarray sum
    if (totalSum == minSum)
        return maxSum;

    //else return the max of maxSubarray or the circular array
    return max(maxSum, totalSum - minSum);
}

// Subarray with given sum
/*
Use 2 pointer approach
i=0,j=0
if currSum<tar -> j++
if currSum>tar -> i++
*/
vector<int> findSubarray(vector<int> &arr, int tar)
{
    //use two pointer approach
    int i = 0, j = 0, currSum = 0;
    while (j < arr.size())
    {
        currSum += arr[j];
        if (currSum > tar)
        {
            while (i <= j && currSum > tar)
                currSum -= arr[i++];
        }
        if (currSum == tar)
            break;
        j++;
    }
    //1 based indexing, so i+1,j+1
    return {i + 1, j + 1};
}
void findSubarray()
{
    int t;
    cin >> t;
    while (t-- > 0)
    {
        int n, tar;
        cin >> n >> tar;

        vector<int> a(n);
        for (int i = 0; i < n; i++)
            cin >> a[i];

        vector<int> ans = findSubarray(a, tar);
        if (ans[1] > n)
        {
            cout << -1 << endl;
            continue;
        }
        cout << ans[0] << " " << ans[1] << endl;
    }
}

// 724. Find Pivot Index (Equilibrium Point)
/*
Note- The Eq point itself is not part of the left or right Sum 
use leftSum, rightSum
find total Sum of array=rsum
start iterating from begin, 
lsum=sum+a[i-1]
rsum=rsum-a[i]
*/
int pivotIndex(vector<int> &a)
{
    int n = a.size();
    if (n == 0)
        return -1;
    int lsum = 0, rsum = 0; //left = 0 , right = total
    for (int i = 0; i < n; i++)
        rsum += a[i];

    if (rsum - a[0] == 0)
        return 0;

    rsum -= a[0];

    for (int i = 1; i < n; i++)
    {
        rsum -= a[i];
        lsum += a[i - 1];
        if (lsum == rsum)
        {
            return i;
        }
    }

    return -1;
}

// Convert array into Zig-Zag fashion
/*
keep a flag to check required condition nextEle > currEle or nextEle < currEle
if condition is false then swap the two elements
*/
void zigzag(vector<int> &a)
{
    int n = a.size();
    int flag = 0; //next element should be 0->inc , 1->dec

    for (int i = 0; i < n - 1; i++)
    {
        if (flag == 0 && a[i + 1] < a[i])
        {
            swap(a[i + 1], a[i]);
        }
        else if (flag == 1 && a[i + 1] > a[i])
        {
            swap(a[i + 1], a[i]);
        }
        flag ^= 1;
    }
}

// Find Pair Given Difference
/*
Sort the array
Use 2 Pointer Approach, i=0,j=1
currDiff=arr[j]-arr[i]
while(cond..){
if(currDiff>tar) i++
if(currDiff<tar) j++
} 
*/
int diffPair(vector<int> &a, int d)
{
    int n = a.size();
    sort(a.begin(), a.end());

    int diff, i = 0, j = 1;
    while (j < n && i < n)
    {
        diff = a[j] - a[i];
        //diff will always be positive, so i<=j is true always
        if (diff > d)
        {
            i++;
        }
        else if (diff < d)
        {
            j++;
        }
        else if (i != j && diff == d)
        {
            return 1;
        }
    }

    return -1;
}

// Chocolate Distribution Problem (no submission option available on G4G)
int minDiff(vector<int> &packets, int children)
{
    sort(packets.begin(), packets.end());

    int minDiff = 1e8;

    // Find the subarray of size m such that
    // difference between last (maximum in case
    // of sorted) and first (minimum in case of
    // sorted) elements of subarray is minimum.
    for (int i = 0; i + children - 1 < packets.size(); i++)
    {
        minDiff = min(minDiff, packets[i + children - 1] - packets[i]);
    }

    return minDiff;
}

// Minimum Number of Platforms Required for a Railway/Bus Station
/*
sort both arrival and departure times
then do it like merging two sorted arrays
*/
int getStations(vector<int> &arrival, vector<int> &departure)
{
    int n = arrival.size();
    int currStations = 0, minStations = 0;
    sort(arrival.begin(), arrival.end());
    sort(departure.begin(), departure.end());

    int i = 0, j = 0;
    //we only need to check for arrival array for size check as all trains will arrive first,
    //so arrival array will finish first always
    while (i < n && j < n)
    {
        //check for <= as arrival and departure times can be same as well and we need a seperate stations in this case
        //train arrives -> currStations++
        if (arrival[i] <= departure[j])
        {
            i++;
            currStations++;
        }
        //train departs -> currStations--
        else
        {
            j++;
            currStations--;
        }
        //update the max stations at a time
        minStations = max(currStations, minStations);
    }

    return minStations;
}
void getStations()
{
    int t;
    cin >> t;
    while (t-- > 0)
    {
        int n;
        cin >> n;
        vector<int> arrival(n), departure(n);
        for (int i = 0; i < n; i++)
            cin >> arrival[i];
        for (int i = 0; i < n; i++)
            cin >> departure[i];

        cout << getStations(arrival, departure) << endl;
    }
}

// 373. Find K Pairs with Smallest Sums
/*
using min heap
the approach is similar to merging K sorted Lists using Priority Queue
*/
struct compareSum
{
    bool operator()(const vector<int> &v1, const vector<int> &v2)
    {
        return v1[0] + v1[1] > v2[0] + v2[1];
    }
};
vector<vector<int>> kSmallestPairs(vector<int> &nums1, vector<int> &nums2, int k)
{
    if (nums1.size() == 0 || nums2.size() == 0 || k == 0)
        return {};

    vector<vector<int>> ans;

    priority_queue<vector<int>, vector<vector<int>>, compareSum> pq; // {nums1[i], nums2[j], j}

    //form pairs with the first element of num2 with all elements in num1
    //so we add all (0,0), (1,0), (2,0)........(n-1,0) into PQ
    for (int i = 0; i < nums1.size() && i < k; i++)
    {
        pq.push({nums1[i], nums2[0], 0});
    }

    //for each pair removed from PQ push the pair with next largest sum into it
    //we have a choice to add either (i,j+1) or (i+1,j), but we only add (i,j+1)
    //as the (i+1,j) we already be in the PQ or not qualified.
    //PQ is of size K
    while (k-- > 0 && pq.size() != 0)
    {
        vector<int> pr = pq.top();
        pq.pop();

        //add the top of PQ into the answer
        ans.push_back({pr[0], pr[1]});

        if (pr[2] == nums2.size() - 1)
            continue;

        //next largest sum will be given by the next element of nums2 from the element removed from PQ
        pq.push({pr[0], nums2[pr[2] + 1], pr[2] + 1});
    }

    return ans;
}

// 26. Remove Duplicates from Sorted Array
/*
Use two pointer approach O(n)
if(a[i]==a[j]) then only j++
if(!=) a[i+1]=a[j]; 
        i++; 
        j++;
*/
int removeDuplicates(vector<int> &nums)
{
    if (nums.size() == 0)
        return 0;
    int i = 0, j = 1;
    while (j < nums.size())
    {
        if (nums[i] != nums[j])
        {
            nums[i + 1] = nums[j];
            i++;
        }
        j++;
    }

    return i + 1;
}

// 153. Find Minimum in Rotated Sorted Array
/*
use binary search
if arr[mid]<arr[hi], then smallest element i.e. pivot lies in left region i.e. region<=mid
else it lies in right region i.e. in region>mid
*/
int findMin(vector<int> &nums)
{
    int lo = 0, hi = nums.size() - 1;
    long mid;
    while (lo < hi)
    {
        mid = (lo + hi) / 2;
        if (nums[mid] < nums[hi])
            hi = mid;
        else
            lo = mid + 1;
    }

    return nums[lo];
}

// 154. Find Minimum in Rotated Sorted Array II (Duplicates allowed)
/*
if arr[mid]==arr[hi], then we have no choice but to do a linear search in whole array,
as we cannot decide which part pivot lies in.
So in that case we just reduce the hi by 1
*/
int findMin(vector<int> &nums)
{
    //Worst case: not rotated array, O(n)
    int lo = 0, hi = nums.size() - 1;
    long mid;
    while (lo < hi)
    {
        mid = (lo + hi) / 2;
        if (nums[mid] < nums[hi])
            hi = mid;
        else if (nums[mid] > nums[hi])
            lo = mid + 1;
        else
        {
            /*if (i-1)th is greater that means (i)th will be the pivot
            eg- 1 1 1 1 2 1 1
            here at i=5 , a[i-1]>a[i] and i is the pivot index
            */
            if (nums[hi - 1] > nums[hi])
                return nums[hi];
            hi--;
        }
    }
    return nums[lo];
}

// 33. Search in Rotated Sorted Array
int search_01(vector<int> &nums, int target)
{
    int lo = 0, hi = nums.size() - 1;
    long mid;
    //find smallest element(pivot)
    while (lo < hi)
    {
        mid = (lo + hi) / 2;
        if (nums[mid] < nums[hi])
            hi = mid;
        else
            lo = mid + 1;
    }

    //find which region target lies in
    int pivot = lo;
    lo = 0, hi = nums.size() - 1;
    if (target <= nums[hi])
    { //right half of pivot
        lo = pivot;
    }
    else
    { //left half of pivot
        hi = pivot - 1;
    }

    //normal binary search in that region
    while (lo <= hi)
    {
        mid = (lo + hi) / 2;
        if (target < nums[mid])
            hi = mid - 1;
        else if (target > nums[mid])
            lo = mid + 1;
        else
            return mid;
    }
    return -1;
}
int search_02(vector<int> &nums, int target)
{
    int lo = 0, hi = nums.size() - 1;
    long mid;
    //find smallest element(pivot)
    while (lo < hi)
    {
        mid = (lo + hi) / 2;
        if (nums[mid] < nums[hi])
            hi = mid;
        else
            lo = mid + 1;
    }

    //binary search accounting for rotation
    int rot = lo, n = nums.size();
    lo = 0;
    hi = n - 1;
    while (lo <= hi)
    {
        int mid = (lo + hi) / 2;
        int realmid = (mid + rot) % n;
        if (nums[realmid] == target)
            return realmid;
        if (nums[realmid] < target)
            lo = mid + 1;
        else
            hi = mid - 1;
    }
    return -1;
}

// Given a sorted and rotated array, find if there is a pair with a given sum(Not available for submission)
bool pairSum(vector<int> &nums, int target)
{
    //find pivot
    int lo = 0, hi = nums.size() - 1;
    long mid;
    while (lo < hi)
    {
        mid = (lo + hi) / 2;
        if (nums[mid] < nums[hi])
            hi = mid;
        else
            lo = mid + 1;
    }

    //use 2 pointer method(meet in the middle) using mod to keep index in range
    int pivot = lo, n = nums.size();
    int i = lo, j = pivot - 1;
    //i goes pivot -> n-1 , j goes pivot-1 -> 0
    while (i != j)
    {
        if (nums[i] + nums[j] < target)
        {
            //to keep i++ in range
            i = (i + 1) % n;
        }
        else if (nums[i] + nums[j] > target)
        {
            //to keep j-- in range
            j = (n + j - 1) % n;
        }
        else
            return true;
    }

    return false;
}

// 525. Contiguous Array
/*
we do count-- for 0s, and count++ for 1s
keep a map for count:index values
two index where count is equal, will have equal number of 0s and 1s, so at each index we just check 
if this count value is present in map and update the max subarray length
*/
int findMaxLength(vector<int> &nums)
{
    int count = 0, maxlen = 0;
    unordered_map<int, int> mp; //{count:index}
    mp[0] = -1;
    for (int i = 0; i < nums.size(); i++)
    {
        if (nums[i] == 0)
            count--;
        else
            count++;
        if (mp.find(count) != mp.end())
            maxlen = max(maxlen, i - mp[count]);
        else
            mp[count] = i;
    }
    return maxlen;
}

// 121. Best Time to Buy and Sell Stock
/*
Just find the max diff in array
Keep a min so far for each element, and update the maxProfit with the max diff i.e arr[i]-minSoFar
*/
int maxProfit(vector<int> &prices)
{
    if (prices.size() == 0)
        return 0;
    int maxp = 0, minSoFar = prices[0];
    for (int i = 0; i < prices.size(); i++)
    {
        minSoFar = min(minSoFar, prices[i]);
        maxp = max(maxp, prices[i] - minSoFar);
    }

    return maxp;
}

// 122. Best Time to Buy and Sell Stock II
/*
Approach 1-
while the price keeps increasing keep adding it to profit,
when it decreases, dont include it
*/
int maxProfit(vector<int> &prices)
{
    if (prices.size() == 0)
        return 0;
    int totalProfit = 0;
    for (int i = 0; i < prices.size() - 1; i++)
    {
        if (prices[i + 1] > prices[i])
            totalProfit += (prices[i + 1] - prices[i]);
    }
    return totalProfit;
}
/*
Approach 2-
keep a currMin (for the current rise in price),
the moment the price drops, you sell the stock and add arr[i]-currMin into the total profit
*/
int maxProfit(vector<int> &prices)
{
    if (prices.size() == 0)
        return 0;
    int currMin = prices[0], totalProfit = 0;
    for (int i = 0; i < prices.size() - 1; i++)
    {
        if (prices[i] > prices[i + 1])
        {
            totalProfit += prices[i] - currMin;
            currMin = prices[i + 1];
        }
    }
    totalProfit += prices[prices.size() - 1] - currMin;

    return totalProfit;
}

// 123. Best Time to Buy and Sell Stock III
/*Approach 1-
Make a prefix and suffix array
Prefix- Stock is sold on this day
Suffix- Stock is bought on this day

Then the max profit the max sum of prefix and suffix arrays
*/
int maxProfit(vector<int> &prices)
{
    if (prices.size() == 0)
        return 0;

    int n = prices.size();
    vector<int> prefix(n), suffix(n);

    //we can maintain the max Profit so far seperatly,
    /*
    int minBuying = prices[0], maxProfitSoFar = 0;
    for (int i = 0; i < n; i++)
    {
        minBuying = min(prices[i], minBuying);
        maxProfitSoFar = max(maxProfitSoFar, prices[i] - minBuying);
        prefix[i] = maxProfitSoFar;
    }

    int maxSelling = prices[n - 1];
    maxProfitSoFar = 0;
    for (int i = n - 1; i >= 0; i--)
    {
        maxSelling = max(maxSelling, prices[i]);
        maxProfitSoFar = max(maxProfitSoFar, maxSelling - prices[i]);
        suffix[i] = maxProfitSoFar;
    }
    */

    //or we can just use the value in i-1 index, as it already has the max value upto that point
    int minBuying = prices[0];
    prefix[0] = 0;
    for (int i = 1; i < n; i++)
    {
        minBuying = min(prices[i], minBuying);
        prefix[i] = max(prefix[i - 1], prices[i] - minBuying);
    }

    int maxSelling = prices[n - 1];
    suffix[n - 1] = 0;
    for (int i = n - 2; i >= 0; i--)
    {
        maxSelling = max(maxSelling, prices[i]);
        suffix[i] = max(suffix[i + 1], maxSelling - prices[i]);
    }

    int maxp = 0;
    for (int i = 0; i < n; i++)
    {
        maxp = max(maxp, prefix[i] + suffix[i]);
    }

    return maxp;
}

//Approach 2 (faster)-
int maxProfit(vector<int> &prices)
{
    int sell1 = 0, sell2 = 0, buy1 = 1e8, buy2 = 1e8;
    for (int i = 0; i < prices.size(); i++)
    {
        buy1 = min(buy1, prices[i]);
        sell1 = max(sell1, prices[i] - buy1);
        buy2 = min(buy2, prices[i] - sell1);
        sell2 = max(sell2, prices[i] - buy2);
    }
    return sell2;
}

// 188. Best Time to Buy and Sell Stock IV (Not Complete)
int maxp; //max profit
//shares==1 -> have an extra share, so can only sell
//shares==0 -> dont have a share so can ony buy
//Recursion
//void type
void maxProfit_rec1(int shares, int profit, int K, int idx, vector<int> &prices)
{
    maxp = max(maxp, profit);
    if (K == 0 || idx == prices.size())
    {
        return;
    }

    //buy
    if (shares == 0)
    {
        maxProfit_rec1(1, profit - prices[idx], K, idx + 1, prices);
    }

    //sell
    else
    {
        maxProfit_rec1(0, profit + prices[idx], K - 1, idx + 1, prices);
    }

    //do nothing
    maxProfit_rec1(shares, profit, K, idx + 1, prices);
}

//return type
int maxProfit_rec2(int shares, int profit, int K, int idx, vector<int> &prices)
{

    if (K == 0 || idx == prices.size())
    {
        return profit;
    }

    int maxp = 0;
    //buy
    if (shares == 0)
    {
        maxp = max(maxp, maxProfit_rec2(1, profit - prices[idx], K, idx + 1, prices));
    }

    //sell
    else
    {
        maxp = max(maxp, maxProfit_rec2(0, profit + prices[idx], K - 1, idx + 1, prices));
    }

    //do nothing
    maxp = max(maxp, maxProfit_rec2(shares, profit, K, idx + 1, prices));

    return maxp;
}
int maxProfit(int k, vector<int> &prices)
{
    return maxProfit_rec2(0, 0, k, 0, prices);
}

// 714. Best Time to Buy and Sell Stock with Transaction Fee
int maxProfit(vector<int> &prices, int fee)
{
    if (prices.size() == 0)
        return 0;

    int maxAfterBuy = -prices[0]; //current cash in hand after buying
    int maxAfterSell = 0;         //current cash in hand after selling
    for (int i = 1; i < prices.size(); i++)
    {
        maxAfterBuy = max(maxAfterBuy, maxAfterSell - prices[i]);
        maxAfterSell = max(maxAfterSell, maxAfterBuy + prices[i] - fee);
    }

    return maxAfterSell;
}

// 309. Best Time to Buy and Sell Stock with Cooldown
int maxProfit(vector<int> &prices)
{
    if (prices.size() == 0)
        return 0;

    int maxAfterBuy = -prices[0]; //current cash in hand after buying
    int maxAfterSell = 0;         //current cash in hand after selling
    int prevMaxSell = 0;          //max cash after cooldown(stores the max sell just before cooldown)
    for (int i = 1; i < prices.size(); i++)
    {
        int prevMaxBuy = maxAfterBuy;
        maxAfterBuy = max(maxAfterBuy, prevMaxSell - prices[i]);
        prevMaxSell = maxAfterSell;
        maxAfterSell = max(maxAfterSell, prevMaxBuy + prices[i]);
    }

    return maxAfterSell;
}

// Maximum Difference (given that second element is greater than first element)
int maxDiff(vector<int> &prices)
{
    if (prices.size() == 0)
        return 0;
    int maxp = 0, minSoFar = prices[0];
    for (int i = 0; i < prices.size(); i++)
    {
        minSoFar = min(minSoFar, prices[i]);
        maxp = max(maxp, prices[i] - minSoFar);
    }

    //if no secondEle > firstEle then return -1
    return maxp == 0 ? -1 : maxp;
}
int maxDiff()
{
    int t;
    cin >> t;
    while (t-- > 0)
    {
        int n;
        cin >> n;
        vector<int> arr(n);
        for (int i = 0; i < n; i++)
            cin >> arr[i];
        cout << maxDiff(arr) << endl;
    }
}

// 56. Merge Intervals
vector<vector<int>> merge(vector<vector<int>> &intervals)
{
    if (intervals.size() == 0)
        return {};

    vector<vector<int>> res;
    //sort the array
    sort(intervals.begin(), intervals.end());
    res.push_back(intervals[0]);
    int j = 0;

    //for each interval check if it overlaps with last one, and add it accordingly
    for (int i = 1; i < intervals.size(); i++)
    {
        //if it overlaps, merge the two intervals
        if (intervals[i][0] <= res[j][1])
        {
            res[j][0] = min(res[j][0], intervals[i][0]); //start will min of both starts
            res[j][1] = max(res[j][1], intervals[i][1]); //end will be max of both ends
        }
        //if it does not overlap, just add it
        else
        {
            res.push_back(intervals[i]);
            j++;
        }
    }

    return res;
}

// 57. Insert Interval
vector<vector<int>> insert(vector<vector<int>> &intervals, vector<int> &newInterval)
{
    if (intervals.size() == 0)
        return {newInterval};
    vector<vector<int>> res;
    int idx = 0;
    //push all intervals less than newInterval into res
    while (idx < intervals.size() && intervals[idx][1] < newInterval[0])
        res.push_back(intervals[idx++]);

    //find and push the merged Interval into res
    int mergedStart = newInterval[0], mergedEnd = newInterval[1];
    while (idx < intervals.size() && newInterval[1] >= intervals[idx][0])
    {
        mergedStart = min(mergedStart, intervals[idx][0]);
        mergedEnd = max(mergedEnd, intervals[idx][1]);
        idx++;
    }
    res.push_back({mergedStart, mergedEnd});

    //push the remaining intervals into res
    while (idx < intervals.size())
        res.push_back(intervals[idx++]);

    return res;
}

// 75. Sort Colors (3 Way Partition)
void sortColors(vector<int> &nums)
{
    int lt = 0, gt = nums.size() - 1, i = 0;
    while (i <= gt)
    {
        //left region: i++, lt++ (after swapping element will definitely belong to left region)
        if (nums[i] < 1)
        {
            swap(nums[lt++], nums[i++]);
        }
        //right region: only gt-- (because after swapping we could get an element that may belong to left or middle region)
        else if (nums[i] > 1)
        {
            swap(nums[gt--], nums[i]);
        }
        //middle region: i++
        else
            i++;
    }
}

// Three way partitioning
vector<int> threeWayPartition(vector<int> nums, int a, int b)
{
    int lt = 0, gt = nums.size() - 1, i = 0;
    while (i <= gt)
    {
        if (nums[i] < a)
        {
            swap(nums[lt++], nums[i++]);
        }
        else if (nums[i] > b)
        {
            swap(nums[gt--], nums[i]);
        }
        else
            i++;
    }

    return nums;
}

// 912. Sort an Array
void mergeSort(vector<int> &arr, int si, int ei)
{
    if (si == ei)
        return;

    long mid = (si + ei) / 2;
    mergeSort(arr, si, mid);
    mergeSort(arr, mid + 1, ei);

    merge(arr, si, mid, ei);
}
void merge(vector<int> &arr, int si, int mid, int ei)
{
    //make a temp res array, to store sorted list, then copy it into original array
    //size of res= ei-si+1, but it will be initialized by zero if you provide size,
    // so will have to use an index instead of push_back
    vector<int> res;
    int i = si, j = mid + 1;     //start index of both halfs
    int n = mid + 1, m = ei + 1; //end index of both halfs

    while (i < n && j < m)
    {
        if (arr[i] < arr[j])
            res.push_back(arr[i++]);
        else
            res.push_back(arr[j++]);
    }

    while (i < n)
        res.push_back(arr[i++]);
    while (j < m)
        res.push_back(arr[j++]);

    for (int k = si; k < ei + 1; k++)
        arr[k] = res[k - si];
}
vector<int> sortArray(vector<int> &nums)
{
    mergeSort(nums, 0, nums.size() - 1);
    return nums;
}

// Inversion of array
long merge_(vector<int> &arr, int si, int mid, int ei)
{
    vector<int> res;
    int i = si, j = mid + 1;
    int n = mid + 1, m = ei + 1;
    long count = 0;

    while (i < n && j < m)
    {
        if (arr[i] <= arr[j])
            res.push_back(arr[i++]);
        else
        {
            res.push_back(arr[j++]);
            //when the element in first subarray > element in second subarray, then there is an inversion
            count += n - i;
        }
    }

    while (i < n)
        res.push_back(arr[i++]);
    while (j < m)
        res.push_back(arr[j++]);

    for (int k = si; k < ei + 1; k++)
        arr[k] = res[k - si];

    return count;
}
long mergeSort_(vector<int> &arr, int si, int ei)
{
    if (si == ei)
        return 0;

    long mid = (si + ei) / 2;
    long lcount = mergeSort_(arr, si, mid);
    long rcount = mergeSort_(arr, mid + 1, ei);

    long myCount = merge_(arr, si, mid, ei);

    return lcount + rcount + myCount;
}
long countInversions(vector<int> &arr)
{
    return mergeSort_(arr, 0, arr.size() - 1);
}

// 775. Global and Local Inversions
/*
Local inverions will occur between consecutive elements only, 
so if there is a single inversion not between consecutive elements, 
then there will be more global than local inversions.

All local inversions are global inversions.
If the number of global inversions is equal to the number of local inversions,
it means that all global inversions in permutations are local inversions.
It also means that we can not find A[i] > A[j] with i+2<=j.
meaning for global==local to be true, inversions can only be between i,i+1
*/
bool isIdealPermutation(vector<int> &A)
{
    if (A.size() < 2)
        return true;
    int currMax = 0;
    for (int i = 0; i < A.size() - 2; i++)
    {
        currMax = max(currMax, A[i]);
        //if at any point current max element is greater than the i+2 position element, then global>local
        if (currMax > A[i + 2])
            return false;
    }
    return true;
}

// 560. Subarray Sum Equals K
int subarraySum(vector<int> &nums, int k)
{
    unordered_map<int, int> mp; //prefixSum : count (how many times it has occured till now)
    int currSum = 0, count = 0;

    for (int i = 0; i < nums.size(); i++)
    {
        currSum += nums[i];
        if (currSum == k)
            count++;
        if (mp.find(currSum - k) != mp.end())
            count += mp[currSum - k];
        mp[currSum]++;
    }

    return count;
}

// 152. Maximum Product Subarray
int maxProduct(vector<int> &nums)
{
    if (nums.size() == 0)
        return 0;

    int maxpos = nums[0], minneg = nums[0], omax = nums[0]; //omax=overall max

    for (int i = 1; i < nums.size(); i++)
    {
        if (nums[i] < 0)
            swap(minneg, maxpos);

        maxpos = max(nums[i], maxpos * nums[i]);
        minneg = min(nums[i], minneg * nums[i]);

        omax = max(omax, maxpos);
    }

    return omax;
}

// Minimize the heights
/*
Just remember it (not sure about logic)
*/
int getMinDiff(int arr[], int n, int k)
{
    //sort the heights
    sort(arr, arr + n);

    //find min difference in before changing with +K or -K
    int minDiff = arr[n - 1] - arr[0];

    int minEle = arr[0] + k, maxEle = arr[n - 1] - k;
    if (maxEle < minEle)
        swap(maxEle, minEle);

    //for the whole array, find the min and max you can get
    for (int i = 1; i < n - 1; i++)
    {
        int currLargest = arr[i] + k;
        int currSmallest = arr[i] - k;

        //if curr Ele does not become both currMax and currMin after +k and -k then skip it
        if (currLargest < maxEle || currSmallest > minEle)
            continue;

        //choose whether it gives minDiff after +k or -k
        if (maxEle - currSmallest <= currLargest - minEle)
            minEle = currSmallest;
        else
            maxEle = currLargest;
    }

    //return min of the original mindiff and after +K and -K
    return min(minDiff, maxEle - minEle);
}

// Minimum swaps and K together
/*
Using Sliding window
Find count of all elements which are less than or equals to ‘k’. Let’s say the count is ‘cnt’
Using two pointer technique for window of length ‘cnt’, each time keep track of how many elements in this range are greater than ‘k’.
Let’s say the total count is ‘bad’.
Repeat step 2, for every window of length ‘cnt’ and take minimum of count ‘bad’ among them. This will be the final answer.
*/
long minSwaps(vector<int> &arr, int k)
{
    int n = arr.size();
    long count = 0, currCount = 0, minCount;

    //find total count
    for (int i = 0; i < n; i++)
        if (arr[i] <= k)
            count++;

    //find count for first window
    for (int i = 0; i < count; i++)
    {
        if (arr[i] > k)
            currCount++;
    }

    //after that, just for each window, just check the element that was removed from it (i-count), and added to it(current i)
    minCount = currCount;
    for (int i = count; i < n; i++)
    {
        if (arr[i - count] > k)
            currCount--;
        if (arr[i] > k)
            currCount++;

        minCount = min(minCount, currCount);
    }

    return minCount;
}

// Maximum sum of i*arr[i] among all rotations of a given array
/*
We can take rotations in either direction, I am taking anti clockwise(towards left)
Brute Force- O(n^2)
Calculate all sums for all rotations of array

Efficient-
Calculate sum of all arr[i]*i for original array, total Sum of all elements
after each rotation the change in sum will be due to two things:;
    the 0th element goes to n-1, arr[i]*0 = 0 becomes arr[i]*(n-1)>0 , so we do rotatedSum+=arr[i]*(n-1)
    every other element goes from i to i-1, so we subtract the total Sum from rotatedSum except for the arr[i] element 
so rotatedSum=previousRotatedSum + arr[i]*(n-1) - (totalSum-arr[i])
*/
int max_sum(int arr[], int n)
{
    long long totalSum = 0, rsum = 0, maxSum = -1e8;
    //calculate initial totalSum, rotatedSum
    for (int i = 0; i < n; i++)
    {
        totalSum += arr[i];
        rsum += arr[i] * i;
    }

    //for each rotation, find the new rotatedSum, by using the formula above
    for (int i = 0; i < n; i++)
    {
        rsum += arr[i] * (n - 1) - (totalSum - arr[i]);
        maxSum = max(maxSum, rsum);
    }

    return maxSum;
}

// 41. First Missing Positive
/*
Approach-
Put every arr[i]>0 at i+1 index
then iterate through array, and wherever arr[i]!=i+1, i+1 is the missing number
if none is missing then n+1 is the missng number
*/
int firstMissingPositive(vector<int> &nums)
{
    int n = nums.size();

    for (int i = 0; i < n; i++)
    {
        //if the number is in range of array, and not already at the right position
        //swap it with the number at that position
        while (nums[i] > 0 && nums[i] < n + 1 && nums[nums[i] - 1] != nums[i])
            swap(nums[i], nums[nums[i] - 1]);
    }

    for (int i = 0; i < n; i++)
    {
        if (nums[i] != i + 1)
            return i + 1;
    }

    return n + 1;
}

// 283. Move Zeroes (Move all zeroes to end of array)
//Approach 1-
void moveZeroes_01(vector<int> &nums)
{
    int idx = 0;
    for (int i = 0; i < nums.size(); i++)
    {
        if (nums[i] != 0)
            nums[idx++] = nums[i];
    }
    for (int i = idx; i < nums.size(); i++)
        nums[i] = 0;
}
//Approach 2- (Optimal)
void moveZeroes_02(vector<int> &nums)
{
    /*
    instead of going through array again and changing to 0,
    we just do it during first iteration, by swapping the non zero element with last position of non zero element
    */
    int idx = 0;
    for (int i = 0; i < nums.size(); i++)
    {
        if (nums[i] != 0)
            swap(nums[idx++], nums[i]);
    }
}

// Largest Sum Subarray of Size at least K
long long largestSum(vector<int> &arr, int k)
{
    int n = arr.size();
    vector<long long> prefixMaxSum(n);

    //make prefix array of all max contiguos sum till that index
    long long maxSoFar = arr[0], currMax = 0;
    for (int i = 0; i < n; i++)
    {
        currMax += arr[i];
        maxSoFar = max(maxSoFar, currMax);
        if (currMax < 0)
            currMax = 0;
        prefixMaxSum[i] = currMax;
    }

    //for each window (i to j),find window sum, and update overall max with max(window sum, windowSum + maxPrefixSum)
    // for the last index(i-1)
    long long currSum = 0, maxSum;
    for (int i = 0; i < k; i++)
        currSum += arr[i];
    maxSum = currSum;
    for (int i = k; i < n; i++)
    {
        //currSum = window sum
        currSum += arr[i] - arr[i - k];
        //prefixMaxSum[i-k] = max Subarray sum for the last index just before the current window
        maxSum = max(maxSum, max(currSum, currSum + prefixMaxSum[i - k]));
    }

    return maxSum;
}

// 179. Largest Number (Arrange given numbers to form the biggest number)
static int myCompare(string x, string y)
{
    /*
    For example, let X and Y be 542 and 60. 
    To compare X and Y, we compare 54260 and 60542. 
    Since 60542 is greater than 54260, we put Y first.
    */
    string xy = x + y;
    string yx = y + x;
    return xy.compare(yx) > 0 ? 1 : 0;
}
string largestNumber(vector<int> &nums)
{
    //convert all numbers to string
    vector<string> sortedNums(nums.size());
    for (int i = 0; i < nums.size(); i++)
        sortedNums[i] = to_string(nums[i]);

    //sort the string array using custom comparator
    sort(sortedNums.begin(), sortedNums.end(), myCompare);

    //edge case: all are 0 in array, so greatest in sorted i.e arr[0]='0'
    if (sortedNums[0] == "0")
        return "0";

    //form the greatest number
    string maxNum = "";
    for (int i = 0; i < sortedNums.size(); i++)
        maxNum += sortedNums[i];

    return maxNum;
}

// Form minimum number from given sequence (G4G, only 1-9 allowed)
/*Brute Force- Backtracking(TLE)
At every step we choice between 1-9, depending on already visited numbers
1.Check current pattern
2.if 'D' : currNum=minimum available
else: currNum=min available whuich is greater than last Available number
3.if the required number is invalid i.e. <1 || >9, backtrack
*/
int ans;
bool minNum(string &pattern, vector<bool> &vis, int idx, int currNum)
{
    if (idx == pattern.size())
    {
        ans = currNum;
        return true;
    }
    int ln = currNum % 10;

    if (pattern[idx] == 'D')
    {
        for (int i = 0; i < 9; i++)
        {
            if (!vis[i] && i + 1 < ln)
            {
                int num = i + 1;
                vis[num - 1] = true;
                if (minNum(pattern, vis, idx + 1, currNum * 10 + num))
                    return true;
                vis[num - 1] = false;
            }
        }
    }
    else
    {
        for (int i = 0; i < 9; i++)
        {
            if (!vis[i] && i + 1 > ln)
            {
                int num = i + 1;
                vis[num - 1] = true;
                if (minNum(pattern, vis, idx + 1, currNum * 10 + num))
                    return true;
                vis[num - 1] = false;
            }
        }
    }

    return false;
}
int minNum(string &pattern)
{
    vector<bool> vis(9, false);
    for (int i = 1; i <= 9; i++)
    {
        vis[i - 1] = true;
        if (minNum(pattern, vis, 0, i))
            return ans;
        vis[i - 1] = false;
    }

    return -1;
}

//884. Find Permutation(Lintcode) (Form minimum number from given sequence)
/*Approach 1(Use Stack)- Time:O(n), Space: O(n)
If it's all just I, then the answer is the numbers in ascending order.
And if there are streaks of D, then just reverse the number streak under each

Eg-
 I I I I I I I I
1 2 3 4 5 6 7 8 9 

 I D D I D D D I
1 4 3 2 8 7 6 5 9

*/
vector<int> findPermutation(string &s)
{
    stack<int> st;
    vector<int> res;
    st.push(-1);
    int lastMin = 1;

    //base case
    if (s[0] == 'I')
    {
        res.push_back(1);
        lastMin++;
    }
    else
    {
        st.push(lastMin++);
    }

    for (int i = 0; i < s.size(); i++)
    {
        //if D, push into stack the lastMin
        if (s[i] == 'D')
        {
            st.push(lastMin++);
        }

        //if I, pop everything from stack and push into res, then push lastMin into stack
        else
        {
            while (st.top() != -1)
            {
                res.push_back(st.top());
                st.pop();
            }
            st.push(lastMin++);
        }
    }

    while (st.top() != -1)
    {
        res.push_back(st.top());
        st.pop();
    }

    return res;
}
//Approach 2: Time: O(n), Space: O(1), same as previous, just emulate stack using 2 pointers
vector<int> findPermutation(string &s)
{
    vector<int> res;
    int lastMin = 1, lastIdx = -1; //lastIdx= index of last 'I'

    //base case
    if (s[0] == 'I')
    {
        lastIdx = 0;
        res.push_back(lastMin);
    }
    else
        lastIdx = -1;

    for (int i = 0; i < s.size(); i++)
    {
        //if D, push into stack the lastMin
        if (s[i] == 'D')
            lastMin++;

        //if I, pop everything from stack and push into res, then push lastMin into stack
        else
        {
            int num = lastMin;
            for (int j = lastIdx; j < i; j++)
                res.push_back(num--);
            lastMin++;
            lastIdx = i;
        }
    }

    //push whatever remains in stack into res
    int num = lastMin;
    for (int j = lastIdx; j < s.size(); j++)
        res.push_back(num--);

    return res;
}
//Approach 2- shorter code(removed the if(s[i]=='D) condition as we didnt do anything in it, and other things)
vector<int> findPermutation(string &s)
{
    vector<int> res;
    int lastMin = 1, lastIdx = -1; //lastIdx= index of last 'I'

    //base case
    if (s[0] == 'I')
    {
        lastIdx = 0;
        res.push_back(lastMin);
    }
    else
        lastIdx = -1;

    for (int i = 0; i < s.size(); i++)
    {
        if (s[i] == 'I')
        {
            int num = lastMin;
            for (int j = lastIdx; j < i; j++)
                res.push_back(num--);
            lastIdx = i;
        }
        lastMin++;
    }

    for (int j = lastIdx; j < s.size(); j++)
        res.push_back(lastMin--);

    return res;
}

// Find the smallest positive integer value that cannot be represented as sum of any subset of a given array
/*
Q)Given a sorted array (sorted in non-decreasing order) of positive numbers,
find the smallest positive integer value that cannot be represented as sum of elements of any subset of given set.

1) We decide that ‘res’ is the final result: 
If arr[i] is greater than ‘res’,
then we found the gap which is ‘res’ because the elements after arr[i] are also going to be greater than ‘res’.
2) The value of ‘res’ is incremented after considering arr[i]: 
The value of ‘res’ is incremented by arr[i] (why? If elements from 0 to (i-1) can represent 1 to ‘res-1’, 
then elements from 0 to i can represent from 1 to ‘res + arr[i] – 1’ be adding ‘arr[i]’ to all subsets that represent 1 to ‘res’)
*/
long long findSmallest(long long arr[], int n)
{
    long long res = 1;
    for (int i = 0; i < n; i++)
    {
        if (arr[i] <= res)
            res += arr[i];
        else
            break;
    }
    return res;
}

// Generate all possible sorted arrays from alternate elements of two given sorted arrays
void possibleArrays(vector<int> &A, vector<int> &B, int i, int j, int flag, vector<int> &smallAns, vector<vector<int>> &ans)
{
    //Add element from A
    if (flag == 0)
    {
        if (i >= A.size())
            return;
        //find index of next element in sorted order
        while (i < A.size() && A[i] < smallAns[smallAns.size() - 1])
            i++;
        //call recursion for every element after that
        for (int k = i; k < A.size(); k++)
        {
            smallAns.push_back(A[k]);
            possibleArrays(A, B, k + 1, j, 1, smallAns, ans);
            smallAns.pop_back();
        }
    }

    //Add element from B
    else
    {
        //same as A, but for B it is added to the ans array
        if (j >= B.size())
            return;
        while (j < B.size() && B[j] < smallAns[smallAns.size() - 1])
            j++;
        for (int k = j; k < B.size(); k++)
        {
            smallAns.push_back(B[k]);
            ans.push_back(smallAns);
            possibleArrays(A, B, i, k + 1, 0, smallAns, ans);
            smallAns.pop_back();
        }
    }
}
vector<vector<int>> possibleArrays(vector<int> &A, vector<int> &B)
{
    vector<int> smallAns;
    vector<vector<int>> ans;

    for (int i = 0; i < A.size(); i++)
    {
        smallAns.push_back(A[i]);
        possibleArrays(A, B, i + 1, 0, 1, smallAns, ans);
        smallAns.pop_back();
    }

    return ans;
}

// 31. Next Permutation
void nextPermutation(vector<int> &nums)
{
    if (nums.size() <= 1)
        return;
    int i = nums.size() - 1, minNum = 1e8, minIdx = i;
    //find the first number where arr[i]>arr[i-1]
    while (i > 0 && nums[i] <= nums[i - 1])
        i--;

    //if array is completely descending
    if (i == 0)
    {
        sort(nums.begin(), nums.end());
        return;
    }

    //find the element after i, that is just greater than arr[i]
    i--;
    for (int k = i; k < nums.size(); k++)
    {
        if (nums[k] > nums[i] && nums[k] < minNum)
        {
            minNum = nums[k];
            minIdx = k;
        }
    }

    //swap those two elements, and sort the sort from i+1 to n
    swap(nums[i], nums[minIdx]);
    sort(nums.begin() + i + 1, nums.end());
}

// 795. Number of Subarrays with Bounded Maximum
/*Approach- T: O(n), S: O(n)
Suppose dp[i] denotes the max number of valid sub-array ending with A[i]. We use following example to illustrate the idea:
A = [2, 1, 4, 2, 3], L = 2, R = 3

2,2,2,1,2,3,2,2,2

if A[i] < L:
For example, i = 1. We can only append A[i] to a valid sub-array ending with A[i-1] to create new sub-array.
So we have dp[i] = dp[i-1] (for i > 0)

if A[i] > R:
For example, i = 2. No valid sub-array ending with A[i] exist. So we have dp[i] = 0.
We also record the position of the invalid number 4 here as prev.

if L <= A[i] <= R
For example, i = 4. In this case any sub-array starts after the previous invalid number to A[i] (A[prev+1..i], A[prev+2..i]) is a new valid sub-array. 
So dp[i] += i - prev

Finally the sum of the dp array is the solution. Meanwhile, notice dp[i] only relies on dp[i-1] (and also prev), we can reduce the space complexity to O(1)
*/
int numSubarrayBoundedMax(vector<int> &A, int L, int R)
{
    int prev = 0;
    vector<int> dp(A.size() + 1, 0);
    for (int i = 1; i <= A.size(); i++)
    {
        if (A[i - 1] < L)
        {
            dp[i] = dp[i - 1];
        }
        else if (A[i - 1] >= L && A[i - 1] <= R)
        {
            dp[i] = i - prev;
        }
        else
        {
            prev = i;
            dp[i] = 0;
        }
    }

    int count = 0;
    for (int i = 0; i < dp.size(); i++)
        count += dp[i];

    return count;
}
//Approach- same as above, O(1) space
int numSubarrayBoundedMax(vector<int> &A, int L, int R)
{
    int prev = 0, prevCount = 0, count = 0;
    for (int i = 1; i <= A.size(); i++)
    {
        if (A[i - 1] < L)
            count += prevCount;
        else if (A[i - 1] >= L && A[i - 1] <= R)
        {
            prevCount = i - prev;
            count += prevCount;
        }
        else
        {
            prev = i;
            prevCount = 0;
        }
    }

    return count;
}

// Number of Sub-arrays of Size K
/*Approach-
simple recursion using subsequence method
We can avoid duplicate subarrays by using a set
OR
We can avoid duplicates by adding following two additional things to above code.
1) Add code to sort the array before calling combinationUtil() in printCombination()
2) Add following lines in the recursion
    // Since the elements are sorted, all occurrences of an element
    // must be together
        while (arr[i] == arr[i+1])
             i++; 
*/
int K;
unordered_set<vector<int>> res;
void numOfSubarrays(vector<int> &arr, int k, int idx, vector<int> &smallAns)
{
    if (k == 0)
    {
        res.insert(smallAns);
    }
    if (idx == arr.size())
        return;
    smallAns.push_back(arr[idx]);
    numOfSubarrays(arr, k - 1, idx + 1, smallAns);
    smallAns.pop_back();
    numOfSubarrays(arr, k, idx + 1, smallAns);
}
int numOfSubarrays(vector<int> &arr, int k, int threshold)
{
    K = k;
    vector<int> smallAns;
    numOfSubarrays(arr, k, 0, smallAns);
    return res.size();
}

// 689. Maximum Sum of 3 Non-Overlapping Subarrays
/*
Use Subarray sum array, prefix array, suffix array to get the three non overlapping windows
In prefix and suffix we store the index so that we can get lexicographically smallest
*/
vector<int> maxSumOfThreeSubarrays(vector<int> &nums, int k)
{
    vector<int> subArraySum(nums.size(), 0);
    int currSum = 0;

    //make the subarray for all window sums
    for (int i = 0; i < k; i++)
        currSum += nums[i];
    subArraySum[k - 1] = currSum;
    for (int i = k; i < nums.size(); i++)
    {
        currSum += nums[i] - nums[i - k];
        subArraySum[i] = currSum;
    }

    //make a prefix array for index of max window sum from 0 to till me
    vector<int> maxLeft(nums.size());
    int maxSum = subArraySum[0], maxIdx = 0;
    for (int i = 0; i < nums.size(); i++)
    {
        if (subArraySum[i] > maxSum)
        {
            maxSum = subArraySum[i];
            maxIdx = i;
        }
        maxLeft[i] = maxIdx;
    }

    //make suffix array for index of max window sum from me to n-1
    vector<int> maxRight(nums.size());
    maxSum = subArraySum[subArraySum.size() - 1];
    maxIdx = subArraySum.size() - 1;
    for (int i = nums.size() - 1; i >= 0; i--)
    {
        //we use >= as we need lexicograhically smallest, we need index of leftmost max sum on my right
        if (subArraySum[i] >= maxSum)
        {
            maxSum = subArraySum[i];
            maxIdx = i;
        }
        maxRight[i] = maxIdx;
    }

    //using the three arrays find the lexicographically smallest max sum
    vector<int> ans(3, -1);
    int maxTotalSum = -1e8;
    for (int i = k; i < nums.size() - k; i++)
    {
        //get the index of window sum of left, middle, right
        int l = maxLeft[i - k], m = i, r = maxRight[i + k];
        int currTotalSum = subArraySum[l] + subArraySum[m] + subArraySum[r];

        if (currTotalSum > maxTotalSum)
        {
            /*use (-k+1) to get the starting index of the window,
            as we are storing the last index of each window in all 3 arrays*/
            maxTotalSum = currTotalSum;
            ans[0] = l - k + 1;
            ans[1] = m - k + 1;
            ans[2] = r - k + 1;
        }
    }

    return ans;
}

// 974. Subarray Sums Divisible by K
/*
count the freq of all remainders you can get
then for each remainder, total subarrays that can be formed= count*(count-1)/2
*/
int subarraysDivByK(vector<int> &A, int K)
{
    vector<int> count(K); //array of count of remainders
    int res = 0;
    int sum = 0;
    //get all the counts
    for (int i = 0; i < A.size(); i++)
    {
        sum += A[i];
        //if the current sum%k==0 then increase ans by 1
        if (sum % K == 0)
            res++;
        //[(sum % K + K) % K] gives the proper index, we use (sum % K + K) % K because the sum can be negative,
        // so sum % K will give negative index
        count[(sum % K + K) % K]++;
    }

    for (int i = 0; i < count.size(); i++)
        res += count[i] * (count[i] - 1) / 2;

    return res;
}

// Find minimum number of merge operations to make an array palindrome
int minMerge(vector<int> &arr)
{
    int n = arr.size();
    int count = 0;
    int i = 0, j = n - 1;
    while (i < j)
    {
        // If right element is greater, then
        // we merge left two elements
        if (arr[i] < arr[j])
        {
            arr[i + 1] = arr[i + 1] + arr[i];
            count++;
            i++;
        }
        // If left element is greater, then
        // we merge right two elements
        else if (arr[i] > arr[j])
        {
            arr[j - 1] = arr[j - 1] + arr[j];
            j--;
            count++;
        }
        else
        {
            i++;
            j--;
        }
    }

    return count;
}

// Reorder an array according to given indexes
/*
Approach 1- O(n) space
just take an auxilliary array, and put everything at its right place in it

Approach 2- Sorting - Time: O(nlogn), Space: O(1)
When all are in their correct positions then the index array will be sorted.
So, we can just sort the index array, and sort the corresponding elements of the element array along with it
So for eg-
Input:  arr[]   = [50, 40, 70, 60, 90]
        index[] = [3,  0,  4,  1,  2]
Output: arr[]   = [40, 60, 90, 50, 70]
        index[] = [0,  1,  2,  3,   4] 

When 0 in index array comes to 0th index, with it 40 will too.
So after sorting they all will be in correct positions

Approach 3- Time: O(n) , Space: O(1)
for every index{
    while(the value at this index is not supposed to be at this index i.e index[i]!=i)
        we swap the value at this index with the value at its correct position

for eg-
at index 0-> 50,3, so we swap with values at index 3 i.e 60,1.
Now 60,1 should not be at 0th index, so swap it with 1th index values i.e. 40,0
Now 40,0 are at correct index and so we move on. 
}

*/
void reorder(vector<int> &arr, vector<int> &index)
{
    //Approach 3-
    for (int i = 0; i < arr.size(); i++)
    {
        while (index[i] != i)
        {
            swap(arr[i], arr[index[i]]);
            swap(index[i], index[index[i]]);
        }
    }
}

// Rearrange an array in maximum minimum form
/*
All even index have max element, and odd have min element
We encode each element with arr[i] += (arr[max_index] % max_element * max_element) for even index
and encode each element with arr[i] += (arr[min_index] % max_element * max_element) for odd index
How does expression “arr[i] += arr[max_index] % max_element * max_element” work ?
The purpose of this expression is to store two elements at index arr[i].
arr[max_index] is stored as multiplier and “arr[i]” is stored as remainder.
For example in {1 2 3 4 5 6 7 8 9}, max_element is 10 and we store 91 at index 0.
With 91, we can get original element as 91%10 and new element as 91/10.
*/
void rearrange(vector<long> &arr)
{
    int n = arr.size();
    long maxIdx = n - 1, minIdx = 0, maxEle = arr[n - 1] + 1;
    for (int i = 0; i < n; i++)
    {
        //encode for even index with maxIdx
        //arr[maxIdx] % maxEle gives the original element at maxIdx,
        //then we encode it with arr[i] += originalElement * maxEle
        if (i % 2 == 0)
        {
            arr[i] += (arr[maxIdx] % maxEle) * maxEle;
            maxIdx--;
        }
        //encode for even index with minIdx
        else
        {
            arr[i] += (arr[minIdx] % maxEle) * maxEle;
            minIdx++;
        }
    }

    //Recover the array
    for (int i = 0; i < n; i++)
        arr[i] /= maxEle;

    for (int i = 0; i < n; i++)
        cout << arr[i] << " ";
    cout << endl;
}

// 442. Find All Duplicates in an Array
/*
For every number make the element at arr[arr[i]] -ve
if it is already -ve, then it is repeating
*/
vector<int> findDuplicates(vector<int> &nums)
{
    int n = nums.size();
    vector<int> res;

    //we can mark nums[nums[i]-1] as -ve, so we wont have to check for n seperately
    // int flag = 0; //to check if the max element i.e. n has been seen already
    for (int i = 0; i < n; i++)
    {
        // if (abs(nums[i]) == n)
        // {
        //     if (flag == 1)
        //         res.push_back(n);
        //     else
        //         flag = 1;
        // }
        if (nums[abs(nums[i]) - 1] < 0)
            res.push_back(abs(nums[i]));
        else
            nums[abs(nums[i]) - 1] *= -1;
    }

    return res;
}

// Find duplicates in O(n) time and O(1) extra space (same element can ve present more than 2 times)
/*
There is a problem in the above approach.
It prints the repeated number more than once. For example: {1, 6, 3, 1, 3, 6, 6} it will give output as : 1 3 6 6.
In below set, another approach is discussed that prints repeating elements only once.
Approach: The basic idea is to use a HashMap to solve the problem. But there is a catch, the numbers in the array are from 0 to n-1,
and the input array has length n. So, the input array can be used as a HashMap. While Traversing the array,
if an element ‘a’ is encountered then increase the value of a%n‘th element by n. The frequency can be retrieved by dividing the a % n’th element by n.
Algorithm:
Traverse the given array from start to end.
For every element in the array increment the arr[i]%n‘th element by n.
Now traverse the array again and print all those indexes i for which arr[i]/n is greater than 1. Which guarantees that the number n has been added to that index
This approach works because all elements are in the range from 0 to n-1 and arr[i] would be greater than n only if a value “i” has appeared more than once.
*/
vector<int> duplicates(int arr[], int n)
{
    vector<int> res;
    for (int i = 0; i < n; i++)
        arr[arr[i] % n] += n;

    for (int i = 0; i < n; i++)
    {
        //if it occurs once, then it will only be arr[i]+=n once, so if it is >n*2 then there are duplicates
        if (arr[i] >= n * 2)
            res.push_back(i);
    }

    if (res.size() == 0)
        return {-1};
    return res;
}

// 239. Sliding Window Maximum
/*Approach 1- Using PQ O(nlogn)
if we remove the out of bound element each time, and maintain a PQ of size K
then time is O(nlogk)
but for that we have to implement our own PQ
*/
vector<int> maxSlidingWindow(vector<int> &nums, int k)
{
    vector<int> res;
    priority_queue<vector<int>, vector<vector<int>>> pq; //max heap

    int i = 0;
    //add elements of first window to PQ
    for (i = 0; i < k; i++)
        pq.push({nums[i], i});
    //top of PQ will contain the window max
    res.push_back(pq.top()[0]);

    for (; i < nums.size(); i++)
    {
        //push current element to PQ
        pq.push({nums[i], i});

        //while the top of PQ i.e. max element is out of current window, remove it
        while (pq.top()[1] <= i - k)
            pq.pop();

        //push the current max element into res
        res.push_back(pq.top()[0]);
    }

    return res;
}
//Approach 2- Using Deque O(n), or any double ended queue like data structure
/*
Add the index in queue
at every index, if the front of queue has element out of current window i.e. <=(i-k), then remove it
before adding the new element at the end of queue, remove all elements less than current element from the back
so the front of queue will always have the max element of the current window
*/
vector<int> maxSlidingWindow(vector<int> &nums, int k)
{
    vector<int> res;
    list<int> que;
    int i;

    //add elements of first window
    for (i = 0; i < k; i++)
    {
        //remove all elements from end that are less than current element
        while (que.size() > 0 && nums[que.back()] < nums[i])
            que.pop_back();
        que.push_back(i);
    }

    res.push_back(nums[que.front()]);

    for (; i < nums.size(); i++)
    {
        //remove first element if its index is out of current window
        if (que.size() > 0 && que.front() <= (i - k))
            que.pop_front();

        //remove all elements from end that are less than current element
        while (que.size() > 0 && nums[que.back()] < nums[i])
            que.pop_back();

        //add the current element
        que.push_back(i);

        //front of que is max element of window
        res.push_back(nums[que.front()]);
    }

    return res;
}

// 1375. Bulb Switcher III
int numTimesAllBlue(vector<int> &light)
{
    int n = light.size();
    int blueBulbs = 0, onBulbs = 0, count = 0;
    vector<int> bulbs(n, 0); //1: ON, 2: blue, 0: OFF

    for (int i = 0; i < light.size(); i++)
    {
        int bulbNo = light[i] - 1;
        if (bulbNo == 0)
        {
            bulbs[0] = 2;
            blueBulbs++;
            int b = bulbNo + 1;
            while (b < n && bulbs[b] == 1)
            {
                bulbs[b] = 2;
                b++;
                blueBulbs++;
            }
        }
        else if (bulbs[bulbNo - 1] == 2)
        {
            bulbs[bulbNo] = 2;
            blueBulbs++;
            int b = bulbNo + 1;
            while (b < n && bulbs[b] == 1)
            {
                bulbs[b] = 2;
                b++;
                blueBulbs++;
            }
        }
        else
            bulbs[bulbNo] = 1;

        onBulbs++;

        if (blueBulbs == onBulbs)
            count++;
    }

    return count;
}
//Better Approach - O(n), O(1)
/*
Only moment when all the bulbs are blue is if,
The rightmost turned on bulb is blue, as this means that every bulb before it is blue.
the condition for this is (rightMostBulb == i+1)
Because at any index i, if the rightmost ON bulb is more than i+1(because of 1 based index)
then all bulbs are not blue, as rightMost bulb>i+1,
which means total ON bulbs is less than rightMost ON bulb.
*/
int numTimesAllBlue(vector<int> &light)
{
    int rightMostBulb = light[0], count = 0;

    for (int i = 0; i < light.size(); i++)
    {
        rightMostBulb = max(light[i], rightMostBulb);
        if (rightMostBulb == i + 1)
            count++;
    }

    return count;
}

// 769. Max Chunks To Make Sorted
/*
The smallest chunk is where for range (a....b), all numbers in range [a,b] are within that chunk
So if the maxTillNow==i, this means everything in range (0...maxTillNow) has occured
So this gives the smallest chunk.
For eg - [0,4,3,1,2]
Here there are two chunks : (0) , (4,3,1,2)

Because for increasing order, at any index i, if ith element is not maxtillNow, that means everything till now is not sorted
So they are all part of one chunk.
*/
int maxChunksToSorted(vector<int> &arr)
{
    int n = arr.size();
    int maxTillNow = 0, count = 0;

    for (int i = 0; i < arr.size(); i++)
    {
        maxTillNow = max(arr[i], maxTillNow);
        if (maxTillNow == i)
            count++;
    }

    return count;
}

void solve()
{
}
int main()
{
    // g++ Arrays.cpp -o out && ./out < input.txt > output.txt
    solve();
    return 0;
}