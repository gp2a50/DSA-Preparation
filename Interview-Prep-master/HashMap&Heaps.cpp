#include <iostream>
#include <vector>
#include <list>
#include <unordered_map>
#include <algorithm>
#include <queue>

using namespace std;

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

    for (int i = 0; i < nums1.size(); i++)
    {
        pq.push({nums1[i], nums2[0], 0});
    }

    while (k-- > 0 && pq.size() != 0)
    {
        vector<int> pr = pq.top();
        pq.pop();

        ans.push_back({pr[0], pr[1]});

        if (pr[2] == nums2.size() - 1)
            continue;

        pq.push({pr[0], nums2[pr[2] + 1], pr[2] + 1});
    }

    return ans;
}

// 692. Top K Frequent Words
struct myCompare
{
    bool operator()(const pair<int, string> &pr1, const pair<int, string> &pr2)
    {
        //if freq of words is equal, sort lexicographically
        if (pr1.first == pr2.first)
            return pr1.second > pr2.second;

        //else sort according to freq
        return pr1.first < pr2.first;
    }
};
vector<string> topKFrequent(vector<string> &words, int k)
{
    vector<string> res;
    unordered_map<string, int> freq;
    priority_queue<pair<int, string>, vector<pair<int, string>>, myCompare> pq;

    //make freq map for words
    for (int i = 0; i < words.size(); i++)
        freq[words[i]]++;

    //put all (freq,words) into PQ using own comparator
    for (auto ele : freq)
    {
        pq.push({ele.second, ele.first});
    }

    //We can only store the top K elements in PQ and save space,
    //then we will have to reverse the result at the end, as it is in increasing order and we want decreasing
    //Reverse the signs in comparator too
    // bool operator()(const pair<int, string> &pr1, const pair<int, string> &pr2)
    // {
    //     if (pr1.first == pr2.first)
    //     {
    //         return pr1.second < pr2.second;
    //     }

    //     return pr1.first > pr2.first;
    // }
    // for (auto ele : freq)
    // {
    //     pq.push({ele.second, ele.first});
    //     if(pq.size()>k){
    //         pq.pop();
    //     }
    // }

    //add the top K elements of PQ into result
    while (k-- > 0)
    {
        string word = pq.top().second;
        pq.pop();
        res.push_back(word);
    }

    // reverse(res.begin(),res.end());

    return res;
}


void solve()
{
}
int main()
{
    // g++ Heaps.cpp -o out && ./out < input.txt > output.txt
    solve();
    return 0;
}