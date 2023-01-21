#include <iostream>
#include <vector>
#include <list>
#include <algorithm>

using namespace std;
// Friends Pairing Problem
long MOD = 1000000007;
int countFriendsPairings_Mem(int n, vector<int> &dp)
{
    if (n <= 1)
        return dp[n] = 1;

    if (dp[n] != 0)
        return dp[n];

    int single = (countFriendsPairings_Mem(n - 1, dp)) % MOD;
    int pairUp = (((countFriendsPairings_Mem(n - 2, dp)) % MOD) * ((n - 1) % MOD)) % MOD;

    return dp[n] = ((single % MOD) + (pairUp % MOD)) % MOD;
}
int countFriendsPairings_DP(int N, vector<int> &dp)
{
    //space can be O(1), if we keep track of just last two values i.e. n-1,n-2
    for (int n = 0; n <= N; n++)
    {
        if (n <= 1)
        {
            dp[n] = 1;
            continue;
        }

        int single = dp[n - 1] % MOD;
        int pairUp = ((dp[n - 2] % MOD) * ((n - 1) % MOD)) % MOD;

        dp[n] = ((single % MOD) + (pairUp % MOD)) % MOD;
    }

    return dp[N];
}
int countFriendsPairings_01(int N)
{
    //space can be O(1), if we keep track of just last two values i.e. n-1,n-2
    if (N <= 1)
        return 1;

    int a = 1, b = 1, c; //a=n-2, b=n-1
    for (int n = 2; n <= N; n++)
    {
        // c = b + a * (n - 1); //dp[n - 2] * (n - 1)
        c = (b % MOD + (a % MOD * (n - 1) % MOD) % MOD) % MOD;
        a = b % MOD; //dp[n - 2]
        b = c % MOD; //dp[n - 1]
    }

    return c;
}
int countFriendsPairings(int n)
{
    vector<int> dp(n + 1, 0);
    // return countFriendsPairings_Rec(n);
    // return countFriendsPairings_Mem(n, dp);
    return countFriendsPairings_DP(n, dp);
}

// 1219. Path with Maximum Gold
vector<vector<int>> dirA = {{0, 1}, {1, 0}, {-1, 0}, {0, -1}};
int multiPath(int sr, int sc, int er, int ec, vector<vector<int>> &board)
{
    if (sc == ec && sr == er)
    {
        return board[er][ec];
    }
    board[sr][sc] *= -1;
    int maxgold = 0;
    for (int i = 0; i < dirA.size(); i++)
    {

        int x = sr + dirA[i][0];
        int y = sc + dirA[i][1];
        if (x >= 0 && y >= 0 && x < board.size() && y < board[0].size() && board[x][y] > 0)
        {
            int gold = multiPath(x, y, er, ec, board);
            maxgold = max(gold, maxgold);
        }
    }
    board[sr][sc] *= -1;
    return maxgold + board[sr][sc];
}
int getMaximumGold(vector<vector<int>> &grid)
{
    int currmax = 0;
    for (int i = 0; i < grid.size(); i++)
    {
        for (int j = 0; j < grid[0].size(); j++)
        {
            if (grid[i][j] != 0)
            {
                int currGold = multiPath(i, j, grid.size(), grid[0].size(), grid);
                currmax = max(currGold, currmax);
            }
        }
    }
    return currmax;
}

//Strings=========================================================================================
// 44. Wildcard Matching
/*
Approach 1- DP
if s[i]==p[j] => n-1, m-1
if s[i]!=p[j]{
    if '?' => n-1, m-1
    if '* => n-1,m and n-1,m-1 and n,m-1(for "" string)
    also in '* wildcard, can be subtituted for empty "" string in begin and end
}

Approach 2- Greedy Solution(not done)

*/
bool isMatch_01(string s, string p, int n, int m)
{
    //both strings are over
    if (n < 0 && m < 0)
        return true;
    if (n < 0)
    {
        //for case: s="", p="***", i.e pattern just has wildcards left
        //then all wildcards are substituted for empty string
        while (m >= 0 && p[m] == '*')
            m--;

        if (n < 0 && m < 0)
            return true;
        else
            return false;
    }

    bool res = false;
    if (s[n] == p[m])
    {
        res = res || isMatch_01(s, p, n - 1, m - 1);
    }
    else
    {
        if (p[m] == '?')
        {
            res = res || isMatch_01(s, p, n - 1, m - 1);
        }
        else if (p[m] == '*')
        {
            res = res || isMatch_01(s, p, n - 1, m) || isMatch_01(s, p, n - 1, m - 1) || isMatch_01(s, p, n, m - 1);
        }
    }

    return res;
}
int isMatch_Mem(string &s, string &p, int n, int m, vector<vector<int>> &dp)
{
    //both strings are over
    if (n < 0 && m < 0)
        return 1;
    if (m < 0)
        return 0;
    if (n < 0)
    {
        //for case: s="", p="***", i.e pattern just has wildcards left
        //then all wildcards are substituted for empty string
        while (m >= 0 && p[m] == '*')
            m--;

        if (n < 0 && m < 0)
            return 1;
        else
            return 0;
    }

    if (dp[n][m] != -1)
        return dp[n][m];

    bool res = false;
    if (s[n] == p[m])
        res = res || isMatch_Mem(s, p, n - 1, m - 1, dp);
    else
    {
        if (p[m] == '?')
            res = res || isMatch_Mem(s, p, n - 1, m - 1, dp);
        else if (p[m] == '*')
            res = res || isMatch_Mem(s, p, n - 1, m, dp) || isMatch_Mem(s, p, n - 1, m - 1, dp) || isMatch_Mem(s, p, n, m - 1, dp);
    }

    return dp[n][m] = (res ? 1 : 0);
}

bool isMatch(string s, string p)
{
    int n = s.size();
    int m = p.size();

    vector<vector<int>> dp(n, vector<int>(m, -1));
    return isMatch_Mem(s, p, n-1, m-1, dp);
}

void solve()
{
}
int main()
{
    solve();
    return 0;
}