#include <iostream>
#include <vector>
#include <list>
#include <algorithm>

using namespace std;

// 79. Word Search
int dir[4][4] = {{0, 1}, {1, 0}, {-1, 0}, {0, -1}};
bool exist(vector<vector<char>> &board, int sr, int sc, string &target, int idx)
{
    if (idx == target.size())
        return true;

    bool res = false;

    for (int d = 0; d < 4; d++)
    {
        int x = sr + dir[d][0];
        int y = sc + dir[d][1];

        if (x >= 0 && y >= 0 && x < board.size() && y < board[0].size() && board[x][y] != '#' && target[idx] == board[x][y])
        {
            char val = board[x][y];
            board[x][y] = '#'; //mark visited
            res = res || exist(board, x, y, target, idx + 1);
            board[x][y] = val; //unmark
        }
    }

    return res;
}
bool exist(vector<vector<char>> &board, string word)
{
    bool res = false;
    int n = board.size();
    int m = board[0].size();
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            char val = board[i][j];
            if (board[i][j] == word[0])
            {
                board[i][j] = '#';
                if (exist(board, i, j, word, 1))
                    return true;
                board[i][j] = val;
            }
        }
    }

    return false;
}

// 212. Word Search II (Passed all cases, but took too long)
int dir[4][4] = {{0, 1}, {1, 0}, {-1, 0}, {0, -1}};
bool exist_(vector<vector<char>> &board, int sr, int sc, string &target, int idx)
{
    if (idx == target.size())
        return true;

    bool res = false;

    for (int d = 0; d < 4; d++)
    {
        int x = sr + dir[d][0];
        int y = sc + dir[d][1];

        if (x >= 0 && y >= 0 && x < board.size() && y < board[0].size() && board[x][y] != '#' && target[idx] == board[x][y])
        {
            char val = board[x][y];
            board[x][y] = '#'; //mark visited
            res = res || exist_(board, x, y, target, idx + 1);
            board[x][y] = val; //unmark
        }
    }

    return res;
}
bool exist_(vector<vector<char>> &board, string &word)
{
    bool res = false;
    int n = board.size();
    int m = board[0].size();
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            char val = board[i][j];
            if (board[i][j] == word[0])
            {
                board[i][j] = '#';
                if (exist_(board, i, j, word, 1))
                {
                    board[i][j] = val;
                    return true;
                }
                board[i][j] = val;
            }
        }
    }

    return false;
}
vector<string> findWords(vector<vector<char>> &board, vector<string> &words)
{
    vector<string> res;
    for (int i = 0; i < words.size(); i++)
    {
        if (exist_(board, words[i]))
            res.push_back(words[i]);
    }
    return res;
}

void solve()
{
}
int main()
{
    solve();
    return 0;
}