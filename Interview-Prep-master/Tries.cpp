#include <iostream>
#include <vector>
#include <unordered_map>
#include <queue>

using namespace std;

//Trie Implementation===============================================================
class Node
{
public:
    int wordsEnd;
    vector<Node *> childs;

    Node()
    {
        this->wordsEnd = 0;
        this->childs.assign(26, nullptr);
    }
};

Node *root = new Node();

void insert(string word)
{
    Node *curr = root;
    for (int i = 0; i < word.length(); i++)
    {
        int ch = word[i] - 'a';

        if (curr->childs[ch] == nullptr)
            curr->childs[ch] = new Node();
        curr = curr->childs[ch];
    }

    curr->wordsEnd++;
}

bool search(string word)
{
    Node *curr = root;

    for (int i = 0; i < word.length(); i++)
    {
        int ch = word[i] - 'a';

        if (curr->childs[ch] == nullptr)
            return false;

        curr = curr->childs[ch];
    }

    return curr->wordsEnd > 0;
}

//To do
void deleteWord(string word);

//LeetCode Questions==========================================================================

// 208. Implement Trie (Prefix Tree)
class Trie
{
public:
    class Node
    {
    public:
        int wordsEnd;
        vector<Node *> childs;

        Node()
        {
            this->wordsEnd = 0;
            this->childs.assign(26, nullptr);
        }
    };

    Node *root;
    /** Initialize your data structure here. */
    Trie()
    {
        root = new Node();
    }

    /** Inserts a word into the trie. */
    void insert(string word)
    {
        Node *curr = root;
        for (int i = 0; i < word.length(); i++)
        {
            int idx = word[i] - 'a';

            if (curr->childs[idx] == nullptr)
                curr->childs[idx] = new Node();

            curr = curr->childs[idx];
        }

        curr->wordsEnd++;
    }

    /** Returns if the word is in the trie. */
    bool search(string word)
    {
        Node *curr = root;

        for (int i = 0; i < word.length(); i++)
        {
            int idx = word[i] - 'a';

            if (curr->childs[idx] == nullptr)
                return false;

            curr = curr->childs[idx];
        }

        return curr->wordsEnd > 0;
    }

    /** Returns if there is any word in the trie that starts with the given prefix. */
    bool startsWith(string prefix)
    {
        Node *curr = root;

        for (int i = 0; i < prefix.length(); i++)
        {
            int idx = prefix[i] - 'a';

            if (curr->childs[idx] == nullptr)
                return false;

            curr = curr->childs[idx];
        }

        return true;
    }
};

// 211. Design Add and Search Words Data Structure
class WordDictionary
{
public:
    class Node
    {
    public:
        int wordsEnd;
        vector<Node *> childs;

        Node()
        {
            this->wordsEnd = 0;
            this->childs.assign(26, nullptr);
        }
    };

    Node *root;
    /** Initialize your data structure here. */
    WordDictionary()
    {
        root = new Node();
    }

    /** Adds a word into the data structure. */
    void addWord(string word)
    {
        Node *curr = root;
        for (int i = 0; i < word.length(); i++)
        {
            int idx = word[i] - 'a';

            if (curr->childs[idx] == nullptr)
                curr->childs[idx] = new Node();

            curr = curr->childs[idx];
        }

        curr->wordsEnd++;
    }

    /** Returns if the word is in the data structure. A word could contain the dot character '.' to represent any one letter. */
    bool search(string &word, Node *node, int idx)
    {
        if (idx == word.length())
            return node->wordsEnd != 0;

        bool res = false;

        if (word[idx] == '.')
        {
            for (int i = 0; i < 26; i++)
            {
                if (node->childs[i] != nullptr)
                    res = res || search(word, node->childs[i], idx + 1);
            }
        }
        else
        {
            if (node->childs[word[idx] - 'a'] == nullptr)
                return false;

            res = res || search(word, node->childs[word[idx] - 'a'], idx + 1);
        }

        return res;
    }

    bool search(string word)
    {
        return search(word, root, 0);
    }
};

// 212. Word Search II
class Node
{
public:
    bool wordsEnd;
    string word;
    vector<Node *> childs;

    Node()
    {
        this->wordsEnd = false;
        this->childs.assign(26, nullptr);
        this->word = "";
    }
};
Node *root = nullptr;

vector<string> res;
int dir[4][2] = {{0, 1}, {0, -1}, {1, 0}, {-1, 0}};

void findWords_(vector<vector<char>> &board, int sr, int sc, Node *node)
{
    if (node->wordsEnd == true)
    {
        res.push_back(node->word);
        node->wordsEnd = false;
    }

    int n = board.size();
    int m = board[0].size();
    char ch = board[sr][sc];
    board[sr][sc] = '#';
    for (int d = 0; d < 4; d++)
    {
        int x = sr + dir[d][0];
        int y = sc + dir[d][1];

        if (x >= 0 && y >= 0 && x < n && y < m && board[x][y] != '#' && node->childs[board[x][y] - 'a'] != nullptr)
        {
            findWords_(board, x, y, node->childs[board[x][y] - 'a']);
        }
    }
    board[sr][sc] = ch;
}
void insert(string word)
{
    Node *curr = root;
    for (int i = 0; i < word.length(); i++)
    {
        int ch = word[i] - 'a';

        if (curr->childs[ch] == nullptr)
            curr->childs[ch] = new Node();
        curr = curr->childs[ch];
    }

    curr->wordsEnd = true;
    curr->word = word;
}
vector<string> findWords(vector<vector<char>> &board, vector<string> &words)
{
    if (board.size() == 0 || board[0].size() == 0 || words.size() == 0)
        return {};

    root = new Node();

    int n = board.size();
    int m = board[0].size();

    for (int i = 0; i < words.size(); i++)
    {
        insert(words[i]);
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            if (root->childs[board[i][j] - 'a'] != nullptr)
            {
                findWords_(board, i, j, root->childs[board[i][j] - 'a']);
            }
        }
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