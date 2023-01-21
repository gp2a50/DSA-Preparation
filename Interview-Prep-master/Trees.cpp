#include <iostream>
#include <vector>
#include <algorithm>
#include <list>
#include <math.h>
#include <queue>
#include <stack>
#include <unordered_map>

using namespace std;

class TreeNode
{
public:
    int val;
    TreeNode *left, *right;
    TreeNode(int val)
    {
        this->val = val;
        this->left = this->right = NULL;
    }
};

struct ListNode
{
    int val;
    ListNode *next;
    ListNode() : val(0), next(nullptr) {}
    ListNode(int x) : val(x), next(nullptr) {}
    ListNode(int x, ListNode *next) : val(x), next(next) {}
};

struct Node
{
    int data;
    Node *left;
    Node *right;

    Node(int val)
    {
        data = val;
        left = right = NULL;
    }
};

class Node
{
public:
    int val;
    int data; //no use, just to remove errors due to two Node classes
    Node *left;
    Node *right;
    Node *next;

    Node() : val(0), left(NULL), right(NULL), next(NULL) {}

    Node(int _val) : val(_val), left(NULL), right(NULL), next(NULL) {}

    Node(int _val, Node *_left, Node *_right, Node *_next)
        : val(_val), left(_left), right(_right), next(_next) {}
};
//Traversal type=================================================================================

//Views==========================================================================================
//878. Boundary of Binary Tree (LintCode)
int leftMin = 0;
int rightMax = 0;
void width(TreeNode *root, int level)
{
    if (root == nullptr)
        return;

    leftMin = min(leftMin, level);
    rightMax = max(rightMax, level);

    width(root->left, level - 1);
    width(root->right, level + 1);
}
void getLeftNodes(TreeNode *root, int level, vector<int> &leftNodes)
{
    if (root == nullptr)
        return;

    leftNodes.push_back(root->val);
    if (level == leftMin && root->left == nullptr && root->right == nullptr)
        return;
    if (root->left == nullptr)
        getLeftNodes(root->right, level + 1, leftNodes);
    else
        getLeftNodes(root->left, level - 1, leftNodes);
}
void getrightNodes(TreeNode *root, int level, vector<int> &rightNodes)
{
    if (root == nullptr)
        return;

    rightNodes.push_back(root->val);
    if (level == rightMax && root->left == nullptr && root->right == nullptr)
        return;
    if (root->right == nullptr)
        getrightNodes(root->left, level - 1, rightNodes);
    else
        getrightNodes(root->right, level + 1, rightNodes);
}
void getLeafNodes(TreeNode *root, vector<int> &leafNodes)
{
    if (root == nullptr)
        return;

    if (root->left == nullptr && root->right == nullptr)
    {
        leafNodes.push_back(root->val);
        return;
    }

    getLeafNodes(root->left, leafNodes);
    getLeafNodes(root->right, leafNodes);
}
vector<int> boundaryOfBinaryTree(TreeNode *root)
{
    leftMin = 0;
    rightMax = 0;
    width(root, 0);

    vector<int> leafNodes;
    vector<int> leftNodes;
    vector<int> rightNodes;
    vector<int> boundary;

    getLeafNodes(root, leafNodes);
    if (root->left == nullptr)
        leftNodes.push_back(root->val);
    else
        getLeftNodes(root, 0, leftNodes);
    if (root->right == nullptr)
        rightNodes.push_back(root->val);
    else
        getrightNodes(root, 0, rightNodes);

    for (int ele : leftNodes)
        boundary.push_back(ele);

    if (leftNodes[leftNodes.size() - 1] == leafNodes[0])
        boundary.pop_back();
    for (int i = 0; i < leafNodes.size(); i++)
        boundary.push_back(leafNodes[i]);

    if (rightNodes[rightNodes.size() - 1] == leafNodes[leafNodes.size() - 1])
        boundary.pop_back();
    for (int i = rightNodes.size() - 1; i > 0; i--)
        boundary.push_back(rightNodes[i]);

    return boundary;
}

//Construct Type=================================================================================
// 606. Construct String from Binary Tree
string tree2str(TreeNode *t)
{
    if (t == nullptr)
        return "";
    if (t->left == nullptr && t->right == nullptr)
        return to_string(t->val);

    string ans = to_string(t->val);
    ans += '(';
    ans += tree2str(t->left);
    ans += ')';
    if (t->right != nullptr)
    {
        ans += '(';
        ans += tree2str(t->right);
        ans += ')';
    }

    return ans;
}

// 880. Construct Binary Tree from String (LintCode) (TLE)
TreeNode *str2tree(string &s)
{
    if (s.size() == 0)
        return nullptr;

    stack<TreeNode *> st;
    int i = 0;
    if (s[0] == '-')
    {
        i++;
        st.push(new TreeNode((s[i++] - '0') * -1));
    }
    else
        st.push(new TreeNode(s[i++] - '0'));

    while (i < s.size())
    {
        if (s[i] == '(')
        {
            i++;
            int nodeVal = 0;
            bool isPos = true;
            if (s[i] == '-')
            {
                i++;
                isPos = false;
            }

            while (s[i] != '(' && s[i] != ')')
            {
                nodeVal = nodeVal * 10 + (s[i] - '0');
                i++;
            }
            if (!isPos)
                nodeVal *= -1;
            TreeNode *node = new TreeNode(nodeVal);
            st.push(node);
        }
        else if (s[i] == ')')
        {
            TreeNode *child = st.top();
            st.pop();
            TreeNode *root = st.top();
            if (root->left == nullptr)
                root->left = child;
            else
                root->right = child;
            i++;
        }
    }

    return st.top();
}

//Misc.==========================================================================================
//Convert Binary Tree to Circular Linked List (LintCode)
TreeNode *head = nullptr;
TreeNode *tail = nullptr;
void DLL(TreeNode *root)
{
    if (root == nullptr)
        return;

    DLL(root->left);

    if (head == nullptr)
    {
        head = root;
    }
    else
    {
        tail->right = root;
        root->left = tail;
    }

    tail = root;

    DLL(root->right);
}
TreeNode *treeToDoublyList(TreeNode *root)
{
    //change bt to DLL
    DLL(root);

    //make it circular
    tail->right = head;
    head->left = tail;

    return head;
}

//109. Convert Sorted List to Binary Search Tree
ListNode *node = nullptr;
TreeNode *LLtoBST(int n)
{
    if (n == 0)
        return nullptr;

    TreeNode *leftChild = LLtoBST(n / 2);

    TreeNode *root = new TreeNode(node->val);
    root->left = leftChild;
    node = node->next;

    TreeNode *rightChild = LLtoBST(n - n / 2 - 1);

    root->right = rightChild;

    return root;
}
TreeNode *sortedListToBST(ListNode *head)
{
    node = head;
    int size = 0;

    ListNode *curr = head;
    while (curr != nullptr)
    {
        size++;
        curr = curr->next;
    }

    TreeNode *root = LLtoBST(size);
    return root;
}

//Convert DLL to BST (not tested, not submit found), use left and right as next and prev
TreeNode *temp = nullptr;
TreeNode *DLLtoBST(int n)
{
    if (n == 0)
        return nullptr;

    TreeNode *leftChild = DLLtoBST(n / 2);

    TreeNode *root = temp;
    root->left = leftChild;
    temp = temp->right;

    TreeNode *rightChild = DLLtoBST(n - n / 2 - 1);

    root->right = rightChild;

    return root;
}
TreeNode *sortedDLLToBST(TreeNode *head)
{
    temp = head;
    int size = 0;

    TreeNode *curr = head;
    while (curr != nullptr)
    {
        size++;
        curr = curr->right;
    }

    TreeNode *root = DLLtoBST(size);
    return root;
}

//Merge two BSTs(not complete)
TreeNode *head = nullptr;
TreeNode *tail = nullptr;
void DLL(TreeNode *root)
{
    if (root == nullptr)
        return;

    DLL(root->left);

    if (head == nullptr)
    {
        head = root;
    }
    else
    {
        tail->right = root;
        root->left = tail;
    }

    tail = root;

    DLL(root->right);
}
TreeNode *treeToDoublyList(TreeNode *root)
{
    //change bt to DLL
    DLL(root);

    return head;
}
TreeNode *mergeBST(TreeNode *root1, TreeNode *root2)
{
    //convert to DLLs - gives head and tail of both
    TreeNode *head1 = treeToDoublyList(root1);
    TreeNode *head2 = treeToDoublyList(root2);

    //merge two DLLs in sorted order
    TreeNode *nhead = new TreeNode(-1);
    // while();

    //convert merged DLL to BST
    TreeNode *nroot = sortedDLLToBST(nhead);

    return nroot;
}

//1038. Binary Search Tree to Greater Sum Tree
int currSum;
void modify_(TreeNode *root)
{
    if (root == nullptr)
        return;

    modify_(root->right);

    root->val += currSum;
    currSum = root->val;

    modify_(root->left);
}
TreeNode *bstToGst(TreeNode *root)
{
    currSum = 0;
    modify_(root);
    return root;
}

// 623. Add One Row to Tree
TreeNode *addOneRow(TreeNode *root, int v, int d)
{
    if (d == 1)
    {
        TreeNode *nroot = new TreeNode(v);
        nroot->left = root;
        return nroot;
    }

    list<TreeNode *> que;
    que.push_back(root);

    while (que.size() != 0)
    {
        int size = que.size();
        d--;
        while (size-- > 0)
        {
            TreeNode *node = que.front();
            que.pop_front();

            if (d == 1)
            {
                TreeNode *leftNode = new TreeNode(v);
                TreeNode *rightNode = new TreeNode(v);

                leftNode->left = node->left;
                rightNode->right = node->right;

                node->left = leftNode;
                node->right = rightNode;
            }

            else
            {
                if (node->left != nullptr)
                    que.push_back(node->left);
                if (node->right != nullptr)
                    que.push_back(node->right);
            }
        }
        if (d == 0)
            break;
    }

    return root;
}

// Ancestors in Binary Tree
bool findAncestors(Node *root, int tar, vector<int> &ans)
{
    if (root == nullptr)
    {
        return false;
    }
    if (root->data == tar)
    {
        return true;
    }

    bool leftAns = findAncestors(root->left, tar, ans);
    if (leftAns)
    {
        ans.push_back(root->data);
        return true;
    }
    bool rightAns = findAncestors(root->right, tar, ans);
    if (rightAns)
    {
        ans.push_back(root->data);
        return true;
    }
}
vector<int> Ancestors(Node *root, int target)
{
    vector<int> ans;
    findAncestors(root, target, ans);
    return ans;
}

// 637. Average of Levels in Binary Tree
vector<double> averageOfLevels(TreeNode *root)
{
    list<TreeNode *> que;
    vector<double> ans;
    que.push_back(root);

    while (que.size() > 0)
    {
        int size = que.size();
        double s = size;
        double sum = 0;
        while (size-- > 0)
        {
            TreeNode *rnode = que.front();
            que.pop_front();

            sum += rnode->val;
            if (rnode->left != nullptr)
                que.push_back(rnode->left);
            if (rnode->right != nullptr)
                que.push_back(rnode->right);
        }
        ans.push_back(sum / s);
    }
    return ans;
}

// 110. Balanced Binary Tree
bool res;
int isBalanced_(TreeNode *root)
{
    if (root == nullptr)
        return -1;
    if (!res) //to stop unnecessary calls
        return 0;

    int lh = isBalanced_(root->left);
    int rh = isBalanced_(root->right);

    if (abs(lh - rh) > 1)
        res = false;

    return max(lh, rh) + 1;
}
bool isBalanced(TreeNode *root)
{
    res = true;
    isBalanced_(root);
    return res;
}

// BBT counter(G4G submit section answers are wrong, they are different from their editorial)
long MOD = 1000000007;
long long countBBT(int h, vector<long long> &dp)
{
    if (h == 0 || h == 1)
        return dp[h] = 1;

    if (dp[h] != -1)
        return dp[h];

    return dp[h] = ((countBBT(h - 1, dp) % MOD) * (2 * (countBBT(h - 2, dp)) % MOD + countBBT(h - 1, dp) % MOD) % MOD) % MOD;
}
void solve()
{
    int t;
    cin >> t;
    while (t-- > 0)
    {
        int h;
        cin >> h;
        vector<long long> dp(h, -1);
        cout << countBBT(h, dp) << endl;
    }
}

// 595. Binary Tree Longest Consecutive Sequence (LintCode)
int maxLen = 0;
void longestConsecutive(TreeNode *root, int len, int prev)
{
    if (root == nullptr)
    {
        maxLen = max(maxLen, len);
        return;
    }

    if (prev + 1 == root->val)
    {
        longestConsecutive(root->left, len + 1, root->val);
        longestConsecutive(root->right, len + 1, root->val);
    }

    else
    {
        maxLen = max(maxLen, len);
        longestConsecutive(root->left, 1, root->val);
        longestConsecutive(root->right, 1, root->val);
    }
}
int longestConsecutive(TreeNode *root)
{
    maxLen = 0;
    longestConsecutive(root, 1, root->val);
    return maxLen;
}

// 563. Binary Tree Tilt
int tilt;
int findtilt_(TreeNode *root)
{
    if (root == nullptr)
        return 0;

    int leftSum = findtilt_(root->left);
    int rightSum = findtilt_(root->right);

    tilt += abs(leftSum - rightSum);

    return leftSum + rightSum + root->val;
}
int findTilt(TreeNode *root)
{
    tilt = 0;
    findtilt_(root);
    return tilt;
}

// Binary Tree to BST
void traversal(Node *root, vector<int> &inorder)
{
    if (root == nullptr)
        return;

    traversal(root->left, inorder);
    inorder.push_back(root->data);
    traversal(root->right, inorder);
}
int idx = 0;
void makeBST(Node *root, vector<int> &inorder)
{
    if (root == nullptr)
        return;

    makeBST(root->left, inorder);
    root->data = inorder[idx++];
    makeBST(root->right, inorder);
}
Node *binaryTreeToBST(Node *root)
{
    idx = 0;
    vector<int> inorder;

    traversal(root, inorder);
    sort(inorder.begin(), inorder.end());
    makeBST(root, inorder);

    return root;
}

// Change of Key in BST
struct node
{
    int key;
    struct node *left;
    struct node *right;
    node(int x)
    {
        key = x;
        left = NULL;
        right = NULL;
    }
};
struct node *deleteKey(struct node *root, int val)
{
    if (root == nullptr)
        return root;

    if (val < root->key)
        root->left = deleteKey(root->left, val);
    else if (val > root->key)
        root->right = deleteKey(root->right, val);

    else
    {
        if (root->left == nullptr || root->right == nullptr)
            return root->left == nullptr ? root->right : root->left;

        struct node *leftNode = root->left;
        while (leftNode->right != nullptr)
            leftNode = leftNode->right;

        int maxInLeft = leftNode->key;
        root->key = maxInLeft;
        root->left = deleteKey(root->left, maxInLeft);
    }

    return root;
}
struct node *insertKey(struct node *root, int val)
{
    if (root == nullptr)
    {
        return new struct node(val);
    }

    if (val < root->key)
        root->left = insertKey(root->left, val);
    else
        root->right = insertKey(root->right, val);

    return root;
}
struct node *changeKey(struct node *root, int oldVal, int newVal)
{
    deleteKey(root, oldVal);
    insertKey(root, newVal);
    return root;
}

// Check if subtree
bool isSubtree_(TreeNode *root, TreeNode *sub)
{
    if (root == nullptr || sub == nullptr)
    {
        if (sub == nullptr && root == nullptr)
            return true;
        return false;
    }

    return root->val == sub->val && isSubtree_(root->left, sub->left) && isSubtree_(root->right, sub->right);
}
bool isSubtree(TreeNode *root, TreeNode *sub)
{
    return root != nullptr && (isSubtree_(root, sub) || isSubtree(root->left, sub) || isSubtree(root->right, sub));
}

// Brothers From Different Roots
int pairs = 0;
bool search(Node *root, int val)
{
    if (root == nullptr)
        return false;

    if (root->data > val)
        return search(root->left, val);
    else if (root->data < val)
        return search(root->right, val);
    else
        return true;
}
void countPairs_01(Node *root1, Node *root2, int x)
{
    if (root1 == nullptr || root2 == nullptr)
        return;

    if (search(root2, x - (root1->data)))
        pairs++;
    countPairs_01(root1->left, root2, x);
    countPairs_01(root1->right, root2, x);
}
int countPairs(Node *root1, Node *root2, int x)
{
    pairs = 0;
    countPairs_01(root1, root2, x);
    return pairs;
}

// 993. Cousins in Binary Tree
int xdepth, ydepth;
TreeNode *xPar, *yPar;
void find(TreeNode *root, int x, int y, int level, TreeNode *par)
{
    if (root == nullptr)
        return;

    if (root->val == x)
    {
        xdepth = level;
        xPar = par;
    }
    else if (root->val == y)
    {
        ydepth = level;
        yPar = par;
    }

    find(root->left, x, y, level + 1, root);
    find(root->right, x, y, level + 1, root);
}
bool isCousins(TreeNode *root, int x, int y)
{
    find(root, x, y, 0, nullptr);
    return xdepth == ydepth && xPar != yPar;
}

// Children Sum Parent
int isSumProperty(Node *root)
{
    if (root == nullptr)
        return 1;
    if (root->left == nullptr && root->right == nullptr)
        return 1;

    int leftVal = root->left != nullptr ? root->left->data : 0;
    int rightVal = root->right != nullptr ? root->right->data : 0;

    int leftRes = isSumProperty(root->left);
    int rightRes = isSumProperty(root->right);

    if (root->data == leftVal + rightVal && leftRes == 1 && rightRes == 1)
        return 1;
    return 0;
}

// Check whether BST contains Dead End
bool search(Node *root, int val)
{
    if (root == nullptr)
        return false;

    if (root->data > val)
        return search(root->left, val);
    else if (root->data < val)
        return search(root->right, val);
    else
        return true;
}
bool isDeadEnd(Node *root, Node *node)
{
    if (node == nullptr)
        return false;
    if (node->left == nullptr && node->right == nullptr)
    {
        int leafVal = node->data;
        if (leafVal == 1)
            return search(root, leafVal + 1);
        return search(root, leafVal - 1) && search(root, leafVal + 1);
    }

    return isDeadEnd(root, node->left) || isDeadEnd(root, node->right);
}
bool isDeadEnd(Node *root)
{
    return isDeadEnd(root, root);
}

// Check if Tree is Isomorphic
bool isIsomorphic(Node *root1, Node *root2)
{
    if (root1 == nullptr && root2 == nullptr)
        return true;
    if (root1 == nullptr || root2 == nullptr)
        return false;

    return (root1->data == root2->data) &&
           ((isIsomorphic(root1->left, root2->left) && isIsomorphic(root1->right, root2->right)) ||
            (isIsomorphic(root1->left, root2->right) && isIsomorphic(root1->right, root2->left)));
}

// 900. Closest Binary Search Tree Value
int closestValue(TreeNode *root, double target)
{
    int ans = root->val;
    TreeNode *curr = root;

    while (curr != nullptr)
    {
        //fabs (float abs) for abs of double values
        if (fabs(curr->val - target) < fabs(ans - target))
            ans = curr->val;

        curr = (curr->val < target) ? curr->right : curr->left;
    }

    return ans;
}

// 854. Closest Leaf in a Binary Tree
vector<TreeNode *> rootToNodePath(TreeNode *root, int node)
{
    if (root == nullptr)
        return {};
    if (root->val == node)
    {
        vector<TreeNode *> base;
        base.push_back(root);
        return base;
    }

    vector<TreeNode *> leftPath = rootToNodePath(root->left, node);
    if (leftPath.size() != 0)
    {
        leftPath.push_back(root);
        return leftPath;
    }

    vector<TreeNode *> rightPath = rootToNodePath(root->right, node);
    if (rightPath.size() != 0)
    {
        rightPath.push_back(root);
        return rightPath;
    }

    return {};
}
int dist = 1e8, leafVal;
bool inSubtree = false, isLeft = false;
void kdown(TreeNode *node, TreeNode *blockNode, int level, bool lor, bool subtree)
{
    if (node == nullptr || node == blockNode)
        return;

    if (node->left == nullptr && node->right == nullptr)
    {
        if (level < dist)
        {
            leafVal = node->val;
            dist = level;
            isLeft = lor;
            inSubtree = subtree;
        }
        else if (level == dist && inSubtree == false)
        {
            if (isLeft == false && lor == true)
            {
                leafVal = node->val;
                dist = level;
                isLeft = lor;
                inSubtree = subtree;
            }
        }
    }

    kdown(node->left, blockNode, level + 1, true, subtree);
    kdown(node->right, blockNode, level + 1, false, subtree);
}
void getAllLeaves(TreeNode *root, int target)
{
    vector<TreeNode *> path = rootToNodePath(root, target);
    TreeNode *blockNode = nullptr;

    kdown(path[0], blockNode, 0, -1, true);
    for (int i = 1; i < path.size(); i++)
    {
        kdown(path[i], blockNode, i, -1, false);
        blockNode = path[i];
    }
}
int findClosestLeaf(TreeNode *root, int k)
{
    if (root->left == nullptr && root->right == nullptr)
        return root->val;

    getAllLeaves(root, k);
    return leafVal;
}

// Closest Neighbor in BST
int findMaxForN(Node *root, int N, int size)
{
    Node *curr = root;
    int ans = -1e8;

    while (curr != nullptr)
    {
        if (curr->data <= N && N - curr->data < N - ans)
            ans = curr->data;
        if (curr->data < N)
            curr = curr->right;
        else
            curr = curr->left;
    }

    if (ans == -1e8)
        return -1;
    return ans;
}

// Delete nodes greater than k
Node *deleteNode(Node *root, int key)
{
    if (root == nullptr)
        return root;

    if (key > root->data)
        root->right = deleteNode(root->right, key);
    else if (key <= root->data)
    {
        root = deleteNode(root->left, key);
        return root;
    }

    return root;
}

//116. Populating Next Right Pointers in Each Node
Node *connect(Node *root)
{
    if (root == nullptr)
        return root;

    //curr to traverse from one level to other
    Node *curr = root;
    while (curr->left != nullptr)
    {
        //temp variable to traverse through current level
        Node *temp = curr;
        while (temp->next != nullptr)
        {
            temp->left->next = temp->right;
            temp->right->next = temp->next->left;
            temp = temp->next;
        }
        temp->left->next = temp->right;
        curr = curr->left;
    }

    return root;
}

// 226. Invert Binary Tree
TreeNode *invertTree(TreeNode *root)
{
    if (root == nullptr)
        return root;

    TreeNode *leftChild = invertTree(root->left);
    TreeNode *rightChild = invertTree(root->right);

    root->left = rightChild;
    root->right = leftChild;

    return root;
}

// 99. Recover Binary Search Tree
TreeNode *node1, *node2, *pre;
void recoverTree_(TreeNode *root)
{
    if (root == nullptr)
        return;

    recoverTree_(root->left);

    if (pre != nullptr && root->val < pre->val)
    {
        if (node1 == nullptr)
            node1 = pre;
        node2 = root;
    }
    pre = root;

    recoverTree_(root->right);
}
void recoverTree(TreeNode *root)
{
    recoverTree_(root);

    int temp = node1->val;
    node1->val = node2->val;
    node2->val = temp;
}

// 1372. Longest ZigZag Path in a Binary Tree
//<leftPath,rightPath>
int maxLen = 0;
pair<int, int> longestZigZag_(TreeNode *root)
{
    if (root == nullptr)
        return {0, 0};

    pair<int, int> leftPath = longestZigZag_(root->left);
    pair<int, int> rightPath = longestZigZag_(root->right);

    maxLen = max(maxLen, max(leftPath.second, rightPath.first));

    return {leftPath.second + 1, rightPath.first + 1};
}
int longestZigZag(TreeNode *root)
{
    maxLen = 0;
    longestZigZag_(root);
    return maxLen;
}

// 101. Symmetric Tree
bool isSymmetric_01(TreeNode *root1, TreeNode *root2)
{
    if (root1 == nullptr && root2 == nullptr)
        return true;
    if (root1 == nullptr || root2 == nullptr)
        return false;
    if (root1->val != root2->val)
        return false;

    return isSymmetric_01(root1->left, root2->right) && isSymmetric_01(root1->right, root2->left);
}
bool isSymmetric(TreeNode *root)
{
    return isSymmetric_01(root, root);
}

// 951. Flip Equivalent Binary Trees
bool isSymmetric_01(TreeNode *root1, TreeNode *root2)
{
    if (root1 == nullptr && root2 == nullptr)
        return true;
    if (root1 == nullptr || root2 == nullptr)
        return false;
    if (root1->val != root2->val)
        return false;

    return (isSymmetric_01(root1->left, root2->right) && isSymmetric_01(root1->right, root2->left)) || (isSymmetric_01(root1->left, root2->left) && isSymmetric_01(root1->right, root2->right));
}
bool flipEquiv(TreeNode *root1, TreeNode *root2)
{
    return isSymmetric_01(root1, root2);
}

// 958. Check Completeness of a Binary Tree
bool isCompleteTree(TreeNode *root)
{
    //base case-> single node
    if (root->left == nullptr && root->right == nullptr)
        return true;
    TreeNode *curr = root;
    list<TreeNode *> que;
    que.push_back(root);

    while (que.size() != 0)
    {
        int size = que.size();
        int n = size;
        list<TreeNode *> prevlvl = que;
        while (size-- > 0)
        {
            TreeNode *rnode = que.front();
            que.pop_front();

            if (rnode->left != nullptr)
                que.push_back(rnode->left);
            if (rnode->right != nullptr)
                que.push_back(rnode->right);
        }

        //check last level and last second level
        if (que.front()->left == nullptr)
        {
            //check last level for all leaf nodes
            while (que.size() != 0)
            {
                TreeNode *rnode = que.front();
                que.pop_front();
                if (rnode->left != nullptr || rnode->right != nullptr)
                    return false;
            }

            //check second last level to check if last is leftmost (only after last level is found)
            TreeNode *prev = prevlvl.front()->left;
            while (prevlvl.size() != 0)
            {
                TreeNode *rnode = prevlvl.front();
                prevlvl.pop_front();

                if (prev == nullptr && rnode->left != nullptr)
                    return false;
                prev = rnode->left;
                if (prev == nullptr && rnode->right != nullptr)
                    return false;
                prev = rnode->right;
            }
        }

        //check no. of nodes at each level, each level is completely filled(except last level)
        else if (que.size() != n * 2)
            return false;
    }

    return true;
}

// 864. Equal Tree Partition
bool res = false;
TreeNode *oRoot;
int findSum(TreeNode *root)
{
    if (root == nullptr)
        return 0;

    return findSum(root->left) + findSum(root->right) + root->val;
}
int findEdge(TreeNode *root, int totalSum)
{
    if (root == nullptr || res)
        return 0;

    int leftSum = findEdge(root->left, totalSum);
    int rightSum = findEdge(root->right, totalSum);

    int mySum = leftSum + rightSum + root->val;
    if (mySum == totalSum - mySum)
    {
        //cannot break the tree at original root because at root mySum==totalSum
        if (oRoot != root)
        {
            res = true;
            return 0;
        }
    }

    return mySum;
}
bool checkEqualTree(TreeNode *root)
{
    oRoot = root;
    int totalSum = findSum(root);
    if (totalSum % 2 != 0)
        return false;
    findEdge(root, totalSum);
    return res;
}

// 513. Find Bottom Left Tree Value
int findBottomLeftValue(TreeNode *root)
{
    list<TreeNode *> que;
    que.push_back(root);
    int ans;
    while (que.size() != 0)
    {
        ans = que.front()->val;
        int size = que.size();
        while (size-- > 0)
        {
            TreeNode *rnode = que.front();
            que.pop_front();
            if (rnode->left != nullptr)
                que.push_back(rnode->left);
            if (rnode->right != nullptr)
                que.push_back(rnode->right);
        }
    }
    return ans;
}

// 652. Find Duplicate Subtrees
//serialize in preorder, because inorder gets same string for skewed left and right trees
string serializeSubtrees(TreeNode *root, unordered_map<string, int> &mp, vector<TreeNode *> &ans)
{
    if (root == nullptr)
        return "#";

    string str = "";
    str += to_string(root->val);
    str += ',';
    str += serializeSubtrees(root->left, mp, ans);
    str += ',';
    str += serializeSubtrees(root->right, mp, ans);

    mp[str]++;
    if (mp[str] == 2)
        ans.push_back(root);

    return str;
}
vector<TreeNode *> findDuplicateSubtrees(TreeNode *root)
{
    unordered_map<string, int> mp;
    vector<TreeNode *> ans;
    serializeSubtrees(root, mp, ans);
    return ans;
}

// 1123. Lowest Common Ancestor of Deepest Leaves
TreeNode *LCA;
TreeNode *firstLeafNode;
int height(TreeNode *root)
{
    if (root == nullptr)
        return -1;
    return max(height(root->left), height(root->right)) + 1;
}
bool findLCA(TreeNode *root, int depth, int maxDepth)
{
    if (root == nullptr)
        return false;
    if (depth == maxDepth)
    {
        if (firstLeafNode == nullptr)
            firstLeafNode = root;
        return true;
    }

    //cannot return after first time an lca has been found as the lca closest to root is the ans as many paths can give deepest leaves
    bool leftDone = findLCA(root->left, depth + 1, maxDepth);
    bool rightDone = findLCA(root->right, depth + 1, maxDepth);

    //the very last root to be set as the lca will be the ans, as many roots can have both leftDone==true and rightDone==true
    //but there may be more deepest leaves in other subtrees
    if (leftDone && rightDone)
        LCA = root;

    return leftDone || rightDone;
}
TreeNode *lcaDeepestLeaves(TreeNode *root)
{
    //find height of tree -> to know the max depth -> so that deepest leaves can be found
    int maxDepth = height(root);

    //find lca
    findLCA(root, 0, maxDepth);

    //if lca has been found
    if (LCA != nullptr)
        return LCA;

    //if lca has not been found -> only one deepest leaf exists -> it is the lca
    return firstLeafNode;
}

// 979. Distribute Coins in Binary Tree
int moves;
int distributeCoins_(TreeNode *root)
{
    if (root == nullptr)
        return 0;

    int leftCoins = distributeCoins_(root->left);
    int rightCoins = distributeCoins_(root->right);

    moves += abs(leftCoins) + abs(rightCoins);

    return root->val + leftCoins + rightCoins - 1;
}
int distributeCoins(TreeNode *root)
{
    moves = 0;
    distributeCoins_(root);
    return moves;
}

// 1145. Binary Tree Coloring Game
int countNodes(TreeNode *root)
{
    if (root == nullptr)
        return 0;

    return countNodes(root->left) + countNodes(root->right) + 1;
}
TreeNode *findNode(TreeNode *root, int data)
{
    if (root == nullptr)
        return root;
    if (root->val == data)
        return root;

    TreeNode *leftFind = findNode(root->left, data);
    if (leftFind != nullptr)
        return leftFind;

    TreeNode *rightFind = findNode(root->right, data);
    if (rightFind != nullptr)
        return rightFind;

    return nullptr;
}
bool btreeGameWinningMove(TreeNode *root, int n, int x)
{
    int totalNodes = n;
    TreeNode *xNode = findNode(root, x);
    int leftCount = countNodes(xNode->left);
    int rightCount = countNodes(xNode->right);
    int parCount = n - leftCount - rightCount - 1;

    return leftCount > n / 2 || rightCount > n / 2 || parCount > n / 2;
}

// 1325. Delete Leaves With a Given Value
TreeNode *removeLeafNodes(TreeNode *root, int target)
{
    if (root == nullptr)
        return nullptr;
    if (root->val == target && root->left == nullptr && root->right == nullptr)
        return nullptr;

    root->left = removeLeafNodes(root->left, target);
    root->right = removeLeafNodes(root->right, target);

    if (root->val == target && root->left == nullptr && root->right == nullptr)
        return nullptr;

    return root;
}

// 1104. Path In Zigzag Labelled Binary Tree
//TLE(if we construct tree)
bool findPath(TreeNode *root, int tar, vector<int> &path)
{
    if (root == nullptr)
        return false;
    if (root->val == tar)
    {
        path.push_back(root->val);
        return true;
    }

    bool leftFind = findPath(root->left, tar, path);
    if (leftFind)
    {
        path.push_back(root->val);
        return true;
    }

    bool rightFind = findPath(root->right, tar, path);
    if (rightFind)
    {
        path.push_back(root->val);
        return true;
    }

    return false;
}
vector<int> pathInZigZagTree(int label)
{
    list<TreeNode *> que;
    TreeNode *root = new TreeNode(1);
    que.push_back(root);

    //0->fill from front , 1->fill form back
    int prevMax = 1, prevSize = 1, flag = 1;
    while (que.size() != 0)
    {
        int size = que.size();

        //set before each level
        int lb = prevMax + 1;
        prevSize *= 2;
        int ub = lb + prevSize - 1;
        int firstele = flag == 0 ? lb : ub;

        while (size-- > 0)
        {
            TreeNode *rnode = que.front();
            que.pop_front();

            if (flag == 0)
            {
                rnode->left = new TreeNode(firstele);
                rnode->right = new TreeNode(firstele + 1);
                firstele += 2;
            }
            else
            {
                rnode->left = new TreeNode(firstele);
                rnode->right = new TreeNode(firstele - 1);
                firstele -= 2;
            }

            que.push_back(rnode->left);
            que.push_back(rnode->right);
        }

        //update after each level
        prevMax = ub;
        flag ^= 1;
        if (prevMax > label)
            break;
    }

    vector<int> path;
    findPath(root, label, path);
    reverse(path.begin(), path.end());
    return path;
}

//instead we use properties of full binary tree (O(logn))
vector<int> pathInZigZagTree(int label)
{
    vector<int> path;
    int level = 1, temp = label;

    //find the total levels
    while (temp != 1)
    {
        temp /= 2;
        level++;
    }

    //find path from last level to root
    while (label != 1)
    {
        path.push_back(label);
        int maxVal = pow(2, level) - 1;
        int minVal = pow(2, level - 1);
        int par = (minVal + maxVal - label) / 2;
        label = par;
        level--;
    }
    path.push_back(1);
    reverse(path.begin(), path.end());
    return path;
}

// 1448. Count Good Nodes in Binary Tree
int goodNodes(TreeNode *node, int maxVal)
{
    if (node == nullptr)
        return 0;

    //get count of goodNodes on left and right
    int leftCount = goodNodes(node->left, max(maxVal, node->val));
    int rightCount = goodNodes(node->right, max(maxVal, node->val));

    //if this node is also a good node
    if (node->val >= maxVal)
        return leftCount + rightCount + 1;

    //if this node is not a good node
    return leftCount + rightCount;
}
int goodNodes(TreeNode *root)
{
    return goodNodes(root, -1e8);
}

// 1443. Minimum Time to Collect All Apples in a Tree
int minTime(int n, vector<vector<int>> &edges, vector<bool> &hasApple)
{
}

// ===============================================================================================================================
void solve()
{
    pathInZigZagTree(14);
}
int main()
{
    solve();
    return 0;
}