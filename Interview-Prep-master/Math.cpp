#include <iostream>
#include <vector>
#include <list>
#include <algorithm>

using namespace std;

// 171. Excel Sheet Column Number
/*
For eg:
A=1
ZY  = 26*26+25
CY = 26*3 + 25
XDV = (26*24 + 4)*26 + 22 , because int(D)=4, int(V)=21
*/
int titleToNumber(string s)
{
    int colNum = 0;

    for (int i = 0; i < s.size(); i++)
    {
        int num = s[i] - 'A' + 1;
        colNum = colNum * 26 + num;
    }

    return colNum;
}

void solve()
{
}
int main()
{
    solve();
    return 0;
}