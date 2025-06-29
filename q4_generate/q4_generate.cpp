﻿#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>

using namespace std;

int main() {
    const int lp = 9;                     // 网格边长：9x9
    const int nnp = lp * lp;              // 总节点数 81
    const int nel = (lp - 1) * (lp - 1);  // 总单元数 64
    const int ndof = 2;                   // 每节点自由度

    vector<double> x_coords, y_coords;

    // 生成坐标：左边加高，形成梯形
    for (int i = 0; i < lp; ++i) {
        double xi = 2.0 * i / (lp - 1);       // linspace(0,2)
        double yi0 = 0.25 * xi + 0.1;         // 修改底边高度，加高左边

        for (int j = 0; j < lp; ++j) {
            double yj = yi0 + j * (1 - yi0) / (lp - 1); // linspace(y0(i), 1)
            x_coords.push_back(xi);
            y_coords.push_back(yj);
        }
    }

    // 生成 IEN 单元连接数组
    vector<vector<int>> IEN;
    int rowcount = 0;
    for (int e = 0; e < nel; ++e) {
        int n1 = e + rowcount;
        int n2 = n1 + 1;
        int n3 = n1 + lp + 1;
        int n4 = n1 + lp;
        IEN.push_back({ n1 + 1, n2 + 1, n3 + 1, n4 + 1 });
        if ((e + 1) % (lp - 1) == 0) rowcount++;
    }

    // 打开输出文件（仍为原文件名）
    ofstream fout("q4-example.dat");
    if (!fout) {
        cerr << "无法创建输出文件 q4-example.dat" << endl;
        return -1;
    }

    // 写入文件头
    fout << "Q4_example\n";
    fout << nnp << " 1 1 1\n";

    // 写入节点坐标和边界条件
    for (int i = 0; i < nnp; ++i) {
        int node_id = i + 1;
        int bc_x = 0;  // 不固定任何x方向
        int bc_y = (fabs(x_coords[i]) < 1e-8) ? 1 : 0;  // 仅左侧固定y方向

        fout << node_id << " " << bc_x << " " << bc_y << " 1 "
            << scientific << setprecision(5)
            << x_coords[i] << " " << y_coords[i] << " 0\n";
    }

    // 写入集中载荷：节点81 y方向 -1000N
    fout << "1 1\n";
    fout << nnp << " 2 -1000.0\n";
    fout << "2 " << nel << " 1\n";  // 材料编号

    // 材料属性：E, nu, t
    fout << "1 " << scientific << setprecision(1)
        << 2.0e11 << " 0.3 1.0\n";

    // 写入单元连接
    for (int i = 0; i < nel; ++i) {
        fout << i + 1 << " "
            << IEN[i][0] << " "
            << IEN[i][1] << " "
            << IEN[i][2] << " "
            << IEN[i][3] << " "
            << "1\n";
    }

    fout.close();
    cout << "文件已生成: q4-example.dat" << endl;
    return 0;
}




// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门使用技巧: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
