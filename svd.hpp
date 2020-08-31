/*
 * KG: 已经验证这段代码的结果是正确的
 * TODO: 很多计算可以用CMathUtils加速优化
 * 本程序在linux g++下编译通过
 * bool svd(vector<vector<double> > A, int K, vector<vector<double> > &U, vector<double> &S, vector<vector<double> > &V);
 * A: 输入待分解矩阵
 * K: 输入，取前K大奇异值及奇异向量
 * U[0],U[1],...,U[K-1]: 前K大奇异值对应的左奇异向量
 * S[0],S[1],...,S[K-1]: 前K大奇异值 S[0]>=S[1]>=...>=S[K-1]
 * V[0],V[1],...,V[K-1]: 前K大奇异值对应的右奇异向量
*/
#ifndef SVD_HPP
#define SVD_HPP
#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <vector>
#include "sys_util.hpp"
#include "string_util.hpp"

namespace audio_engine {
    using namespace std;
    const int MAX_ITER = 100000;
    const double eps = 0.0000001;

    double get_norm(double *x, int n) {
        double r = 0;
        for (int i = 0; i < n; i++)
            r += x[i] * x[i];
        return sqrt(r);
    }

    double normalize(double *x, int n) {
        double r = get_norm(x, n);
        if (r < eps)
            return 0;
        for (int i = 0; i < n; i++)
            x[i] /= r;
        return r;
    }

    inline double product(double *a, double *b, int n) {
        double r = 0;
        for (int i = 0; i < n; i++)
            r += a[i] * b[i];
        return r;
    }

    void orth(double *a, double *b, int n) {//|a|=1
        double r = product(a, b, n);
        for (int i = 0; i < n; i++)
            b[i] -= r * a[i];

    }

    bool
    svd(vector<vector<double> > A, int K, vector<vector<double> > &U, vector<double> &S, vector<vector<double> > &V) {
        int M = A.size();
        int N = A[0].size();
        U.clear();
        V.clear();
        S.clear();
        S.resize(K, 0);
        U.resize(K);
        for (int i = 0; i < K; i++)
            U[i].resize(M, 0);
        V.resize(K);
        for (int i = 0; i < K; i++)
            V[i].resize(N, 0);


        srand(time(0));
        double *left_vector = new double[M];
        double *next_left_vector = new double[M];
        double *right_vector = new double[N];
        double *next_right_vector = new double[N];
        int col = 0;
        for (int col = 0; col < K; col++) {
            double diff = 1;
            double r = -1;
            while (1) {
                for (int i = 0; i < M; i++)
                    left_vector[i] = (float) rand() / RAND_MAX;
                if (normalize(left_vector, M) > eps)
                    break;
            }

            for (int iter = 0; diff >= eps && iter < MAX_ITER; iter++) {
                memset(next_left_vector, 0, sizeof(double) * M);
                memset(next_right_vector, 0, sizeof(double) * N);
                for (int i = 0; i < M; i++)
                    for (int j = 0; j < N; j++)
                        next_right_vector[j] += left_vector[i] * A[i][j];

                r = normalize(next_right_vector, N);
                if (r < eps) break;
                for (int i = 0; i < col; i++)
                    orth(&V[i][0], next_right_vector, N);
                normalize(next_right_vector, N);

                for (int i = 0; i < M; i++)
                    for (int j = 0; j < N; j++)
                        next_left_vector[i] += next_right_vector[j] * A[i][j];
                r = normalize(next_left_vector, M);
                if (r < eps) break;
                for (int i = 0; i < col; i++)
                    orth(&U[i][0], next_left_vector, M);
                normalize(next_left_vector, M);
                diff = 0;
                for (int i = 0; i < M; i++) {
                    double d = next_left_vector[i] - left_vector[i];
                    diff += d * d;
                }

                memcpy(left_vector, next_left_vector, sizeof(double) * M);
                memcpy(right_vector, next_right_vector, sizeof(double) * N);
            }
            if (r >= eps) {
                S[col] = r;
                memcpy((char *) &U[col][0], left_vector, sizeof(double) * M);
                memcpy((char *) &V[col][0], right_vector, sizeof(double) * N);
            } else {
                cout << r << endl;
                break;
            }
        }
        delete[] next_left_vector;
        delete[] next_right_vector;
        delete[] left_vector;
        delete[] right_vector;

        return true;
    }

    void getRotateAngle(const std::vector<std::vector<double> > &matrix, const int k, double &angle, double &percent) {
        // 计算二维平面上的点的分布大致方向，以及集中程度。matrix中每个分量是一个二维向量，对应二维平面上的一个顶点坐标
        vector<vector<double> > U;
        vector<double> S;
        vector<vector<double> > V;
        svd(matrix, k, U, S, V);

        angle = asin(V[0][0]) * 180.0 / M_PI;
        percent = 1.0 * S[0] / (S[0] + S[1]);
    }

    std::vector<std::vector<double>> ReadMatrixFile(std::string fpath) {
        std::vector<std::vector<double>> matrix;
        const char *fileName = fpath.c_str();
        string fname, ext, fdir;
        cb_tools::splitFilePath(fpath, fname, ext, fdir);
        vector<string> allstrs;
        cb_tools::getStringLinesFromFile(fpath, allstrs);
        for (int i = 0; i < allstrs.size(); i++) {
//        allstrs[i] = cb_tools::trimAll(allstrs[i], " ");
            vector<string> ele;
            cb_tools::split(allstrs[i], " ", ele);
            if (ele.size() < 9) {
                continue;
            }
            std::vector<double> tmp;
            for (auto &item : ele) {
                int t = cb_tools::string2int(item);
                if (t > 0) {
                    tmp.push_back(t);
                }
            }
            if (tmp.size() == 2) {
                matrix.push_back(tmp);
            }
        }
        return matrix;
    }

//int main() {
//    //读入二维矩阵
//    std::vector<std::vector<double>> A = ReadMatrixFile("/Users/ailab002/Desktop/claw_asset/data/matrix.txt");
//    int m = A.size();
//    int n = 2;
//    int k = 2;
//
//    cout << "A=" << endl;
//    for (int i = 0; i < A.size(); i++) {
//        for (int j = 0; j < A[i].size(); j++) {
//            cout << setw(12) << A[i][j] << ' ';
//        }
//        cout << endl;
//    }
//    cout << endl;
//
//    double angle = 0;
//    double percent = 0;
//    getRotateAngle(A, k, angle, percent);
//    std::cout << "angle: " << angle << endl;
//    std::cout << "percent: " << percent << endl;
//
//    vector<vector<double> > U;
//    vector<double> S;
//    vector<vector<double> > V;
//    svd(A, k, U, S, V);
//    cout << "U=" << endl;
//    for (int i = 0; i < U[0].size(); i++) {
//        for (int j = 0; j < U.size(); j++) {
//            cout << setw(12) << U[j][i] << ' ';
//        }
//        cout << endl;
//    }
//    cout << endl;
//    cout << "S=" << endl;
//    for (int i = 0; i < S.size(); i++) {
//        cout << setw(7) << S[i] << ' ';
//    }
//    cout << endl;
//    cout << "V=" << endl;
//    for (int i = 0; i < V[0].size(); i++) {
//        for (int j = 0; j < V.size(); j++) {
//            cout << setw(12) << V[j][i] << ' ';
//        }
//        cout << endl;
//    }
//    return 0;
//}
}
#endif