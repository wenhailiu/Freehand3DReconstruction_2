#ifndef UTILITIES
#define UTILITIES

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <chrono>

namespace Utilities{
	
template<typename T> 
extern void readFromBin(T *Output, int Num_Elements, const std::string FILENAME);

template<typename T>
extern void writeToBin(T *Output, int Num_Elements, const std::string FILENAME);

template<typename T>
void MatrixDeserialization(const std::vector<std::vector<T>>& InputMatrix, T* OutputArray){
	for(int row_it = 0; row_it < InputMatrix.size(); ++row_it){
		for(int col_it = 0; col_it < InputMatrix[0].size(); ++col_it){
			OutputArray[col_it + row_it * InputMatrix[0].size()] = InputMatrix[row_it][col_it];
		}
	}
}

extern bool MatrixS4X4Invert(const float m[16], float invOut[16]);
extern void MatrixMatrixS4X4Multiply(const float m_i_1[16], const float m_i_2[16], float m_o[16]);
extern void MatrixPointS4X4Multiply(const float m_i[16], const float p_i[4], float p_o[4]);
extern void MatrixS4X4Identity(std::vector<float>& InputMatrix);
extern void MatrixS4X4Identity(std::vector<std::vector<float>>& InputMatrix);

class MyTimer{
public:
    MyTimer();
    void tic();
    void toc();
    double Duration(std::ostream& os);
    double Duration();

private:
    std::chrono::time_point<std::chrono::_V2::system_clock, std::chrono::nanoseconds> begin;
    std::chrono::time_point<std::chrono::_V2::system_clock, std::chrono::nanoseconds> end;

    bool BeginSet;
    bool EndSet;
    double duration;
};

std::vector<float> MatrixVector2D_Deserializer(const std::vector<std::vector<float>>& inputMatrix);
void MatrixPrint(const std::vector<float> InputMatrix);
void MatrixPrint(const float *InputMatrix);

}//End namespace: Utilities



#endif