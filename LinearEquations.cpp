#include "Matrix.h"

#include <iostream>
#include <fstream>
#include <utility>
#include <chrono>
#include <tuple>

constexpr std::size_t NUMBER_OF_RESULTS_IN_ROW{ 5U }, F{ 8U };
constexpr double maxNorm{ 10e-9 };

enum class Method {
    Jacobi,
    Gauss,
    LU_Decomposition
};

double* residual(Matrix* A, double* b, double* x, const std::size_t N) {

    double* r = (*A) * x;

    for (std::size_t i{}; i < N; i++) 
        r[i] -= b[i];
    
    return r;
}

double norm(double* v, const std::size_t N) {

    double n{};

    for (std::size_t i{}; i < N; i++) {
        const double currVal{ v[i] };
        n += currVal * currVal;
    }

    return sqrt(n);
}

double* forwardSubstitution(Matrix* A, double* b, const std::size_t N) { //b = U * r, A = -(D + L)

    double* result = new double[N];
    std::fill(result, result + N, 0);

    result[0] = b[0] / A->A[0][0];

    for (int i{1}; i < N; i++) {

        double wart = b[i];

        for (int j{}; j < i; j++) {
            wart -= A->A[i][j] * result[j];
        }

        wart /= A->A[i][i];
        result[i] = wart;
    }

    return result; //wynik to b * A^-1
}

double* backwardSubstitution(Matrix* A, double* y, const std::size_t N) { //y is from forward, A = -(D + U)

    double* result = new double[N];
    std::fill(result, result + N, 0);

    for (int i = N - 1; i >= 0; i--) {

        double wart = y[i];

        for (int j{ i + 1 }; j < N; j++) {
            wart -= A->A[i][j] * result[j];
        }

        wart /=  A->A[i][i];
        result[i] = wart;
    }

    return result;
} //wynik y * A^-1


std::tuple<double, int, double> jacobi(Matrix* A, double* b, double* x, const std::size_t N, std::vector<double>& norms, std::vector<double>& iterations) {

    unsigned int iterationCount{};
    double calcTime{}, rowSum{}, normVal{ 1 };

    double* res = new double[N];
    double* const xPrev = new double[N];

    std::fill(xPrev, xPrev + N, 1); //They have to start from 1

    const auto startTime{ std::chrono::system_clock::now() };

    while (normVal > maxNorm) {

        for (std::size_t i{}; i < N; i++) {

            rowSum = 0;

            for (std::size_t j{}; j < N; j++) 
                if (j != i) //Taking values from L or U, avoiding diagonal values
                    rowSum += A->A[i][j] * xPrev[j];

            //x[i] = rowSum * (-1) / A->A[i][i] + b[i] / A->A[i][i]; //Full formula

            //rowSum is a i-th row of (L + U) * r^k, where r^k is xPrev - vector of previous calculations, and i-th row is calculated by "rowSum += A->A[i][j] * xPrev[j];".
            //We only need to take i-th row of b and rowSum vector since we are multiplying twice by diagonal matrix, the rest of the sums are 0 no matter the other rows of rowSum or b.
            x[i] = (-rowSum + b[i]) / A->A[i][i]; //Simplified formula
        }

        for (std::size_t i{}; i < N; i++) 
            xPrev[i] = x[i]; //store previous vector calculations
        
        delete[] res;

        res = residual(A, b, x, N); 
        normVal = norm(res, N);

        norms.push_back(normVal);
        iterations.push_back(iterationCount);

        iterationCount++;
    }

    const auto endTime{ std::chrono::system_clock::now() };
    calcTime = static_cast<double>( std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count() );

    delete[] res, xPrev;

    return std::make_tuple(calcTime, iterationCount, normVal);
}


std::tuple<double, int, double> gaussSeidl(Matrix* A, double* b, double* x, const std::size_t N, std::vector<double>& norms, std::vector<double>& iterations) {

    int iterationCount{};
    double calcTime{}, tempSum{}, normVal{ 1 };

    double* res = new double[N];
    double* const xPrev = new double[N];

    std::fill(xPrev, xPrev + N, 1); 

    Matrix* DL = new Matrix(N, 0, 0, 0);
    Matrix* U = new Matrix(N, 0, 0, 0);

    for (int i = 0; i < N; i++) {

        for (int j = 0; j < N; j++) {

            if (j <= i)  //taking only D or L values
                DL->A[i][j] += A->A[i][j]; //setting them up
           
            else
                U->A[i][j] += A->A[i][j]; //otherwise it's U
        }
    }

    auto startTime{ std::chrono::system_clock::now() };

    while (normVal > maxNorm) {

        double* Urk = (*U) * xPrev;

        double* firstSum = forwardSubstitution(DL, Urk, N); //first half of the sum equation, fs makes ^-1 out of DL - (D + L) vector.
        double* secondSum = forwardSubstitution(DL, b, N); //second half

        for (int j = 0; j < N; j++) {
            x[j] = (-1) * firstSum[j] + secondSum[j]; //adding them up
            xPrev[j] = x[j]; 
        }

        delete[] res;

        res = residual(A, b, x, N);
        normVal = norm(res, N);

        norms.push_back(normVal);
        iterations.push_back(iterationCount);

        iterationCount++;
    }

    auto endTime{ std::chrono::system_clock::now() };
    calcTime = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count());

    delete[] res, xPrev, U->A, DL->A;
    delete U, DL;

    return std::make_tuple(calcTime, iterationCount, normVal);
}

void printResults(Method method, std::tuple<double, int, double>& additionalResults, double* x, const std::size_t N) {

    std::string_view selectedMethod{};

    switch (method) {

        using enum Method;

        case Jacobi: selectedMethod = "Jacobi"; break;
        case Gauss: selectedMethod = "Gauss"; break;
        case LU_Decomposition: selectedMethod = "LU Decomposition"; break;

    }

    std::cout << '\n' << selectedMethod << ": " << '\n';

    std::cout << "Duration in milliseconds: " << std::get<0>(additionalResults) << '\n';
    std::cout << "Iterations number: " << std::get<1>(additionalResults) << '\n';
    std::cout << "Residual vector norm: " << std::get<2>(additionalResults) << '\n' << '\n';

    std::cout << '\n' << '\n';
}

std::tuple<double, int, double> LU_decomposition(Matrix* A, double* b, double* x, const std::size_t N)
{
    double calcTime{}, normVal{};

    Matrix* U = new Matrix(*A);			//upper triangular matrix
    Matrix* L = new Matrix(N, 1, 0, 0);	//lower triangular matrix

    auto startTime{ std::chrono::system_clock::now() };

    for (int i = 0; i < N - 1; ++i)	
    {
        for (int j = i + 1; j < N; ++j)
        {
            L->A[j][i] = U->A[j][i] / U->A[i][i];
            for (int k = i; k < N; ++k)
                U->A[j][k] = U->A[j][k] - L->A[j][i] * U->A[i][k];
        }
    }
     
    double* y = new double[N];

    y = forwardSubstitution(L, b, N);
    double* resX = backwardSubstitution(U, y, N);

    for (int i = 0; i < N; i++) {
        x[i] = resX[i];
    }

    auto endTime{ std::chrono::system_clock::now() };
    calcTime = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count());

    double* res = residual(A, b, x, N);		//residual vector
    normVal = norm(res, N);

    delete[] L->A, U->A, y, res, resX;
    delete L, U;

    return std::make_tuple(calcTime, 0, normVal);
}

void addToFile(Method method, std::ofstream& plotData, std::vector<double>& sizes, std::vector<double>& calcTimes) {

    std::string_view selectedMethod{};

    switch (method) {

        using enum Method;

        case Jacobi: selectedMethod = "Jacobi"; break;
        case Gauss: selectedMethod = "Gauss"; break;
        case LU_Decomposition: selectedMethod = "LU Decomposition"; break;

    }

    plotData << selectedMethod << '\n';

    for (int i = 0; i < sizes.size(); i++) {

        plotData << " " << calcTimes[i];

    }

    plotData << '\n';

    for (int i = 0; i < sizes.size(); i++) {

        plotData << " " << sizes[i];

    } 

    plotData << '\n';
}


int main()
{
    constexpr std::size_t e{ 7U }, N{ 908U }; //908

    constexpr int a1{ 5 + e }, a2{ -1 }, a3{ -1 };    

    //Task A
    Matrix* A = new Matrix(N, a1, a2, a3);

    double* b = new double[N];
    double* x = new double[N];
    
    for (std::size_t i{}; i < N; i++) 
        b[i] = sin(i * (F + 1));

    //Task B
    std::vector<double> normsJacobi{}, normsGauss{};
    std::vector<double> iterationsJacobi{}, iterationsGauss{};

    std::tuple<double, int, double> additionalResults = jacobi(A, b, x, N, normsJacobi, iterationsJacobi);
    printResults(Method::Jacobi, additionalResults, x, N);

    additionalResults = gaussSeidl(A, b, x, N, normsGauss, iterationsGauss);
    printResults(Method::Gauss, additionalResults, x, N);

    delete[] A->A;
    delete A;

    //Task C
    constexpr int a4 = 3;
    A = new Matrix(N, a4, a2, a3);

    //Task D
    additionalResults = LU_decomposition(A, b, x, N);
    printResults(Method::LU_Decomposition, additionalResults, x, N);

    //Task E
    constexpr std::size_t NUMBER_OF_N{ 5 - 3 }, N_DIFFERENCE{ 1000 }, STARTING_N{ 2000 };
    
    std::vector<double> sizes{ 100.0, 500.0, 1000.0 };

    std::size_t nToPushBack{ STARTING_N };

    for (int i = 0; i < NUMBER_OF_N; i++) {

        sizes.push_back(nToPushBack);

        if (nToPushBack > std::numeric_limits<std::size_t>::max() - N_DIFFERENCE) break;

        nToPushBack += N_DIFFERENCE;
        
    }

    std::cout << '\n' << "Additional results (Task E): " << '\n';

    delete[] A->A, b, x;
    delete A;

    std::vector<double> calcTimesJacobi{}, calcTimesGauss{}, calcTimesLU{};
    std::vector<double> normsJacobi2{}, normsGauss2{};
    std::vector<double> iterationsJacobi2{}, iterationsGauss2{};

    for (const std::size_t& i : sizes) {

        std::cout << '\n' << "For size N = " << i << '\n';

        A = new Matrix(i, a1, a2, a3);

        x = new double[i];
        b = new double[i];

        for (std::size_t j{}; j < i; j++)
            b[j] = sin(j * (F + 1));
      
        additionalResults = jacobi(A, b, x, i, normsJacobi2, iterationsJacobi2);

        calcTimesJacobi.push_back(std::get<0>(additionalResults));
        normsJacobi.push_back(std::get<2>(additionalResults));
        iterationsJacobi.push_back(std::get<1>(additionalResults));

        printResults(Method::Jacobi, additionalResults, x, i);

        additionalResults = gaussSeidl(A, b, x, i, normsGauss2, iterationsGauss2);

        calcTimesGauss.push_back(std::get<0>(additionalResults));
        normsGauss.push_back(std::get<2>(additionalResults));
        iterationsGauss.push_back(std::get<1>(additionalResults));

        printResults(Method::Gauss, additionalResults, x, i);

        additionalResults = LU_decomposition(A, b, x, i);

        calcTimesLU.push_back(std::get<0>(additionalResults));

        printResults(Method::LU_Decomposition, additionalResults, x, i);

        delete[] A->A, b, x;
        delete A;
    }

    std::ofstream plotData{};
    plotData.open("Wyniki/plotData.txt");

    addToFile(Method::Jacobi, plotData, sizes, calcTimesJacobi);
    addToFile(Method::Gauss, plotData, sizes, calcTimesGauss);
    addToFile(Method::LU_Decomposition, plotData, sizes, calcTimesLU);

    std::ofstream plotDataNorm;
    plotDataNorm.open("Wyniki/plotDataNorm.txt");

    addToFile(Method::Jacobi, plotDataNorm, iterationsJacobi, normsJacobi);
    addToFile(Method::Gauss, plotDataNorm, iterationsGauss, normsGauss);

    plotData.close();
    plotDataNorm.close();

    sizes.clear();
    sizes.shrink_to_fit();
    
    return 0;
    
}


