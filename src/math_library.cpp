#include "math_library.hpp"


using namespace Math_lib;

// Matrix class

// Constructor: Initialize with given value
Matrix::Matrix(size_t r, size_t c, double init) : rows(r), cols(c), data(r * c, init){}

// Constructor: Initialize with list of values
Matrix::Matrix(size_t r, size_t c, std::initializer_list<double> values) 
    : rows(r), cols(c), data(values){
    if (values.size() != r * c){
        throw std::invalid_argument("Initializer list size does not match matrix dimensions.");
    }
}

Matrix::Matrix(size_t r, size_t c, const std::vector<double>& values)
    : rows(r), cols(c), data(values){
    if (values.size() != r * c){
        throw std::invalid_argument("Vector size does not match matrix dimensions.");
    }
}

// Copy constructor
Matrix::Matrix(const Matrix& other) : rows(other.rows), cols(other.cols), data(other.data){}

// Move constructor
Matrix::Matrix(Matrix&& other) noexcept : rows(other.rows), cols(other.cols), data(std::move(other.data)){
    other.rows = 0;
    other.cols = 0;
}

// Assignment operator
Matrix& Matrix::operator=(const Matrix& other){
    if (this != &other){
        rows = other.rows;
        cols = other.cols;
        data = other.data;
    }
    return *this;
}

// Move assignment operator
Matrix& Matrix::operator=(Matrix&& other) noexcept{
    if (this != &other){
        rows = other.rows;
        cols = other.cols;
        data = std::move(other.data);
        other.rows = 0;
        other.cols = 0;
    }
    return *this;
}

// Access elements
double& Matrix::operator()(size_t i, size_t j){
    if (i >= rows || j >= cols) throw std::out_of_range("Index out of bounds");
    return data[i * cols + j];
}

double Matrix::operator()(size_t i, size_t j) const{
    if (i >= rows || j >= cols) throw std::out_of_range("Index out of bounds");
    return data[i * cols + j];
}

// Print matrix
void Matrix::print() const{
    for (size_t i = 0; i < rows; ++i){
        for (size_t j = 0; j < cols; ++j){
            std::cout << std::setw(8) << (*this)(i, j) << " ";
        }
        std::cout << "\n";
    }
}

// Transpose matrix
Matrix Matrix::transpose() const{
    Matrix result(cols, rows);
    for (size_t i = 0; i < rows; ++i){
        for (size_t j = 0; j < cols; ++j){
            result(j, i) = (*this)(i, j);
        }
    }
    return result;
}

// Matrix addition
Matrix Matrix::operator+(const Matrix& other) const{
    if (rows != other.rows || cols != other.cols) throw std::invalid_argument("Matrix dimensions do not match");
    Matrix result(rows, cols);
    for (size_t i = 0; i < rows * cols; ++i){
        result.data[i] = data[i] + other.data[i];
    }
    return result;
}

// Matrix multiplication
Matrix Matrix::operator*(const Matrix& other) const{
    if (cols != other.rows) throw std::invalid_argument("Matrix multiplication not possible");
    Matrix result(rows, other.cols);
    for (size_t i = 0; i < rows; ++i){
        for (size_t j = 0; j < other.cols; ++j){
            double sum = 0.0;
            for (size_t k = 0; k < cols; ++k){
                sum += (*this)(i, k) * other(k, j);
            }
            result(i, j) = sum;
        }
    }
    return result;
}

// Scalar multiplication
Matrix Matrix::scalarMultiply(double scalar) const{
    Matrix result(rows, cols);
    for (size_t i = 0; i < data.size(); ++i){
        result.data[i] = data[i] * scalar;
    }
    return result;
}

// Determinant
double Matrix::determinant() const{
    if (rows != cols) throw std::invalid_argument("Determinant is only defined for square matrices");
    if (rows == 1) return data[0];
    if (rows == 2) return data[0] * data[3] - data[1] * data[2];
    double det = 0.0;
    for (size_t i = 0; i < cols; ++i){
        Matrix submatrix(rows - 1, cols - 1);
        for (size_t j = 1; j < rows; ++j){
            for (size_t k = 0; k < cols; ++k){
                if (k < i){
                    submatrix(j - 1, k) = (*this)(j, k);
                } else if (k > i){
                    submatrix(j - 1, k - 1) = (*this)(j, k);
                }
            }
        }
        det += (i % 2 == 0 ? 1 : -1) * (*this)(0, i) * submatrix.determinant();
    }
    return det;
}

// Inverse
Matrix Matrix::inverse() const{
    if (rows != cols) throw std::invalid_argument("Inverse is only defined for square matrices");
    double det = determinant();
    if (det == 0) throw std::invalid_argument("Matrix is singular");
    Matrix result(rows, cols);
    for (size_t i = 0; i < rows; ++i){
        for (size_t j = 0; j < cols; ++j){
            Matrix submatrix(rows - 1, cols - 1);
            for (size_t k = 0; k < rows; ++k){
                if (k == i) continue;
                for (size_t l = 0; l < cols; ++l){
                    if (l == j) continue;
                    submatrix(k < i ? k : k - 1, l < j ? l : l - 1) = (*this)(k, l);
                }
            }
            result(j, i) = (i + j) % 2 == 0 ? 1 : -1;
            result(j, i) *= submatrix.determinant();
            result(j, i) /= det;
        }
    }
    return result;
}

// Identity matrix
Matrix Matrix::identity(size_t size){
    Matrix result(size, size);
    for (size_t i = 0; i < size; ++i){
        result(i, i) = 1.0;
    }
    return result;
}

// Random matrix generator
Matrix Matrix::random(size_t r, size_t c, double min_val, double max_val){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(min_val, max_val);
    
    Matrix result(r, c);
    for (size_t i = 0; i < r * c; ++i){
        result.data[i] = dist(gen);
    }
    return result;
}

namespace Math_lib {
    // Statistics functions
    long double mean(std::vector<long double> data){
        long double sum = 0;
        for (double i : data){
            sum += i;
        }
        return sum / data.size();
    }

    long double variance(std::vector<long double> data){
        long double m = mean(data);
        long double sum = 0;
        for (double i : data){
            sum += (i - m) * (i - m);
        }
        return sum / data.size();
    }

    long double standardDeviation(std::vector<long double> data){
        return sqrt(variance(data));
    }

    long double covariance(std::vector<long double> data1, std::vector<long double> data2){
        if (data1.size() != data2.size()){
            throw std::invalid_argument("Data sizes do not match");
        }
        long double m1 = mean(data1);
        long double m2 = mean(data2);
        long double sum = 0;
        for (size_t i = 0; i < data1.size(); ++i){
            sum += (data1[i] - m1) * (data2[i] - m2);
        }
        return sum / data1.size();
    }

    long double correlation(std::vector<long double> data1, std::vector<long double> data2){
        return covariance(data1, data2) / (standardDeviation(data1) * standardDeviation(data2));
    }

    long double median(std::vector<long double> data){
        std::sort(data.begin(), data.end());
        size_t n = data.size();
        if (n % 2 == 0){
            return (data[n / 2 - 1] + data[n / 2]) / 2;
        } else {
            return data[n / 2];
        }
    }

    long double probabilityDensityFunction(long double x, long double mu, long double sigma){
        return exp(-0.5 * (x - mu) * (x - mu) / (sigma * sigma)) / (sigma * sqrt(2 * M_PI));
    }

    long double cumulativeDistributionFunction(long double x, long double mu, long double sigma){
        return 0.5 * (1 + erf((x - mu) / (sigma * sqrt(2))));
    }

    long double inverseCumulativeDistributionFunction(long double p, long double mu, long double sigma){
        return mu + sigma * sqrt(2) * boost::math::erf_inv(2 * p - 1);
    }
    
    long double normalDistribution(long double x, long double mu, long double sigma){
        return exp(-0.5 * (x - mu) * (x - mu) / (sigma * sigma)) / (sigma * sqrt(2 * M_PI));
    }

    long double inverseNormalDistribution(long double p, long double mu, long double sigma){
        return mu + sigma * sqrt(2) * boost::math::erf_inv(2 * p - 1);
    }

    long double areaUnderNormalDistribution(long double x1, long double x2, long double mu, long double sigma){
        return cumulativeDistributionFunction(x2, mu, sigma) - cumulativeDistributionFunction(x1, mu, sigma);
    }
}