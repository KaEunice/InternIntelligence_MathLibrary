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

// Constructor: Initialize with vector of values
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

// Copy assignment operator
Matrix& Matrix::operator=(const Matrix& other){
    if (this != &other){
        rows = other.rows;
        cols = other.cols;
        data = other.data;
    }
    return *this;
}

// Move assignment operator
Matrix& Matrix::operator=(Matrix&& other) noexcept {
    if (this != &other){
        rows = other.rows;
        cols = other.cols;
        data = std::move(other.data);
        other.rows = 0;
        other.cols = 0;
    }
    return *this;
}

// Element access
double& Matrix::operator()(size_t i, size_t j){
    return data[i * cols + j];
}

double Matrix::operator()(size_t i, size_t j) const{
    return data[i * cols + j];
}

// Print matrix
void Matrix::print() const{
    for (size_t i = 0; i < rows; ++i){
        for (size_t j = 0; j < cols; ++j){
            std::cout << std::setw(10) << (*this)(i, j) << " ";
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
    if (rows != other.rows || cols != other.cols){
        throw std::invalid_argument("Matrix dimensions must match for addition.");
    }
    Matrix result(rows, cols);
    for (size_t i = 0; i < rows * cols; ++i){
        result.data[i] = data[i] + other.data[i];
    }
    return result;
}

// Matrix multiplication
Matrix Matrix::operator*(const Matrix& other) const{
    if (cols != other.rows){
        throw std::invalid_argument("Matrix dimensions must match for multiplication.");
    }
    Matrix result(rows, other.cols);
    for (size_t i = 0; i < rows; ++i){
        for (size_t j = 0; j < other.cols; ++j){
            for (size_t k = 0; k < cols; ++k){
                result(i, j) += (*this)(i, k) * other(k, j);
            }
        }
    }
    return result;
}

// Scalar multiplication
Matrix Matrix::scalarMultiply(double scalar) const{
    Matrix result(rows, cols);
    for (size_t i = 0; i < rows * cols; ++i){
        result.data[i] = data[i] * scalar;
    }
    return result;
}

// Determinant
double Matrix::determinant() const{
    if (rows != cols){
        throw std::invalid_argument("Determinant is only defined for square matrices.");
    }
    if (rows == 1){
        return data[0];
    }
    if (rows == 2){
        return data[0] * data[3] - data[1] * data[2];
    }
    double det = 0.0;
    for (size_t i = 0; i < cols; ++i){
        Matrix submatrix(rows - 1, cols - 1);
        for (size_t j = 1; j < rows; ++j){
            size_t sub_col = 0;
            for (size_t k = 0; k < cols; ++k){
                if (k == i) continue;
                submatrix(j - 1, sub_col++) = (*this)(j, k);
            }
        }
        det += (i % 2 == 0 ? 1 : -1) * data[i] * submatrix.determinant();
    }
    return det;
}

// Inverse matrix
Matrix Matrix::inverse() const{
    if (rows != cols){
        throw std::invalid_argument("Inverse is only defined for square matrices.");
    }
    double det = determinant();
    if (det == 0){
        throw std::runtime_error("Matrix is singular and cannot be inverted.");
    }
    Matrix result(rows, cols);
    for (size_t i = 0; i < rows; ++i){
        for (size_t j = 0; j < cols; ++j){
            Matrix submatrix(rows - 1, cols - 1);
            for (size_t k = 0; k < rows; ++k){
                if (k == i) continue;
                size_t sub_row = k < i ? k : k - 1;
                for (size_t l = 0; l < cols; ++l){
                    if (l == j) continue;
                    size_t sub_col = l < j ? l : l - 1;
                    submatrix(sub_row, sub_col) = (*this)(k, l);
                }
            }
            result(j, i) = ((i + j) % 2 == 0 ? 1 : -1) * submatrix.determinant() / det;
        }
    }
    return result;
}

// Special matrices
Matrix Matrix::identity(size_t size){
    Matrix result(size, size);
    for (size_t i = 0; i < size; ++i){
        result(i, i) = 1.0;
    }
    return result;
}

Matrix Matrix::random(size_t r, size_t c, double min_val, double max_val){
    Matrix result(r, c);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(min_val, max_val);
    for (size_t i = 0; i < r * c; ++i){
        result.data[i] = dis(gen);
    }
    return result;
}

// Statistics functions
namespace Math_lib {
    long double mean(std::vector<long double> data){
        return std::accumulate(data.begin(), data.end(), 0.0L) / data.size();
    }

    long double variance(std::vector<long double> data){
        long double m = mean(data);
        long double sum = 0.0L;
        for (auto& val : data){
            sum += (val - m) * (val - m);
        }
        return sum / data.size();
    }

    long double standardDeviation(std::vector<long double> data){
        return std::sqrt(variance(data));
    }

    long double covariance(std::vector<long double> data1, std::vector<long double> data2){
        if (data1.size() != data2.size()){
            throw std::invalid_argument("Data vectors must have the same size.");
        }
        long double mean1 = mean(data1);
        long double mean2 = mean(data2);
        long double sum = 0.0L;
        for (size_t i = 0; i < data1.size(); ++i){
            sum += (data1[i] - mean1) * (data2[i] - mean2);
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
            return (data[n / 2 - 1] + data[n / 2]) / 2.0L;
        } else {
            return data[n / 2];
        }
    }

    long double normalProbabilityDensity(long double x, long double mean, long double variance){
        long double stddev = std::sqrt(variance);
        return (1.0L / (stddev * std::sqrt(2.0L * M_PI))) * std::exp(-0.5L * std::pow((x - mean) / stddev, 2));
    }

    long double cumulativeDistributionFunction(long double x, long double mean, long double variance){
        long double stddev = std::sqrt(variance);
        return 0.5L * (1.0L + boost::math::erf((x - mean) / (stddev * std::sqrt(2.0L))));
    }

    long double inverseCumulativeDistributionFunction(long double p, long double mean, long double variance){
        long double stddev = std::sqrt(variance);
        return mean + stddev * std::sqrt(2.0L) * boost::math::erf_inv(2.0L * p - 1.0L);
    }

    long double areaUnderNormalDistribution(long double x1, long double x2, long double mean, long double variance){
        return cumulativeDistributionFunction(x2, mean, variance) - cumulativeDistributionFunction(x1, mean, variance);
    }
}  // namespace Math_lib