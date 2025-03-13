#include <iostream>
#include "math_library.hpp"

using namespace Math_lib;

void testMatrix(){
    Matrix m1(2, 2, 1.0);
    Matrix m2(2, 2, {1.0, 2.0, 3.0, 4.0});
    Matrix m3 = m1 + m2;
    Matrix m4 = m2 * m2;
    Matrix m5 = m2.scalarMultiply(2.0);
    double det = m2.determinant();
    Matrix m6 = m2.inverse();
    Matrix m7 = Matrix::identity(3);
    Matrix m8 = Matrix::random(3, 3, 0.0, 1.0);

    std::cout << "Matrix m1:\n";
    m1.print();
    std::cout << "Matrix m2:\n";
    m2.print();
    std::cout << "Matrix m3 (m1 + m2):\n";
    m3.print();
    std::cout << "Matrix m4 (m2 * m2):\n";
    m4.print();
    std::cout << "Matrix m5 (m2 * 2.0):\n";
    m5.print();
    std::cout << "Determinant of m2: " << det << "\n";
    std::cout << "Inverse of m2:\n";
    m6.print();
    std::cout << "Identity matrix of size 3:\n";
    m7.print();
    std::cout << "Random matrix 3x3:\n";
    m8.print();
}

void testNumberTheory(){
    int a = 12, b = 18;
    std::cout << "GCD of " << a << " and " << b << ": " << gcd(a, b) << "\n";
    std::cout << "LCM of " << a << " and " << b << ": " << lcm(a, b) << "\n";
    std::cout << "Is 17 prime? " << (isPrime(17) ? "Yes" : "No") << "\n";
    std::cout << "Is 153 an Armstrong number? " << (isArmstrongNumber(153) ? "Yes" : "No") << "\n";
    std::cout << "Is 28 a perfect number? " << (isPerfectNumber(28) ? "Yes" : "No") << "\n";
    std::cout << "Prime factors of 28: ";
    auto factors = primeFactor(28);
    for (const auto& factor : factors){
        std::cout << "(" << factor[0] << "^" << factor[1] << ") ";
    }
    std::cout << "\n";
    boost::multiprecision::cpp_int n = 116;
    std::cout << "Fibonacci of 116: " << fastFibonacci(n) << "\n";
    std::cout << "Sum of divisors of 28: " << sumDivisors(28) << "\n";
    std::cout << "Number of divisors of 28: " << numberOfDivisors(28) << "\n";
    std::cout << "Euler's totient function of 28: " << totient(28) << "\n";
}

void testDataStructures(){
    Stack<int> stack;
    stack.push(1);
    stack.push(2);
    stack.push(3);
    std::cout << "Stack top: " << stack.top() << "\n";
    stack.pop();
    std::cout << "Stack top after pop: " << stack.top() << "\n";

    Queue<int> queue;
    queue.push(1);
    queue.push(2);
    queue.push(3);
    std::cout << "Queue front: " << queue.front() << "\n";
    queue.pop();
    std::cout << "Queue front after pop: " << queue.front() << "\n";

    LinkedList<int> list;
    list.pushBack(1);
    list.pushBack(2);
    list.pushBack(3);
    std::cout << "LinkedList front: " << list.front() << "\n";
    list.popFront();
    std::cout << "LinkedList front after popFront: " << list.front() << "\n";

    BinaryTree<int> tree;
    tree.insert(5);
    tree.insert(3);
    tree.insert(7);
    tree.insert(2);
    tree.insert(4);
    tree.insert(6);
    tree.insert(8);
    std::cout << "BinaryTree inorder: ";
    for (const auto& val : tree.inorder()){
        std::cout << val << " ";
    }
    std::cout << "\n";
    std::cout << "BinaryTree height: " << tree.height() << "\n";
    std::cout << "BinaryTree search 4: " << (tree.search(4) ? "Found" : "Not Found") << "\n";
    tree.remove(4);
    std::cout << "BinaryTree inorder after removing 4: ";
    for (const auto& val : tree.inorder()){
        std::cout << val << " ";
    }
    std::cout << "\n";
}

void testStatistics(){
    std::vector<long double> data = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::cout << "Mean: " << mean(data) << "\n";
    std::cout << "Variance: " << variance(data) << "\n";
    std::cout << "Standard Deviation: " << standardDeviation(data) << "\n";
    std::vector<long double> data2 = {2.0, 3.0, 4.0, 5.0, 6.0};
    std::cout << "Covariance: " << covariance(data, data2) << "\n";
    std::cout << "Correlation: " << correlation(data, data2) << "\n";
    std::cout << "Median: " << median(data) << "\n";
    std::cout << "PDF at 3.0: " << normalProbabilityDensity(3.0, mean(data), variance(data)) << "\n";
    std::cout << "CDF at 3.0: " << cumulativeDistributionFunction(3.0, mean(data), variance(data)) << "\n";
    std::cout << "Inverse CDF at 0.5: " << inverseCumulativeDistributionFunction(0.5, mean(data), variance(data)) << "\n";
    std::cout << "Area under Normal Distribution between 2.0 and 4.0: " << areaUnderNormalDistribution(2.0, 4.0, mean(data), variance(data)) << "\n";
}

int main(){
    std::cout << "Testing Matrix operations:\n";
    testMatrix();
    std::cout << "\nTesting Number Theory functions:\n";
    testNumberTheory();
    std::cout << "\nTesting Data Structures:\n";
    testDataStructures();
    std::cout << "\nTesting Statistics functions:\n";
    testStatistics();
    return 0;
}