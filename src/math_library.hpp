#ifndef MATH_LIBRARY_HPP
#define MATH_LIBRARY_HPP

#include <iostream>
#include <iomanip>
#include <vector>
#include <random>
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <initializer_list>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/math/special_functions/erf.hpp>

namespace Math_lib{
    // Matrices
    class Matrix{
    private:
        std::vector<double> data;
        size_t rows, cols;

    public:
        // Constructors
        Matrix(size_t r, size_t c, double init = 0.0);
        Matrix(size_t r, size_t c, std::initializer_list<double> values);
        Matrix(size_t r, size_t c, const std::vector<double>& values);
        Matrix(const Matrix& other); // Copy constructor
        Matrix(Matrix&& other) noexcept; // Move constructor

        // Operators
        Matrix& operator=(const Matrix& other);
        Matrix& operator=(Matrix&& other) noexcept;
        double& operator()(size_t i, size_t j);
        double operator()(size_t i, size_t j) const;

        // Utility functions
        size_t rowCount() const{ return rows; }
        size_t colCount() const{ return cols; }
        void print() const;

        // Matrix operations
        Matrix transpose() const;
        Matrix operator+(const Matrix& other) const;
        Matrix operator*(const Matrix& other) const;
        Matrix scalarMultiply(double scalar) const;
        double determinant() const;
        Matrix inverse() const;

        // Special Matrices
        static Matrix identity(size_t size);
        static Matrix random(size_t r, size_t c, double min_val = 0.0, double max_val = 1.0);
    };

    // Number theory functions
    template <typename T>
    typename std::enable_if<std::is_integral<T>::value || std::is_same<T, boost::multiprecision::cpp_int>::value, T>::type
    powerInt(T p, T a){
        T result = 1;
        while (a > 0){
            if (a % 2 == 1){
                result *= p;
            }
            p *= p;
            a /= 2;
        }
        return result;
    }

    template <typename T>
    typename std::enable_if<std::is_integral<T>::value || std::is_same<T, boost::multiprecision::cpp_int>::value, T>::type
    gcd(T a, T b){
        while (b != 0){
            T temp = b;
            b = a % b;
            a = temp;
        }
        return a;
    }

    template <typename T>
    typename std::enable_if<std::is_integral<T>::value || std::is_same<T, boost::multiprecision::cpp_int>::value, T>::type
    lcm(T a, T b){
        return (a / gcd(a, b)) * b;
    }

    template <typename T>
    typename std::enable_if<std::is_integral<T>::value || std::is_same<T, boost::multiprecision::cpp_int>::value, bool>::type
    isPrime(T n){
        if (n < 2) return false;
        if (n == 2) return true;
        if (n % 2 == 0) return false;
        for (T i = 3; i * i <= n; i += 2){
            if (n % i == 0) return false;
        }
        return true;
    }

    template <typename T>
    typename std::enable_if<std::is_integral<T>::value || std::is_same<T, boost::multiprecision::cpp_int>::value, bool>::type
    isArmstrongNumber(T n){
        T sum = 0;
        T temp = n;
        while (temp > 0){
            T digit = temp % 10;
            sum += powerInt(digit, static_cast<T>(std::to_string(n).size()));
            temp /= 10;
        }
        return sum == n;
    }

    template <typename T>
    typename std::enable_if<std::is_integral<T>::value || std::is_same<T, boost::multiprecision::cpp_int>::value, bool>::type
    isPerfectNumber(T n){
        T sum = 1;
        for (T i = 2; i * i <= n; i++){
            if (n % i == 0){
                sum += i;
                if (i * i != n){
                    sum += n / i;
                }
            }
        }
        return sum == n;
    }

    template <typename T>
    typename std::enable_if<std::is_integral<T>::value || std::is_same<T, boost::multiprecision::cpp_int>::value, std::vector<std::vector<T>>>::type
    primeFactor(T n){
        std::vector<std::vector<T>> factors;
        if (n % 2 == 0){
            factors.push_back({2, 0});
            while (n % 2 == 0){
                n /= 2;
                factors.back()[1]++;
            }
        }
        for (T i = 3; i * i <= n; i += 2){
            if (n % i == 0){
                factors.push_back({i, 0});
                while (n % i == 0){
                    n /= i;
                    factors.back()[1]++;
                }
            }
        }
        if (n > 1){
            factors.push_back({n, 1});
        }
        return factors;
    }

    template <typename T>
    typename std::enable_if<std::is_integral<T>::value || std::is_same<T, boost::multiprecision::cpp_int>::value, T>::type
    fastFibonacci(T n){
        if (n == 0) return 0;
        if (n == 1) return 1;

        T a = fibonacci(n / 2);
        T b = fibonacci(n / 2 + 1);

        if (n % 2 == 0){
            return a * (2 * b - a);
        } else{
            return b * b + a * a;
        }
    }

    template <typename T>
    typename std::enable_if<std::is_integral<T>::value || std::is_same<T, boost::multiprecision::cpp_int>::value, T>::type
    sumDivisors(T n){
        if (n == 1) return 1;

        std::vector<std::vector<T>> factors = primeFactor(n);

        T sum = 1;

        for (const auto& factor : factors){
            sum *= ((powerInt(factor[0], factor[1] + 1) - 1) / (factor[0] - 1));
        }

        return sum;
    }

    template <typename T>
    typename std::enable_if<std::is_integral<T>::value || std::is_same<T, boost::multiprecision::cpp_int>::value, T>::type
    numberOfDivisors(T n){
        if (n == 1) return 1;

        std::vector<std::vector<T>> factors = primeFactor(n);

        T divisors = 1;
        for (const auto& factor : factors){
            divisors *= (factor[1] + 1);
        }

        return divisors;
    }

    template <typename T>
    typename std::enable_if<std::is_integral<T>::value || std::is_same<T, boost::multiprecision::cpp_int>::value, T>::type
    totient(T n){
        if (n == 1) return 1;

        std::vector<std::vector<T>> factors = primeFactor(n);

        T result = n;

        for (const auto& factor : factors){
            result = result / factor[0] * (factor[0] - 1);
        }

        return result;
    }

    // Data structures

    template <typename T>
    class Stack{
    private:
        std::vector<T> data;
    public:
        Stack() = default;
        Stack(const Stack& other) = default;
        Stack(Stack&& other) noexcept = default;
        Stack& operator=(const Stack& other) = default;
        Stack& operator=(Stack&& other) noexcept = default;

        void push(const T& value){ data.push_back(value); }
        void pop(){ data.pop_back(); }
        T& top(){ return data.back(); }
        const T& top() const{ return data.back(); }
        bool empty() const{ return data.empty(); }
        size_t size() const{ return data.size(); }
    };

    template <typename T>
    class Queue{
    private:
        std::vector<T> data;
    public:
        Queue() = default;
        Queue(const Queue& other) = default;
        Queue(Queue&& other) noexcept = default;
        Queue& operator=(const Queue& other) = default;
        Queue& operator=(Queue&& other) noexcept = default;

        void push(const T& value){ data.push_back(value); }
        void pop(){ data.erase(data.begin()); }
        T& front(){ return data.front(); }
        const T& front() const{ return data.front(); }
        T& back(){ return data.back(); }
        const T& back() const{ return data.back(); }
        bool empty() const{ return data.empty(); }
        size_t size() const{ return data.size(); }
    };

    template <typename T>
    class LinkedList{
    private:
        struct Node{
            T data;
            Node* next;
            Node* prev;
        };

        Node* head;
        Node* tail;
        size_t length;
    public:
        LinkedList() : head(nullptr), tail(nullptr), length(0){}
        LinkedList(const LinkedList& other) = delete;
        LinkedList(LinkedList&& other) noexcept = delete;
        LinkedList& operator=(const LinkedList& other) = delete;
        LinkedList& operator=(LinkedList&& other) noexcept = delete;
        ~LinkedList(){
            while (head != nullptr){
                Node* temp = head;
                head = head->next;
                delete temp;
            }
        }

        void pushFront(const T& value){
            Node* newNode = new Node{value, head, nullptr};
            if (head != nullptr){
                head->prev = newNode;
            } else{
                tail = newNode;
            }
            head = newNode;
            length++;
        }

        void pushBack(const T& value){
            Node* newNode = new Node{value, nullptr, tail};
            if (tail != nullptr){
                tail->next = newNode;
            } else{
                head = newNode;
            }
            tail = newNode;
            length++;
        }

        void popFront(){
            if (head == nullptr) return;
            Node* temp = head;
            head = head->next;
            if (head != nullptr){
                head->prev = nullptr;
            } else{
                tail = nullptr;
            }
            delete temp;
            length--;
        }

        void popBack(){
            if (tail == nullptr) return;
            Node* temp = tail;
            tail = tail->prev;
            if (tail != nullptr){
                tail->next = nullptr;
            } else{
                head = nullptr;
            }
            delete temp;
            length--;
        }

        T& front(){ return head->data; }
        const T& front() const{ return head->data; }
        T& back(){ return tail->data; }
        const T& back() const{ return tail->data; }
        bool empty() const{ return head == nullptr; }
        size_t size() const{ return length; }
    };

    template <typename T>
    class BinaryTree{
    private:
        struct Node{
            T data;
            Node* left;
            Node* right;
        };

        Node* root;

        Node* insert(Node* node, const T& value){
            if (node == nullptr){
                return new Node{value, nullptr, nullptr};
            }
            if (value < node->data){
                node->left = insert(node->left, value);
            } else{
                node->right = insert(node->right, value);
            }
            return node;
        }

        Node* findMin(Node* node){
            while (node->left != nullptr){
                node = node->left;
            }
            return node;
        }

        Node* remove(Node* node, const T& value){
            if (node == nullptr) return nullptr;
            if (value < node->data){
                node->left = remove(node->left, value);
            } else if (value > node->data){
                node->right = remove(node->right, value);
            } else{
                if (node->left == nullptr){
                    Node* temp = node->right;
                    delete node;
                    return temp;
                } else if (node->right == nullptr){
                    Node* temp = node->left;
                    delete node;
                    return temp;
                }
                Node* temp = findMin(node->right);
                node->data = temp->data;
                node->right = remove(node->right, temp->data);
            }
            return node;
        }

        void inorder(Node* node, std::vector<T>& result){
            if (node == nullptr) return;
            inorder(node->left, result);
            result.push_back(node->data);
            inorder(node->right, result);
        }

        void preorder(Node* node, std::vector<T>& result){
            if (node == nullptr) return;
            result.push_back(node->data);
            preorder(node->left, result);
            preorder(node->right, result);
        }

        void postorder(Node* node, std::vector<T>& result){
            if (node == nullptr) return;
            postorder(node->left, result);
            postorder(node->right, result);
            result.push_back(node->data);
        }

        void levelOrder(Node* node, std::vector<T>& result){
            if (node == nullptr) return;
            std::queue<Node*> q;
            q.push(node);
            while (!q.empty()){
                Node* current = q.front();
                q.pop();
                result.push_back(current->data);
                if (current->left != nullptr) q.push(current->left);
                if (current->right != nullptr) q.push(current->right);
            }
        }

        int height(Node* node){
            if (node == nullptr) return 0;
            return 1 + std::max(height(node->left), height(node->right));
        }

        bool search(Node* node, const T& value){
            if (node == nullptr) return false;
            if (node->data == value) return true;
            if (value < node->data){
                return search(node->left, value);
            } else{
                return search(node->right, value);
            }
        }

        void clear(Node* node){
            if (node == nullptr) return;
            clear(node->left);
            clear(node->right);
            delete node;
        }

    public:
        BinaryTree() : root(nullptr){}
        BinaryTree(const BinaryTree& other) = delete;
        BinaryTree(BinaryTree&& other) noexcept = delete;
        BinaryTree& operator=(const BinaryTree& other) = delete;
        BinaryTree& operator=(BinaryTree&& other) noexcept = delete;
        ~BinaryTree(){
            clear(root);
        }

        void insert(const T& value){
            root = insert(root, value);
        }

        void remove(const T& value){
            root = remove(root, value);
        }

        std::vector<T> inorder(){
            std::vector<T> result;
            inorder(root, result);
            return result;
        }

        std::vector<T> preorder(){
            std::vector<T> result;
            preorder(root, result);
            return result;
        }

        std::vector<T> postorder(){
            std::vector<T> result;
            postorder(root, result);
            return result;
        }

        std::vector<T> levelOrder(){
            std::vector<T> result;
            levelOrder(root, result);
            return result;
        }

        int height(){
            return height(root);
        }

        bool search(const T& value){
            return search(root, value);
        }

        bool empty() const{
            return root == nullptr;
        }

        void clear(){
            clear(root);
            root = nullptr;
        }
    };

    // Statistics functions
    long double mean(std::vector<long double> data);
    long double variance(std::vector<long double> data);
    long double standardDeviation(std::vector<long double> data);
    long double covariance(std::vector<long double> data1, std::vector<long double> data2);
    long double correlation(std::vector<long double> data1, std::vector<long double> data2);
    long double median(std::vector<long double> data);
    long double probabilityDensityFunction(long double x, long double mean, long double variance);
    long double cumulativeDistributionFunction(long double x, long double mean, long double variance);
    long double inverseCumulativeDistributionFunction(long double p, long double mean, long double variance);
    long double normalDistribution(long double x, long double mean, long double variance);
    long double inverseNormalDistribution(long double p, long double mean, long double variance);
    long double areaUnderNormalDistribution(long double x1, long double x2, long double mean, long double variance);
}   // namespace Math_lib

#endif // MATH_LIBRARY_HPP