# Intel Intelligence Projects

## Overview
This project is a comprehensive math library implemented in C++. It includes a variety of mathematical functions, data structures, and statistical operations. The library is designed to be modular and easy to use, making it suitable for a wide range of applications.

## Features
- **Matrix Operations**: Addition, multiplication, scalar multiplication, determinant, inverse, transpose, identity matrix, and random matrix generation.
- **Number Theory Functions**: GCD, LCM, prime checking, Armstrong number checking, perfect number checking, prime factorization, Fibonacci, sum of divisors, number of divisors, and Euler's totient function.
- **Data Structures**: Stack, Queue, LinkedList, and BinaryTree with various operations.
- **Statistical Functions**: Mean, variance, standard deviation, covariance, correlation, median, probability density function, cumulative distribution function, inverse cumulative distribution function, normal distribution, and area under normal distribution.

## Installation
To use this library, you need to have a C++ compiler and the Boost library installed on your system.

### Install Boost Library
On Ubuntu, you can install the Boost library using the following command:
```sh
sudo apt-get install libboost-all-dev
```

## Usage
Include the `math_library.hpp` header file in your project to use the library functions.

### Example
Here is an example of how to use the library in your project:

#### User's Source Code (`user_code.cpp`)
```cpp
#include <iostream>
#include "math_library.hpp"

using namespace Math_lib;

int main() {
    // Example usage of the library
    Matrix m1(2, 2, 1.0);
    Matrix m2(2, 2, {1.0, 2.0, 3.0, 4.0});
    Matrix m3 = m1 + m2;
    m3.print();

    int a = 12, b = 18;
    std::cout << "GCD of " << a << " and " << b << ": " << gcd(a, b) << "\n";

    Stack<int> stack;
    stack.push(1);
    stack.push(2);
    std::cout << "Stack top: " << stack.top() << "\n";

    std::vector<long double> data = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::cout << "Mean: " << mean(data) << "\n";

    return 0;
}
```

### Compilation
To compile the library and the main program, use the following command:

```sh
g++ -o user_program user_code.cpp src/math_library.cpp -lboost_system -lboost_filesystem
```

### Running the Program
After compiling, you can run the resulting executable:

```sh
./user_program
```

## Project Structure
```
Intel_Inteligence_Projects/
├── LICENSE
├── README.md
└── src/
    ├── main.cpp
    ├── math_library.cpp
    └── math_library.hpp
```

## License
This project is licensed under the MIT License.