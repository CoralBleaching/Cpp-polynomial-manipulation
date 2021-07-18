# Polynomial manipulation in C++

This module implements practical representations and operations on polynomials.

### General description

A polynomial is represented by a list of its coefficients in the usual order, e.g. _x<sup>2</sup> - 5x + 6_ is represented by _{ 1, - 5, 6 }_. More specifically, such list is of the form of an ordered list of doubles, stored in a container of type `vector<double>` inside the class `Polynomial`. A polynomial can be instantiated such as `Polynomial p_x = { 1, -5, 6 };`.

`Polynomial` Defines an object representing a polynomial exactly like a `vector<double>` whose elements can be accessed via brackets (`polynomial[i]`). It contains a single private field named `list` of type `vector<double>` which stores the list of its coefficients in order of decreasing degree. It has many methods and features of a standard vector, such as `push_back`, `pop_back`, an iterator and iterator methods , such as `begin()`, `end()`, `ìnsert()`, `erase()`, and `back()`. It can be evaluated at _x_ via the `()` operator. 

For instance:

	Polynomial p = { 1, 2, 3 };
    cout << p(10);
Output:

	123

### Constructors:

- `Polynomial()` Default constructor.
- `Polynomial(uint n)` Creates a "zero" polynomial of size `n`.
- `Polynomial(std::vector<double>& v)` Initializes the polynomial with vector `v`.
- `Polynomial(const Polynomial& p)` Copy constructor.
- `Polynomial(std::initializer_list<double> l)` Initializes the polynomial object via a list argument.

### Notable methods:

- `to_string()` Creates a string representing the list of coefficients.
- `derivative()` Returns the derivative of the polynomial in the form of another polynomial (of 1 lesser degree).
- `size()` Returns the size of the polynomial (size of the underlying list, *not* its degree).
- `empty()` Returns `true` if the polynomials list of coefficients has no elements, otherwise `false`.
- The aforementioned `begin()`, `end()`, `ìnsert()`, `erase()`, and `back()`.
- `degree(Polynomial)` Returns the degree of the polynomial as an unsigned integer.
- `number of roots(Polynomial p, double a, double b)` Gives the number of real roots located within _[a, b]_.
- `roots(Polynomial)` Gives all the real roots of a polynomial that can be found numerically.

### Overloaded operators:

- `+` Addition is defined for two polynomials or for polynomials with scalars. 
- `-` Subtraction is defined likewise.
- `*` Multiplication is defined likewise.
- `/` Division is defined likewise, however, a scalar cannot be divided by a polynomial. If a polynomial is divided by another polynomial, the quotient is returned.
- `%` Modulo operation, returns the remainder of polynomial division. It's not defined for scalars.
