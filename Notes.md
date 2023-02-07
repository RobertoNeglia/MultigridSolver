## Some Notes ##
- Overall good progamming with templates and polymorphism. You have created your own matrix and vector classes. OK, since it is an exercise for learning to program, but in general better use well established linear algebra libraries
- Code compiles and runs smoothly, you provide interesting examples in the readme
- Some comments in the code help readability.... here there are a few, but very little to describe the function and function arguments.

- You have used templates in other parts of the code. Why not make the class for multigrid a template with template argument the smoother? So using different smoothers
would just mean change the template parameter. Or I am missing something?

- Some constructors are defined but do nothing.  Example: `SparseMatrix() {}`. If the synthetic operator is fine you do not need to define your own. Moreover, try to initialize POD variables (you can do it in-class), don't assume they are initialised by 0. In `SparseMatrix` better replace

```
  int          nnz;
  unsigned int n_rows;
  unsigned int n_cols;

```

with
```
  int          nnz=0;
  unsigned int n_rows=0;
  unsigned int n_cols=0;

```

- In `SparseMatrix` you define a custom copy-constructor and assignment operator, but they achieve nothing different from the default ones.

-Possible memory leaks in many points!

```cpp
static void *
  operator new(std::size_t count) {
    return ::operator new(count);
  }
//...
Vector<double> *Av = new Vector<double>();
//...
SparseMatrix *AB = new SparseMatrix;
//...
SparseMatrix *At = new SparseMatrix;
// in jacobi_r.hpp
SparseMatrix *D = new SparseMatrix;
// in multigrid_solver.hpp
SparseMatrix *A_2h_ = new SparseMatrix;
```

where is the corresponding `delete`? Who is in charge of the resource? In modern C++ you should use smart pointers for owning pointers, and avoid possible terrible 
headaches caused by a memory leak!. Remember RAII and the single responsibility principle! An object should be responsible of properly managing its own resources, and have possible a single responsibility. Give the responsibility of memory handling to smart pointers and you live happy. This applies also elsewhere in your code.

- No need to specify as const the value returned by a function. It is useless. `std::pair<double, bool> const coeff()` should e just  `std::pair<double, bool> coeff()` (it's only a warning though,the compiler discards the const). 

- Why `static` in 

```
    template <typename T>
    static double
    norm(const Vector<T> &a) {
```
A static function is visible only in the translation unit that defines it. But this is a function template in a header file. I do not understand the point of making it static. 
Moreover, you are performing a reduction operator, so either you make the update atomic, or you specify `reduction(+:n)` in the pragma.



## Small stuff ##
- `coeff_ref` returns a pair with a pointer and a boolean. The reason for the boolean I assume is to signal whether the element is not present.
It is fine, but it makes the use of this function a bit cumbersome. Indeed, the fact that the element is not present is already indicated by the null pointer. You do not need the bool.

- When you have a fatal error better throw an exception than just print 
on the standard output. It is not so complicate if you use standard exceptions:
```
    if (n_cols != B.rows()) {
        std::cout << "ERROR: INCOMPATIBLE SIZES" << std::endl;
        return *AB;
```

becomes
```
    if (n_cols != B.rows())
      throw std::runtime_error("ERROR: INCOMPATIBLE SIZES");
```
That's it.

- I do not know if you have written it or some IDE you are using, but calling preprocessor macros as `__SPARSE_MATRIX_R_H__` is dangerous. The pattern `__Name__` and also `__Name` or `_Name` is used internally for system variables, so there is the risk of a name clash. A good rule is to avoid `_` at the beginning of a name altogether. `SPARSE_MATRIX_R_H__` is fine.

- No need to use `return` in a function declared `void` when you are already at the end of the function 

```
  void
  mul(std::unique_ptr<Vector<double>> &res, const Vector<double> &v) const
  {
  ...
  return; // useless!
  }
  
```

- In C++ `time.h` should be replaced by `<ctime>`.

- In `jacobi_r.hpp` you could have stored the preconditioner in a `Vector` since it is diagonal, and when you apply it you just make an element-wise multiplication between vectors.

- In `unsigned int M = std::floor((N - 1) / 2);` the function `std::floor` may just cause errors: the operator `/` between two integers is already the integer division, when you call `std::floor` it outputs a floating point number that is then implicitly cast to unsigned integer.
