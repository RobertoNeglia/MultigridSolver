# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt


# exact solution
def u(x, y):
    return np.sin(2.0 * np.pi * x) * np.sin(2.0 * np.pi * y)


# force
def f(x, y):
    return 8.0 * np.pi * np.pi * np.sin(2.0 * np.pi * x) * np.sin(2.0 * np.pi * y)


def plot(X, Y, W, U):
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="3d")
    # Plot a basic wireframe.
    ax.plot_wireframe(X, Y, W, color="r")
    ax.plot_wireframe(X, Y, U, color="b")

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("w")
    plt.show()


def jacobi(A, b, toll, max_iter, x=None):
    # if a initial guess is given use it, otherwise we use all zeros
    x = np.zeros_like(b) if x is None else x
    i = 0
    D = np.diag(A)
    R = A - np.diagflat(D)
    err = np.linalg.norm(A @ x - b)
    while err > toll and max_iter > i:
        i += 1
        x = (b - R @ x) / D
        err = np.linalg.norm(A @ x - b)
    return i, x


# one iteration of the two levels method
def mg_2lvls(A, IhH, b, x, pre, sub, post, toll):
    # coarse level matrix
    AH = IhH @ A @ IhH.T
    # pre-smoothing
    j1, x = jacobi(A, b, toll, pre, x)
    # residual
    r = b - A @ x
    # error correction
    j2, e = jacobi(AH, IhH @ r, toll, sub)
    x += IhH.T @ e
    # post smoothing
    j3, x = jacobi(A, b, toll, post, x)
    # return number of jacobi calls and solution
    return j1 + j2 + j3, x


# multigrid solver
def mg(A, IhH, b, pre, sub, post, toll, max_iter):
    x = np.zeros_like(b)
    i = 0  # number of two-levels calls
    j = 0  # number of jacobi calls

    err = np.linalg.norm(A @ x - b)
    while err > toll and max_iter > i:
        i += 1
        ji, x = mg_2lvls(A, IhH, b, x, pre, sub, post, toll)
        j += ji
        err = np.linalg.norm(A @ x - b)

    print("MG | two lvls iter:", i, "| total jacobi calls:", j)
    return x


def build_matrix(N):
    N2 = (N - 1) * (N - 1)
    A = np.zeros((N2, N2))
    # Diagonal
    for i in range(0, N - 1):
        for j in range(0, N - 1):
            A[i + (N - 1) * j, i + (N - 1) * j] = -4

    # LOWER DIAGONAL
    for i in range(1, N - 1):
        for j in range(0, N - 1):
            A[i + (N - 1) * j, i + (N - 1) * j - 1] = 1
    # UPPPER DIAGONAL
    for i in range(0, N - 2):
        for j in range(0, N - 1):
            A[i + (N - 1) * j, i + (N - 1) * j + 1] = 1

    # LOWER IDENTITY MATRIX
    for i in range(0, N - 1):
        for j in range(1, N - 1):
            A[i + (N - 1) * j, i + (N - 1) * (j - 1)] = 1
    # UPPER IDENTITY MATRIX
    for i in range(0, N - 1):
        for j in range(0, N - 2):
            A[i + (N - 1) * j, i + (N - 1) * (j + 1)] = 1
    return A


def solve(N, u, f):
    h = 1 / N
    x = np.linspace(0, 1.0, N + 1)
    y = np.linspace(0, 1.0, N + 1)
    X, Y = np.meshgrid(x, y)
    Xi, Yi = np.meshgrid(x[1:-1], y[1:-1])

    # w is the full solution
    w = np.zeros((N + 1, N + 1))

    # fix the boundary conditions
    w[0:N, 0] = u(x[0:N], 0.0)  # left Boundary
    w[0:N, N] = u(x[0:N], 1.0)  # Right Boundary
    w[0, 0:N] = u(0.0, y[0:N])  # Lower Boundary
    w[N, 0:N] = u(1.0, y[0:N])  # Upper Boundary

    print("Building matrix")
    A = build_matrix(N)

    print("Building interpolator")
    # Interpolator (assume N even)
    NI, MI = (N - 2) // 2, N - 1
    IhH = np.zeros((NI * NI, MI * MI))

    for il in range(NI):
        i = 1 + il * 2
        for jl in range(NI):
            ijl = il + jl * NI
            j = 1 + jl * 2

            IhH[ijl, i + j * MI] = 1.0 / 4.0
            IhH[ijl, i - 1 + j * MI] = 1.0 / 8.0
            IhH[ijl, i + 1 + j * MI] = 1.0 / 8.0
            IhH[ijl, i + (j - 1) * MI] = 1.0 / 8.0
            IhH[ijl, i + (j + 1) * MI] = 1.0 / 8.0
            IhH[ijl, i + 1 + (j - 1) * MI] = 1.0 / 16.0
            IhH[ijl, i + 1 + (j + 1) * MI] = 1.0 / 16.0
            IhH[ijl, i - 1 + (j - 1) * MI] = 1.0 / 16.0
            IhH[ijl, i - 1 + (j + 1) * MI] = 1.0 / 16.0

    print("Building RHS")
    # RHS
    b = -f(Xi.flatten(), Yi.flatten()) * h * h
    # Add to RHS the boundary conditions contribute
    b[0 : N - 1] -= u(x[1:N], 0.0)  # Bottom Boundary
    b[(N - 1) * (N - 2) : (N - 1) * (N - 2) + N - 1] -= u(x[1:N], 1.0)  # Top Boundary

    b[0 : (N - 1) * (N - 1) : (N - 1)] -= u(0.0, y[1:N])  # Left Boundary
    b[N - 2 : N - 2 + (N - 1) * (N - 1) : (N - 1)] -= u(1.0, y[1:N])  # Right Boundary

    print("Solving linear system")
    toll = 1e-8
    x0 = mg(A, IhH, b, 32, 4, 32, toll, 1e4)
    # compare with jacobi
    j, x0 = jacobi(A, b, toll, 1e4)
    print("plain jacobi iter needed:", j)
    print("res:", np.linalg.norm(A @ x0 - b))
    # put the internal node values in the full solution
    w[1:N, 1:N] = x0.reshape((N - 1, N - 1))
    return X, Y, w


def convergenge(Ns):
    e = []
    for N in Ns:
        X, Y, w = solve(N, u, f)
        U = u(X, Y)
        e.append(np.abs(w - U).max())
        print("err:", e[-1])
        print("=================")
    plt.plot(Ns, e, "o-")
    # check that we converge with order 2 as by theory
    plt.plot(Ns, 5 / np.array(Ns) ** 2, "k--")
    plt.yscale("log")
    plt.xscale("log")


if __name__ == "__main__":
    convergenge([4, 8, 16, 32])
