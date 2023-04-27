#include <iostream>
#include <vector>
#include <cmath>




using namespace std;

template<typename T>
class Matrix {
protected:
    int n, m;
    vector<vector<T>> matrix;
public:
    Matrix(int n, int m) {
        this->n = n;
        this->m = m;
        matrix.resize(n);
        for (int i = 0; i < n; i++) {
            matrix[i].resize(m);
            for (int j = 0; j < m; ++j) {
                matrix[i][j] = 0;
            }
        }
    }

    T get(int n, int m) const {
        return matrix[n][m];
    }

    void set(int n, int m, T value) {
        matrix[n][m] = value;
    }

    int *getSize() const {
        return new int[2]{n, m};
    }

    Matrix<T> &operator=(const Matrix<T> &other) {
        this->n = other.n;
        this->m = other.m;
        this->matrix = other.matrix;
        return *this;
    }

    Matrix<T> operator-(Matrix<T> &other) {
        if (this->n != other.n || this->m != other.m) {
            throw std::invalid_argument("Error: the dimensional problem occurred\n");
        }
        Matrix<T> output = Matrix(n, m);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                output.matrix[i][j] = this->matrix[i][j] - other.matrix[i][j];
            }
        }
        return output;
    }

    Matrix<T> operator*(Matrix<T> &other) {
        if (this->m != other.n) {
            throw std::invalid_argument("Error: the dimensional problem occurred\n");
        }
        Matrix<T> output = Matrix<T>(n, other.m);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < other.m; ++j)
                for (int k = 0; k < m; ++k) {
                    output.matrix[i][j] += this->matrix[i][k] * other.matrix[k][j];
                }
        return output;
    }

    Matrix<T> transpose() {
        Matrix<T> output = Matrix(m, n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < m; ++j) {
                output.matrix[j][i] = matrix[i][j];
            }
        return output;
    }

    Matrix<T> operator+(Matrix<T> &other) {
        if (this->n != other.n || this->m != other.m) {
            throw std::invalid_argument("Error: the dimensional problem occurred\n");
        }
        Matrix output = Matrix<T>(n, m);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                output.matrix[i][j] = this->matrix[i][j] + other.matrix[i][j];
            }
        }
        return output;
    }


    // Helper function to find the index of max absolute value in the column
    int find_max_abs(int col) {
        int max_abs_idx = col;
        double max_abs_val = abs(matrix[col][col]);
        for (int i = col + 1; i < n; i++) {
            double abs_val = abs(matrix[i][col]);
            if (abs_val > max_abs_val) {
                max_abs_val = abs_val;
                max_abs_idx = i;
            }
        }
        return max_abs_idx;
    }

    // Helper function to swap rows i and j
    void swap_rows(int i, int j) {
        swap(matrix[i], matrix[j]);
    }

    // Helper function to print matrix
    void print_mat(string title) {
        cout << title << endl;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                printf("%.4f ", matrix[i][j]);
            }
            cout << endl;
        }
    }
};

template<typename T>
class SquareMatrix : public Matrix<T> {
public:
    SquareMatrix(int n) : Matrix<T>(n, n) {};
};

template<typename T>
class IdentityMatrix : public SquareMatrix<T> {
public:
    IdentityMatrix(int n) : SquareMatrix<T>(n) {
        for (int i = 0; i < n; i++) {
            this->matrix[i][i] = 1;
        }
    };
};

template<typename T>
class EliminationMatrix : public IdentityMatrix<T> {
public:
    EliminationMatrix(Matrix<T> matrix, int row, int col) : IdentityMatrix<T>(matrix.getSize()[0]) {
        double factor = matrix.get(row, col) / matrix.get(col, col);
        if (matrix.get(col, col) != 0)
            this->matrix[row][col] = -factor;
    };
};

template<typename T>
class PermutationMatrix : public IdentityMatrix<T> {
public:
    PermutationMatrix(int n, int i, int j) : IdentityMatrix<T>(n) {
        swap(this->matrix[i], this->matrix[j]);
    };
};


template<typename T>
ostream &operator<<(ostream &out, const Matrix<T> &m) {
    for (int i = 0; i < m.getSize()[0]; i++) {
        for (int j = 0; j < m.getSize()[1]; j++) {
            out << m.get(i, j) << " ";
        }
        out << endl;
    }
    return out;
}

template<typename T>
istream &operator>>(istream &in, Matrix<T> &m) {
    int value;
    for (int i = 0; i < m.getSize()[0]; i++) {
        for (int j = 0; j < m.getSize()[1]; j++) {
            in >> value;
            m.set(i, j, value);
        }
    }
    return in;
}


// Calculate the determinant of the matrix A using Gaussian elimination
template<typename T>
double calc_det(Matrix<T> A) {
    double det = 1.0;
    int step = 1;
    for (int col = 0; col < A.getSize()[0]; col++) {
        // Find row with max absolute value in col
        int max_abs_idx = A.find_max_abs(col);
        if (max_abs_idx != col) {
            A.swap_rows(max_abs_idx, col);
            det *= -1.0; // Flip sign of determinant
            A.print_mat("step #" + to_string(step++) + ": permutation");
        }

        // Check if determinant is zero
        if (A.get(col, col) == 0.0) {
            det = 0.0;
            break;
        }

        // Gaussian Elimination
        for (int row = col + 1; row < A.getSize()[0]; row++) {
            EliminationMatrix<T> E = EliminationMatrix<T>(A, row, col);
            A = E * A;
            A.print_mat("step #" + to_string(step++) + ": elimination");
        }
    }

    // Calculate determinant
    for (int i = 0; i < A.getSize()[0]; i++) {
        det *= A.get(i, i);
    }

    // Print final matrix and determinant
    cout << "result:" << endl;
    printf("%.2f\n", det);

    return det;
}

template<typename T>
void print_aug(Matrix<T> A, Matrix<T> I, string title) {
    cout << title << endl;
    int n = A.getSize()[0];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%.2f ", A.get(i, j));
        }
        for (int j = 0; j < n; j++) {
            printf("%.2f ", I.get(i, j));
        }
        cout << endl;
    }
}

template<typename T>
Matrix<T> calc_inverse(Matrix<T> matrix) {
    int n = matrix.getSize()[0];
    Matrix<T> I = IdentityMatrix<T>(n);

    int step = 1;
    for (int col = 0; col < n; col++) {
        // Find row with max absolute value in col
        int max_abs_idx = matrix.find_max_abs(col);
        if (max_abs_idx != col) {
            matrix.swap_rows(max_abs_idx, col);
            I.swap_rows(max_abs_idx, col);
        }

        if (matrix.get(col, col) == 0.0) {
            break;
        }

        // Gaussian Elimination
        for (int row = col + 1; row < matrix.getSize()[0]; row++) {
            if (matrix.get(row, col) == 0.0)
                continue;
            EliminationMatrix<T> E = EliminationMatrix<T>(matrix, row, col);
            matrix = E * matrix;
            I = E * I;
        }
    }

    // Backward substitution
    for (int col = n - 1; col >= 0; --col) {
        for (int row = col - 1; row >= 0; --row) {
            double factor = matrix.get(row, col) / matrix.get(col, col);
            for (int k = 0; k < n; k++) {
                matrix.set(row, k, matrix.get(row, k) - factor * matrix.get(col, k));
                I.set(row, k, I.get(row, k) - factor * I.get(col, k));
            }
            matrix.set(row, col, 0.0);
        }
    }


    // Diagonal normalization:
    for (int i = 0; i < n; ++i) {
        double cur = matrix.get(i, i);
        for (int j = 0; j < n; ++j) {
            matrix.set(i, j, matrix.get(i, j) / cur);
            I.set(i, j, I.get(i, j) / cur);
        }
    }

    return I;
}

int main() {
    int m;
    cin >> m;
    vector<double> A_val = vector<double>();
    Matrix<double> B = Matrix<double>(m, 1);
    for (int i = 0; i < m; ++i) {
        double a, b;
        cin >> a >> b;
        A_val.push_back(a);
        B.set(i, 0, b);
    }

    int n;
    cin >> n;
    Matrix<double> A = Matrix<double>(m, n+1);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j <= n; ++j) {
            A.set(i, j, pow(A_val.at(i), j));
        }
    }

    A.print_mat("A:");
    Matrix<double> AtA = A.transpose() * A;
    AtA.print_mat("A_T*A:");
    Matrix<double> AtA_inverse = calc_inverse(AtA);
    AtA_inverse.print_mat("(A_T*A)^-1:");
    Matrix<double> AtB = A.transpose() * B;
    AtB.print_mat("A_T*b:");
    Matrix<double> x = AtA_inverse * AtB;
    x.print_mat("x~:");
    return 0;
}
