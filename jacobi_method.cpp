#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

class ColumnVector;

class Matrix {
protected:
    int n;
    int m;
public:
    Matrix() : n(0), m(0), data(vector<vector<double>>()) {}
    Matrix(int n, int m) : n(n), m(m), data(n, vector<double>(m, 0)) {}

    [[nodiscard]] int getN() const {
        return n;
    }
    [[nodiscard]] int getM() const {
        return m;
    }

    [[nodiscard]] vector<vector<double>> getData() const {
        return data;
    }

    [[nodiscard]] double getValue(int i, int j) const {
        return data[i][j];
    }

    void setValue(int i, int j, double value) {
        data[i][j] = value;
    }

    void setData(const vector<vector<double>>& values) {
        this->data = values;
    }

    void setN(int newN) {
        this->n = newN;
    }

    void setM(int newM) {
        this->m = newM;
    }

    Matrix transpose() const {
        Matrix transposed(m, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                transposed.setValue(j, i, getValue(i, j));
            }
        }
        return transposed;
    }

    virtual Matrix operator+(const Matrix& other) const {
        if (n != other.n || m != other.m) {
            cout << "Error: the dimensional problem occurred" << endl;
            return {};
        }

        Matrix resultAddition(n, m);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                resultAddition.data[i][j] = data[i][j] + other.data[i][j];
            }
        }
        return resultAddition;
    }

    virtual Matrix operator-(const Matrix& other) const {
        if (n != other.n || m != other.m) {
            cout << "Error: the dimensional problem occurred" << endl;
            return {};
        }

        Matrix resultSubtraction(n, m);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                resultSubtraction.data[i][j] = data[i][j] - other.data[i][j];
            }
        }
        return resultSubtraction;
    }

    virtual Matrix operator*(const Matrix& other) const {
        if (m != other.n) {
            cout << "Error: the dimensional problem occurred" << endl;
            return {};
        }

        Matrix resultMultiplication(n, other.m);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < other.m; j++) {
                for (int k = 0; k < m; k++) {
                    resultMultiplication.data[i][j] += data[i][k] * other.data[k][j];
                }
            }
        }
        return resultMultiplication;
    }

    Matrix& operator=(const Matrix& other) {
        if (this == &other) {
            return *this;
        }

        data.clear();
        n = other.n;
        m = other.m;

        data = other.data;
        return *this;
    }

    friend istream& operator>>(istream& in, Matrix& matrix) {
        in >> matrix.n;
        in >> matrix.m;
        matrix.data.resize(matrix.n, vector<double>(matrix.m));
        for (int i = 0; i < matrix.n; i++) {
            for (int j = 0; j < matrix.m; j++) {
                in >> matrix.data[i][j];
            }
        }
        return in;
    }

    friend ostream& operator<<(ostream& out, const Matrix& matrix) {
        for (int i = 0; i < matrix.n; i++) {
            for (int j = 0; j < matrix.m; j++) {
                out << fixed << setprecision(4) << matrix.data[i][j] << " ";
            }
            out << endl;
        }
        return out;
    }
    vector<vector<double>> data;
};

class SquareMatrix : public Matrix{
public:
    SquareMatrix() : Matrix() {}
    explicit SquareMatrix(int n) : Matrix(n, n) {}

    SquareMatrix transpose() const {
        SquareMatrix transposed(n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                transposed.setValue(j, i, getValue(i, j));
            }
        }
        return transposed;
    }

    Matrix operator+(const Matrix& other) const override {
        if (dynamic_cast<const SquareMatrix*>(&other)) {
            const auto& squareMatrix = dynamic_cast<const SquareMatrix&>(other);
            return *this + squareMatrix;
        } else {
            return Matrix::operator+(other);
        }
    }

    Matrix operator-(const Matrix& other) const override {
        if (dynamic_cast<const SquareMatrix*>(&other)) {
            const auto& squareMatrix = dynamic_cast<const SquareMatrix&>(other);
            return *this - squareMatrix;
        } else {
            return Matrix::operator-(other);
        }
    }

    Matrix operator*(const Matrix& other) const override {
        if (dynamic_cast<const SquareMatrix*>(&other)) {
            const auto& squareMatrix = dynamic_cast<const SquareMatrix&>(other);
            return *this * squareMatrix;
        } else {
            return Matrix::operator*(other);
        }
    }

    SquareMatrix& operator=(const SquareMatrix& other) {
        if (this == &other) {
            return *this;
        }

        data.clear();
        n = other.n;
        m = other.m;

        data = other.data;
        return *this;
    }

    friend istream& operator>>(istream& in, SquareMatrix& matrix) {
        in >> matrix.n;
        matrix.m = matrix.n;
        matrix.data.resize(matrix.n, vector<double>(matrix.m));

        for (int i = 0; i < matrix.n; i++) {
            for (int j = 0; j < matrix.m; j++) {
                in >> matrix.data[i][j];
            }
        }
        return in;
    }

    friend ostream& operator<<(ostream& out, const SquareMatrix& matrix) {
        for (int i = 0; i < matrix.n; i++) {
            for (int j = 0; j < matrix.n; j++) {
                out << fixed << setprecision(4) << matrix.data[i][j] << " ";
            }
            out << endl;
        }
        return out;
    }
};

class Identity : public SquareMatrix {
public:
    Identity() : SquareMatrix() {}
    explicit Identity(int n) : SquareMatrix(n) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                if (i == j) {
                    data[i][j] = 1;
                } else {
                    data[i][j] = 0;
                }
            }
        }
    }

    Identity transpose() {
        return *this;
    }

    virtual SquareMatrix operator*(const SquareMatrix& other) const {
        if (n != other.getN()) {
            cout << "Error: the dimensional problem occurred" << endl;
            return {};
        }
        return other;
    }

    Matrix operator*(const Matrix& other) const {
        if (n != other.getN()) {
            cout << "Error: the dimensional problem occurred" << endl;
            return {};
        }
        return other;
    }

    friend istream& operator>>(istream& in, Identity& matrix) {
        in >> matrix.n;
        matrix.m = matrix.n;
        matrix.data.resize(matrix.n, vector<double>(matrix.m));

        for (int i = 0; i < matrix.n; i++) {
            for (int j = 0; j < matrix.m; j++) {
                if (i == j) {
                    matrix.data[i][j] = 1;
                } else {
                    matrix.data[i][j] = 0;
                }
            }
        }
        return in;
    }

    friend ostream& operator<<(ostream& out, const Identity& matrix) {
        for (int i = 0; i < matrix.n; i++) {
            for (int j = 0; j < matrix.n; j++) {
                out << fixed << setprecision(4) << matrix.data[i][j] << " ";
            }
            out << endl;
        }
        return out;
    }
};

class ColumnVector : public Matrix {
public:
    ColumnVector() : Matrix() {}
    explicit ColumnVector(int n) : Matrix(n, 1) {}

    ColumnVector& operator=(Matrix other) {
        if (this != &other) {
            this->data = other.getData();
        }
        return *this;
    }

    friend istream& operator>>(istream& in, ColumnVector& matrix) {
        in >> matrix.n;
        matrix.m = 1;
        matrix.data.resize(matrix.n, vector<double>(matrix.m));

        for (int i = 0; i < matrix.n; i++) {
            in >> matrix.data[i][0];
        }
        return in;
    }

    friend ostream& operator<<(ostream& out, const ColumnVector& matrix) {
        for (int i = 0; i < matrix.n; i++) {
            out << matrix.data[i][0] << endl;
        }
        return out;
    }

};

class EliminationMatrix : public Identity {
protected:
    int row, column;
public:
    EliminationMatrix() : Identity() {row = 0, column = 0;}
    EliminationMatrix(int n, int row, int column) : Identity(n) {
        this->row = row - 1;
        this->column = column - 1;
    }
    EliminationMatrix(int n, int row, int column, const Matrix& matrix) : Identity(n) {
        this->row = row - 1;
        this->column = column - 1;
        this->getEliminationMatrix(matrix);
    }

    void getEliminationMatrix(const Matrix& matrix) {
        double value = 0;
        vector<double> rowElimination = data[row];
        vector<double> columnMatrix;
        columnMatrix.resize(n);
        for (int i = 0; i < n; i++) {
            columnMatrix[i] = matrix.getData()[i][column];
        }
        for (int i = 0; i < n; i++) {
            if (column != i) {
                value += rowElimination[i] * columnMatrix[i];
            }
        }
        if (columnMatrix[column] != 0) {
            double result = - value / columnMatrix[column];
            data[row][column] = result;
        } else {
            return;
        }
    }

    ColumnVector operator*(const ColumnVector& other) const {
        if (n != other.getN()) {
            cout << "Error: the dimensional problem occurred" << endl;
            return {};
        }

        ColumnVector resultElimination(n);
        for (int i = 0; i < n; i++) {
            for (int k = 0; k < m; k++) {
                resultElimination.setValue(i, 0, resultElimination.getValue(i, 0) + data[i][k] * other.getValue(k, 0));
            }
        }
        return resultElimination;
    }

    SquareMatrix operator*(const SquareMatrix& other) const {
        if (n != other.getN()) {
            cout << "Error: the dimensional problem occurred" << endl;
            return {};
        }

        SquareMatrix resultElimination(n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < other.getM(); j++) {
                for (int k = 0; k < m; k++) {
                    resultElimination.setValue(i, j, resultElimination.getValue(i, j) + data[i][k] * other.getValue(k, j));
                }
            }
        }
        return resultElimination;
    }

    friend istream& operator>>(istream& in, EliminationMatrix& matrix) {
        in >> matrix.n;
        in >> matrix.row >> matrix.column;
        matrix.m = matrix.n;
        matrix.data.resize(matrix.n, vector<double>(matrix.m));
        matrix.row -= 1;
        matrix.column -= 1;

        for (int i = 0; i < matrix.n; i++) {
            for (int j = 0; j < matrix.m; j++) {
                if (i == j) {
                    matrix.data[i][j] = 1;
                } else {
                    matrix.data[i][j] = 0;
                }
            }
        }
        return in;
    }

    friend ostream& operator<<(ostream& out, const EliminationMatrix& matrix) {
        for (int i = 0; i < matrix.n; i++) {
            for (int j = 0; j < matrix.n; j++) {
                out << fixed << setprecision(4) << matrix.data[i][j] << " ";
            }
            out << endl;
        }
        return out;
    }
};

class PermutationMatrix : public EliminationMatrix {
public:
    PermutationMatrix() : EliminationMatrix() {}
    explicit PermutationMatrix(int n, int row, int column) : EliminationMatrix(n, row, column) {
        this->swap();
    }

    void swap() {
        vector<double> temp = data[row];
        data[row] = data[column];
        data[column] = temp;
    }

    friend istream& operator>>(istream& in, PermutationMatrix& matrix) {
        in >> matrix.n;
        in >> matrix.row >> matrix.column;
        matrix.m = matrix.n;
        matrix.data.resize(matrix.n, vector<double>(matrix.m));

        for (int i = 0; i < matrix.n; i++) {
            for (int j = 0; j < matrix.m; j++) {
                if (i == j) {
                    matrix.data[i][j] = 1;
                } else {
                    matrix.data[i][j] = 0;
                }
            }
        }
        matrix.swap();
        return in;
    }

    friend ostream& operator<<(ostream& out, const PermutationMatrix& matrix) {
        for (int i = 0; i < matrix.n; i++) {
            for (int j = 0; j < matrix.n; j++) {
                out << matrix.data[i][j] << " ";
            }
            out << endl;
        }
        return out;
    }
};

class RowVector : public Matrix {
public:
    RowVector() : Matrix() {}
    explicit RowVector(int m) : Matrix(1, m) {}

    friend istream& operator>>(istream& in, RowVector& matrix) {
        in >> matrix.m;
        matrix.n = 1;
        matrix.data.resize(matrix.n, vector<double>(matrix.m));

        for (int i = 0; i < matrix.n; i++) {
            for (int j = 0; j < matrix.m; j++) {
                in >> matrix.data[i][j];
            }
        }
        return in;
    }
};

class AugmentedMatrixInverse : public SquareMatrix {
private:
    SquareMatrix inverse;
public:
    AugmentedMatrixInverse() : SquareMatrix(), inverse() {}
    explicit AugmentedMatrixInverse(int n) : SquareMatrix(n), inverse(n) {
        Identity identity(n);
        inverse.setData(identity.getData());
    }

    friend istream& operator>>(istream& in, AugmentedMatrixInverse& matrix) {
        in >> matrix.n;
        matrix.m = matrix.n;
        matrix.data.resize(matrix.n, vector<double>(matrix.m));
        Identity identity(matrix.n);
        matrix.inverse.setData(identity.getData());

        for (int i = 0; i < matrix.n; i++) {
            for (int j = 0; j < matrix.m; j++) {
                in >> matrix.data[i][j];
            }
        }
        return in;
    }

    friend ostream& operator<<(ostream& out, const AugmentedMatrixInverse& matrix) {
        for (int i = 0; i < matrix.n; i++) {
            for (int j = 0; j < 2 * matrix.n; j++) {
                if (j >= matrix.n) {
                    out << fixed << setprecision(4) << matrix.inverse.getValue(i, j - matrix.n) << " ";
                } else {
                    out << fixed << setprecision(4) << matrix.data[i][j] << " ";
                }
            }
            out << endl;
        }
        return out;
    }

    AugmentedMatrixInverse& operator=(const SquareMatrix& other) {
        if (this != &other) {
            this->data = other.getData();
        }
        return *this;
    }

    void setRow(int j, double ratio) {
        for (int i = 0; i < n; i++) {
            inverse.setValue(j, i, inverse.getValue(j, i) * ratio);
        }
    }

    SquareMatrix getInverse() {
        return inverse;
    }

    void setInverse(const SquareMatrix& matrix) {
        this->inverse = matrix;
    }

    double determinant() {
        double d = 1;
        int step = 0;

        SquareMatrix m(n);
        vector<vector<double>> temp(data);
        m.setData(temp);
        for (int i = 0; i < n; i++) {
            int pivot = i;
            for (int j = i + 1; j < n; j++) {
                if (abs(m.getValue(j, i) ) > abs(m.getValue(pivot, i) )) {
                    pivot = j;
                }
            }

            if (pivot != i) {
                PermutationMatrix p(n, pivot + 1, i + 1);
                m = p * m;
                d *= (-1);
            }

            for (int j = i + 1; j < n; j++) {
                EliminationMatrix e(n, j + 1, i + 1, m);
                m = e * m;
            }
        }

        for (int i = 0; i < n; i++) {
            d *= m.getValue(i, i);
        }
        stringstream ss;
        ss << fixed << setprecision(4) << d;
        double result;
        ss >> result;

        return result;
    }
};

SquareMatrix inverseWithoutPrint(const Matrix& matrix) {
    int n = matrix.getN();
    AugmentedMatrixInverse m(n);
    m.setData(matrix.getData());

    if (m.determinant() == 0) {
        cout << "Error: matrix A is singular" << endl;
        return {};
    }

    for (int i = 0; i < n; i++) {
        int pivot = i;
        for (int j = i + 1; j < n; j++) {
            if (abs(m.getValue(j, i) ) > abs(m.getValue(pivot, i) )) {
                pivot = j;
            }
        }

        if (pivot != i) {
            PermutationMatrix p(n, pivot + 1, i + 1);
            m = p * m;
            m.setInverse(p * m.getInverse());
        }

        for (int j = i + 1; j < n; j++) {
            if (m.getValue(j, i) != 0) {
                EliminationMatrix e(n, j + 1, i + 1, m);
                m = e * m;
                m.setInverse(e * m.getInverse());
            }
        }
    }

    for (int j = n - 1; j >= 0; j--) {
        for (int i = j - 1; i >= 0; i--) {
            if (m.getValue(i, j) != 0) {
                EliminationMatrix e(n, i + 1, j + 1, m);
                m = e * m;
                m.setInverse(e * m.getInverse());
            }
        }
    }

    for (int i = 0; i < n; i++) {
        double ratio = 1 / m.getValue(i, i);
        m.setValue(i, i, m.getValue(i, i) * ratio);
        m.setRow(i, ratio);
    }

    return m.getInverse();
}

class AugmentedMatrixJM : public Matrix {
private:
    ColumnVector b;
public:
    AugmentedMatrixJM() : Matrix(), b() {}
    explicit AugmentedMatrixJM(int n) : Matrix(n, 2), b(ColumnVector(n)) {}

    friend istream& operator>>(istream& in, AugmentedMatrixJM& matrix) {
        in >> matrix.n;
        matrix.m = matrix.n;
        matrix.data.resize(matrix.n, vector<double>(matrix.m));
        matrix.b = ColumnVector(matrix.n);
        matrix.b.setN(matrix.n);
        matrix.b.setM(1);

        for (int i = 0; i < matrix.n; i++) {
            for (int j = 0; j < matrix.m; j++) {
                in >> matrix.data[i][j];
            }
        }
        int newN;
        in >> newN;
        for (int i = 0; i < matrix.n; i++) {
            in >> matrix.b.data[i][0];
        }
        return in;
    }

    friend ostream& operator<<(ostream& out, const AugmentedMatrixJM& matrix) {
        for (int i = 0; i < matrix.getN(); i++) {
            for (int j = 0; j < matrix.getM(); j++) {
                out << fixed << setprecision(4) << matrix.data[i][j] << " ";
            }
            out << endl;
        }
        return out;
    }

    ColumnVector getB() const {
        return b;
    }

    bool diagonalPredominance() {
        for (int i = 0; i < n; i++) {
            double diagonalValue = 0;
            double rowSum = 0;
            for (int j = 0; j < m; j++) {
                if (i == j) {
                    diagonalValue = data[i][j];
                } else {
                    rowSum += data[i][j];
                }
            }
            if (diagonalValue <= rowSum) {
                return false;
            }
        }
        return true;
    }
};

void jacobiMethod(AugmentedMatrixJM A, double eps) {
    cout << A;
    if (!A.diagonalPredominance()) {
        cout << "The method is not applicable" << endl;
        return;
    }

    AugmentedMatrixJM alpha = A;
    for (int i = 0; i < A.getN(); i++) {
        for (int j = 0; j < A.getM(); j++) {
            if (i == j) {
                alpha.setValue(i, j, 0);
            } else {
                alpha.setValue(i, j, - A.getValue(i, j) / A.getValue(i, i));
            }
        }
    }
    cout << "alpha:" << endl;
    cout << alpha;
    ColumnVector beta = A.getB();
    for (int i = 0; i < A.getN(); i++) {
        double first = A.getValue(i, i);
        beta.setValue(i, 0, A.getB().getValue(i, 0) / first);
    }
    cout << "beta:" << endl;
    cout << beta;

    int i = 1;
    ColumnVector x_0 = beta;
    ColumnVector x_i = x_0;
    while (true) {
        cout << "x(" << i << "):" << endl;
        x_i = alpha * x_0 + beta;
        cout << x_i;
        double sm = 0;
        for (int i = 0; i < A.getN(); i++) {
            sm += pow(x_i.getValue(i, 0) - x_0.getValue(i, 0), 2);
        }

        sm = sqrt(sm);
        if (sm < eps) {
            cout << "e: " << setprecision(4) << sm << endl;
            break;
        }
        cout << "e: " << setprecision(4) << sm << endl;
        x_0 = x_i;
        i += 1;
    }
    cout << "x~:" << endl;
    cout << x_i;
}

int main() {
    AugmentedMatrixJM A;
    double eps;
    cin >> A >> eps;
    jacobiMethod(A, eps);

    return 0;
}
