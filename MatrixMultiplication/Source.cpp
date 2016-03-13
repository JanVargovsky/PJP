#include <iostream>
#include <omp.h>
#include <ctime>
#include <string>
#include <cstring>
#include <sstream>

using namespace std;

int USE_STRASSEN = 1;
int USE_PARALLEL = 1;

class Matrix
{
private:
	const int MIN_SIZE_TO_USE_STRASSEN = 128;

	const int rows, columns;
	const int size;
	double * data;

	Matrix * SubMatrixSingle(int rowFrom, int rowTo, int colFrom, int colTo) const
	{
		// TODO: add validation
		Matrix * result = new Matrix(rowTo - rowFrom, colTo - colFrom);

		for (int row = rowFrom, i = 0; row < rowTo; row++, i++)
			for (int col = colFrom, j = 0; col < colTo; col++, j++)
				result->Set(i, j, Get(row, col));

		return result;
	}

	Matrix * SubMatrixParallel(int rowFrom, int rowTo, int colFrom, int colTo) const
	{
		// TODO: add validation
		Matrix * result = new Matrix(rowTo - rowFrom, colTo - colFrom);


		for (int row = rowFrom, i = 0; row < rowTo; row++, i++)
#pragma omp parallel for
			for (int col = colFrom, j = 0; col < colTo; col++, j++)
			{
				result->Set(i, j, Get(row, col));
			}

		return result;
	}

	Matrix * CombineSubMatrices(Matrix & a11, Matrix & a12, Matrix & a21, Matrix & a22) const
	{
		// TODO: add validation
		Matrix * result = new Matrix(a11.rows * 2, a11.columns * 2);
		int rowShift = a11.rows;
		int colShift = a11.columns;

		for (int row = 0; row < a11.rows; row++)
			for (int col = 0; col < a11.columns; col++)
			{
				result->Set(row, col, a11.Get(row, col));
				result->Set(row, col + colShift, a12.Get(row, col));
				result->Set(row + rowShift, col, a21.Get(row, col));
				result->Set(row + rowShift, col, a22.Get(row, col));
			}

		return result;
	}

public:
	Matrix(double *data, int rows, int columns)
		: data(data), rows(rows), columns(columns), size(rows * columns)
	{
	}

	Matrix(int rows, int columns)
		: Matrix(new double[columns * rows]{}, rows, columns)
	{
	}

	~Matrix()
	{
		delete[] this->data;
	}

#pragma region Operators
	Matrix* operator + (Matrix & other) const
	{
		if (USE_PARALLEL)
			return this->AddParallel(other);
		else
			return this->AddSingle(other);
	}

	Matrix* operator - (Matrix & other) const
	{
		if (USE_PARALLEL)
			return this->SubtractParallel(other);
		else
			return this->SubtractParallel(other);
	}

	Matrix* operator * (Matrix & other) const
	{
		if (USE_STRASSEN)
			//cerr << "Strassen used" << endl;
			if (USE_PARALLEL)
				return this->StrassenMultiplyParallel(other);
			else
				return this->StrassenMultiplySingle(other);
		else
			//cerr << "normal used" << endl;
			if (USE_PARALLEL)
				return this->NormalMultiplyParallel(other);
			else
				return this->NormalMultiplySingle(other);
	}
#pragma endregion

#pragma region Getters, setters and helper functions
	int GetRows() const
	{
		return rows;
	}

	int GetColumns() const
	{
		return columns;
	}

	double const Get(int row, int col) const
	{
		return data[Index(row, col)];
	}

	void Set(int row, int col, double value)
	{
		data[Index(row, col)] = value;
	}

	inline int Index(int row, int col) const
	{
		return row * columns + col;
	}
#pragma endregion

#pragma region Add
	Matrix * AddSingle(Matrix & other) const
	{
		if (this->rows != other.rows || this->columns != other.columns)
			throw new exception("matrices has to have same size");

		Matrix * result = new Matrix(this->rows, this->columns);

		for (int row = 0; row < result->rows; row++)
			for (int col = 0; col < result->columns; col++)
				result->Set(row, col, this->Get(row, col) + other.Get(row, col));

		return result;
	}

	Matrix * AddParallel(Matrix & other) const
	{
		if (this->rows != other.rows || this->columns != other.columns)
			throw new exception("matrices has to have same size");

		Matrix * result = new Matrix(this->rows, this->columns);

#pragma omp parallel for
		for (int row = 0; row < result->rows; row++)
#pragma omp parallel for
			for (int col = 0; col < result->columns; col++)
				result->Set(row, col, this->Get(row, col) + other.Get(row, col));

		return result;
	}
#pragma endregion

#pragma region Subtract
	Matrix * SubtractSingle(Matrix & other) const
	{
		if (this->rows != other.rows || this->columns != other.columns)
			throw new exception("matrices has to have same size");

		Matrix * result = new Matrix(this->rows, this->columns);

		for (int row = 0; row < result->rows; row++)
			for (int col = 0; col < result->columns; col++)
				result->Set(row, col, this->Get(row, col) - other.Get(row, col));

		return result;
	}

	Matrix * SubtractParallel(Matrix & other) const
	{
		if (this->rows != other.rows || this->columns != other.columns)
			throw new exception("matrices has to have same size");

		Matrix * result = new Matrix(this->rows, this->columns);

#pragma omp parallel for
		for (int row = 0; row < result->rows; row++)
#pragma omp parallel for
			for (int col = 0; col < result->columns; col++)
				result->Set(row, col, this->Get(row, col) - other.Get(row, col));

		return result;
	}
#pragma endregion

#pragma region NormalMultiply
	Matrix * NormalMultiplySingle(Matrix &other) const
	{
		if (this->columns != other.rows)
			throw new exception("not valid matrix to multiply");

		Matrix * result = new Matrix(this->rows, other.columns);

		for (int row = 0; row < this->rows; row++)
			for (int col = 0; col < other.columns; col++)
			{
				double tmp = 0;
				for (int i = 0; i < this->columns; i++)
					tmp += this->Get(row, i) * other.Get(i, col);

				result->Set(row, col, tmp);
			}

		return result;
	}

	Matrix * NormalMultiplyParallel(Matrix &other) const
	{
		if (this->columns != other.rows)
			throw new exception("not valid matrix to multiply");

		Matrix * result = new Matrix(this->rows, other.columns);

#pragma omp parallel for
		//schedule(static, this->rows / omp_get_num_threads())
		for (int row = 0; row < this->rows; row++)
#pragma omp parallel for
			//schedule(static, other->columns / omp_get_num_threads())
			for (int col = 0; col < other.columns; col++)
			{
				double tmp = 0;
#pragma omp parallel for reduction(+:tmp)
				for (int i = 0; i < this->columns; i++)
					tmp += this->Get(row, i) * other.Get(i, col);

				result->Set(row, col, tmp);
			}

		return result;
	}

#pragma endregion

#pragma region StrassenMultiply
	Matrix * StrassenMultiplySingle(Matrix & other) const
	{
		// TODO: add validation

		int n = other.rows;
		if (n <= MIN_SIZE_TO_USE_STRASSEN)
			return this->NormalMultiplySingle(other);

		int halfN = n / 2;

		auto a11 = SubMatrixSingle(0, halfN, 0, halfN);
		auto a12 = SubMatrixSingle(0, halfN, halfN, n);
		auto a21 = SubMatrixSingle(halfN, n, 0, halfN);
		auto a22 = SubMatrixSingle(halfN, n, halfN, n);

		auto b11 = SubMatrixSingle(0, halfN, 0, halfN);
		auto b12 = SubMatrixSingle(0, halfN, halfN, n);
		auto b21 = SubMatrixSingle(halfN, n, 0, halfN);
		auto b22 = SubMatrixSingle(halfN, n, halfN, n);

		// TODO: Fix memory leaks !!!
		Matrix* m[7] = {
			(*a11 + *a22)->StrassenMultiplySingle(*(*b11 + *b22)), // m1
			(*a21 + *a22)->StrassenMultiplySingle(*b11), // m2
			(a11)->StrassenMultiplySingle(*(*b12 - *b22)), // m3
			(a22)->StrassenMultiplySingle(*(*b21 - *b11)), // m4
			(*a11 + *a12)->StrassenMultiplySingle(*b22), // m5
			(*a21 - *a11)->StrassenMultiplySingle(*(*b11 + *b12)), // m6
			(*a12 - *a22)->StrassenMultiplySingle(*(*b21 + *b22)), // m7
		};

		// TODO: Fix memory leaks !!!
		auto c11 = *(*(*m[0] + *m[3]) - *m[4]) + *m[6];
		auto c12 = *m[2] + *m[4];
		auto c21 = *m[1] + *m[3];
		auto c22 = *(*(*m[0] - *m[1]) + *m[2]) + *m[5];

		auto result = CombineSubMatrices(*c11, *c12, *c21, *c22);

		delete a11;
		delete a12;
		delete a21;
		delete a22;
		delete b11;
		delete b12;
		delete b21;
		delete b22;
		for (int i = 0; i < 7; i++)
			delete m[i];
		delete c11;
		delete c12;
		delete c21;
		delete c22;

		return result;
	}

	Matrix * StrassenMultiplyParallel(Matrix & other) const
	{
		// TODO: add validation

		int n = other.rows;
		if (n <= MIN_SIZE_TO_USE_STRASSEN)
			return this->NormalMultiplyParallel(other);

		int halfN = n / 2;

		auto a11 = SubMatrixParallel(0, halfN, 0, halfN);
		auto a12 = SubMatrixParallel(0, halfN, halfN, n);
		auto a21 = SubMatrixParallel(halfN, n, 0, halfN);
		auto a22 = SubMatrixParallel(halfN, n, halfN, n);

		auto b11 = SubMatrixParallel(0, halfN, 0, halfN);
		auto b12 = SubMatrixParallel(0, halfN, halfN, n);
		auto b21 = SubMatrixParallel(halfN, n, 0, halfN);
		auto b22 = SubMatrixParallel(halfN, n, halfN, n);

		// m1
		auto m1_a = (*a11 + *a22);
		auto m1_b = (*b11 + *b22);
		// m2
		auto m2_a = *a21 + *a22;
		// m3
		auto m3_b = (*b12 - *b22);
		// m4
		auto m4_b = (*b21 - *b11);
		// m5
		auto m5_a = *a11 + *a12;
		// m6
		auto m6_a = *a21 - *a11;
		auto m6_b = *b11 + *b12;
		// m7
		auto m7_a = *a12 - *a22;
		auto m7_b = *b21 + *b22;

		Matrix* m[7] = {
			m1_a->StrassenMultiplyParallel(*m1_b), // m1
			m2_a->StrassenMultiplyParallel(*b11), // m2
			a11->StrassenMultiplyParallel(*m3_b), // m3
			a22->StrassenMultiplyParallel(*m4_b), // m4
			m5_a->StrassenMultiplyParallel(*b22), // m5
			m6_a->StrassenMultiplyParallel(*m6_b), // m6
			m7_a->StrassenMultiplyParallel(*m7_b), // m7
		};

		auto c11_l = (*m[0] + *m[3]);
		auto c11_r = (*m[4] + *m[6]);

		auto c22_l = (*m[0] - *m[1]);
		auto c22_r = (*m[2] + *m[5]);

		auto c11 = *c11_l - *c11_r;
		auto c12 = *m[2] + *m[4];
		auto c21 = *m[1] + *m[3];
		auto c22 = *c22_l + *c22_r;


		auto result = CombineSubMatrices(*c11, *c12, *c21, *c22);

		delete a11;
		delete a12;
		delete a21;
		delete a22;
		delete b11;
		delete b12;
		delete b21;
		delete b22;

		delete  m1_a;
		delete  m1_b;
		delete  m2_a;
		delete  m3_b;
		delete  m4_b;
		delete  m5_a;
		delete  m6_a;
		delete  m6_b;
		delete  m7_a;
		delete  m7_b;

		for (int i = 0; i < 7; i++)
			delete m[i];

		delete c11_l;
		delete c11_r;
		delete c11;
		delete c12;
		delete c21;
		delete c22_l;
		delete c22_r;
		delete c22;

		return result;
	}
#pragma endregion
};

ostream & operator<<(ostream &out, const Matrix& m)
{
	for (int row = 0; row < m.GetRows(); row++)
	{
		for (int col = 0; col < m.GetColumns(); col++)
			out << m.Get(row, col) << " ";
		out << endl;
	}
	return out;
}


Matrix * RandomMatrix(int rows, int cols)
{
	srand(0);
	auto matrix = new Matrix(rows, cols);

#pragma omp parallel for shared(matrix) schedule(static, cols / omp_get_num_threads())
	for (int row = 0; row < rows; row++)
		for (int col = 0; col < cols; col++)
			//matrix->Set(row, col, (rand() % 100) + (rand() % 10) / 10.);
			matrix->Set(row, col, col);

	return matrix;
}

void PrintMeasureResult(string message, double start, double end)
{
	cout << message << ": " << (end - start) << "s" << endl;
}

void RunTest(Matrix * a, Matrix * b, int strassen, int parallel)
{
	USE_STRASSEN = strassen;
	USE_PARALLEL = parallel;
	auto start = omp_get_wtime();
	auto c = *a * *b;
	auto end = omp_get_wtime();
	stringstream ss;
	ss << "Strassen = " << strassen << " Parallel = " << parallel << " TIME";
	PrintMeasureResult(ss.str(), start, end);
	delete c;
}

int main()
{
#ifdef _DEBUG
	FILE *stream;
	freopen_s(&stream, "in.txt", "r", stdin);
	freopen_s(&stream, "out.txt", "w", stdout);
#endif
	omp_set_num_threads(omp_get_max_threads());

	const int SIZE = 4096;
	double start = omp_get_wtime();
	auto a = RandomMatrix(SIZE, SIZE);
	double end = omp_get_wtime();
	PrintMeasureResult("Generating A", start, end);

	start = omp_get_wtime();
	auto b = RandomMatrix(SIZE, SIZE);
	end = omp_get_wtime();
	PrintMeasureResult("Generating B", start, end);

	RunTest(a, b, 1, 1); // parallel strassen
	//RunTest(a, b, 1, 0); // normal strassen
	//RunTest(a, b, 0, 1); // parallel multiplication
	//RunTest(a, b, 0, 0); // normal multiplication

	delete a;
	delete b;
	cin.get();
	return 0;
}