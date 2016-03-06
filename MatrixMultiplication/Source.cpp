#include <iostream>
#include <omp.h>
#include <ctime>
#include <string>
#include <cstring>

using namespace std;

class Matrix
{
private:
	int rows, columns;
	double * data;
public:
	Matrix(int rows, int columns)
		: rows(rows), columns(columns)
	{
		this->data = new double[columns * rows];
	}

	Matrix(double *data, int rows, int columns)
		: data(data), rows(rows), columns(columns)
	{
	}

	~Matrix()
	{
		delete[] this->data;
	}

	int GetRows() { return rows; }
	int GetColumns() { return columns; }
	double Get(int row, int col) { return data[Index(row, col)]; }
	void Set(int row, int col, double value) { data[Index(row, col)] = value; }

	inline int Index(int row, int col)
	{
		return row * columns + col;
	}

	Matrix * NormalMultiply(Matrix *other)
	{
		if (this->columns != other->rows)
			throw new exception("not valid matrix to multiply");

		Matrix * result = new Matrix(this->rows, other->columns);

#pragma omp parallel for
		//schedule(static, this->rows / omp_get_num_threads())
		for (int row = 0; row < this->rows; row++)
#pragma omp parallel for
		//schedule(static, other->columns / omp_get_num_threads())
			for (int col = 0; col < other->columns; col++)
			{
				double tmp = 0;
				#pragma omp parallel for reduction(+:tmp)
				for (int i = 0; i < this->columns; i++)
					tmp += this->Get(row, i) * other->Get(i, col);

				result->Set(row, col, tmp);
			}

		return result;
	}

	Matrix* operator * (Matrix * other)
	{
		return this->NormalMultiply(other);
	}

	void Write(ostream & out)
	{
		for (int row = 0; row < rows; row++)
		{
			for (int col = 0; col < columns; col++)
				out << Get(row, col) << " ";
			out << endl;
		}
	}
};

Matrix * RandomMatrix(int rows, int cols)
{
	auto matrix = new Matrix(rows, cols);
	srand(0);
#pragma omp parallel for shared(matrix) schedule(static, cols / omp_get_num_threads())
	for (int row = 0; row < rows; row++)
		for (int col = 0; col < cols; col++)
			matrix->Set(row, col, (rand() % 100) + (rand() % 10) / 10.);

	return matrix;
}

void PrintMeasureResult(string message, double start, double end)
{
	cout << message << ": " << (end - start) << "s" << endl;
}

int main()
{
#ifdef _DEBUG
	//FILE *stream;
	//freopen_s(&stream, "in.txt", "r", stdin);
	//freopen_s(&stream, "out.txt", "w", stdout);
#endif
	omp_set_num_threads(omp_get_max_threads());

	const int SIZE = 1500;
	double start = omp_get_wtime();
	auto a = RandomMatrix(SIZE, SIZE);
	double end = omp_get_wtime();
	PrintMeasureResult("Generating A", start, end);

	start = omp_get_wtime();
	auto b = RandomMatrix(SIZE, SIZE);
	end = omp_get_wtime();
	PrintMeasureResult("Generating B", start, end);


	start = omp_get_wtime();
	auto c = a->NormalMultiply(b);
	end = omp_get_wtime();
	PrintMeasureResult("A*B", start, end);

	//cout << "A" << endl;
	//a->Write(cout);
	//cout << "B" << endl;
	//b->Write(cout);
	//cout << "A*B" << endl;
	//c->Write(cout);
	delete a;
	delete b;
	delete c;
	cin.get();
	return 0;
}