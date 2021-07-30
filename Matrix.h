#pragma once
#include <iostream>
#include <type_traits>
#include <functional>

template<typename T, size_t N, size_t M>
class Matrix
{
public:
	Matrix();
	Matrix(const Matrix& m);
	Matrix(Matrix&& m);
	template<typename U, size_t K, size_t L>
	Matrix& operator=(const Matrix<U, K, L>& m);
	template<typename U, size_t K, size_t L>
	Matrix& operator=(Matrix<U, K, L>&& m);
	~Matrix();
	T& operator()(size_t row, size_t colum) const;
	Matrix& operator++();
	Matrix& operator--();
	Matrix& operator-();
	void fillValues (T value);
	void print() const;
	void setValueAt(size_t row, size_t colum, T value);
	T getValueAt(size_t row, size_t colum) const;
	size_t getRowsNumber() const { return m_rows; }
	size_t getColumsNumber() const { return m_colums; }
	void setRowsNumber(size_t n) {m_rows = n; }
	void setColumsNumber(size_t n) {m_colums = n; }
	T** getMatrix() const { return m_matrix; }
	void setMatrixNullptr() { m_matrix = nullptr; }
	T** getMoveMatrix() { return std::move(m_matrix); }

protected:
	void matrixAllocation();
	void matrixFree();
private:
	size_t m_rows = N;
	size_t m_colums = M;
	T** m_matrix;
};

template<typename T, size_t N, size_t M>
Matrix<T, N, M>::Matrix()
{
	matrixAllocation();
	fillValues(0);
}

template<typename T, size_t N, size_t M>
Matrix<T, N, M>::Matrix(const Matrix& m):
	m_rows(m.m_rows),
	m_colums(m.m_colums)
{
	matrixAllocation();
	for (int i = 0; i < m_rows; ++i)
	{
		for (int j = 0; j < m_colums; ++j)
			m_matrix[i][j] = m.getMatrix()[i][j];
	}
}

template<typename T, size_t N, size_t M>
Matrix<T, N, M>::Matrix(Matrix&& m):
	m_rows(m.m_rows),
	m_colums(m.m_colums),
	m_matrix(std::move(m.m_matrix))
{
	m.m_matrix = nullptr;
}

template<typename T, size_t N, size_t M>
template<typename U, size_t K, size_t L>
Matrix<T, N, M>& Matrix<T, N, M>::operator=(const Matrix<U, K, L>& m)
{
	matrixFree();
	setRowsNumber(m.getRowsNumber());
	setColumsNumber(m.getColumsNumber()); 
	matrixAllocation();
	// fill new matrix
	for (int i = 0; i < m_rows; ++i)
	{
		for (int j = 0; j < m_colums; ++j)
			m_matrix[i][j] = m.getMatrix()[i][j];
	}
	return *this;
}

template<typename T, size_t N, size_t M>
template<typename U, size_t K, size_t L>
Matrix<T, N, M>& Matrix<T, N, M>::operator=(Matrix<U, K, L>&& m)
{
	matrixFree();
	setRowsNumber(m.getRowsNumber());
	setColumsNumber(m.getColumsNumber());
	m_matrix = m.getMoveMatrix();
	m.setMatrixNullptr();
	return *this;
}

template<typename T, size_t N, size_t M>
Matrix<T, N, M>::~Matrix()
{
	matrixFree();
}

template<typename T, size_t N, size_t M>
T& Matrix<T, N, M>::operator()(size_t row, size_t colum) const
{	
	if (m_rows >= row && m_colums >= colum )
		return getMatrix()[row - 1][colum - 1];
	else 
		throw std::out_of_range("Matrices sizes mismatch");
}

template<typename T, size_t N, size_t M>
Matrix<T, N, M>& Matrix<T, N, M>::operator++()
{
	for (int i = 0; i < m_rows; ++i)
	{
		for (int j = 0; j < m_colums; ++j)
		{
			++m_matrix[i][j];
		}
	}
	return *this;
}

template<typename T, size_t N, size_t M>
Matrix<T, N, M>& Matrix<T, N, M>::operator--()
{
	for (int i = 0; i < m_rows; ++i)
	{
		for (int j = 0; j < m_colums; ++j)
		{
			--m_matrix[i][j];
		}
	}
	return *this;
}

template<typename T, size_t N, size_t M>
Matrix<T, N, M>& Matrix<T, N, M>::operator-()
{
	for (int i = 0; i < m_rows; ++i)
	{
		for (int j = 0; j < m_colums; ++j)
		{
			m_matrix[i][j] *= -1;
		}
	}
	return *this;
}


template<typename T, size_t N, size_t M>
void Matrix<T, N, M>::fillValues(T value)
{
	for (int i = 0; i < m_rows; ++i)
	{
		for (int j = 0; j < m_colums; ++j)
			m_matrix[i][j] = value;
	}
}

template<typename T, size_t N, size_t M>
void Matrix<T, N, M>::print() const
{
	if (m_matrix)
	{
		for (int i = 0; i < m_rows; ++i)
		{
			for (int j = 0; j < m_colums; ++j)
				std::cout << m_matrix[i][j] << " ";
			std::cout << std::endl;
		}
	}
}

template<typename T, size_t N, size_t M>
void Matrix<T, N, M>::setValueAt(size_t row, size_t colum, T value)
{
	if (m_rows >= row  && m_colums >= colum)
		m_matrix[row - 1][colum - 1] = value;
	else
		throw std::out_of_range("Matrices sizes mismatch");
}

template<typename T, size_t N, size_t M>
T Matrix<T, N, M>::getValueAt(size_t row, size_t colum) const
{
	if (m_rows >= row && m_colums >= colum)
		return m_matrix[row - 1][colum - 1];
	else
		throw std::out_of_range("Matrices sizes mismatch");
}

template<typename T, size_t N, size_t M>
void Matrix<T, N, M>::matrixAllocation()
{
	m_matrix = new T * [m_rows];
	for (int i = 0; i < m_rows; ++i)
		m_matrix[i] = new T[m_colums];
}

template<typename T, size_t N, size_t M>
void Matrix<T, N, M>::matrixFree()
{
	if (m_matrix)
	{
		for (int i = 0; i < m_rows; ++i)
		{
			if (m_matrix[i])
				delete[] m_matrix[i];
		}
		delete[] m_matrix;
	}
}

template<typename T, typename U, size_t N, size_t M, size_t K, size_t L>
auto operator+(const Matrix<T, N, M>& m1, const Matrix<U, K, L>& m2)
{	
	if (m1.getRowsNumber() == m2.getRowsNumber() && m1.getColumsNumber() == m2.getColumsNumber())
	{
		auto t = m1.getValueAt(1, 1) + m2.getValueAt(1, 1);
		Matrix<decltype(t), N, M> res;
		for (int i = 0; i < res.getRowsNumber(); ++i)
		{
			for (int j = 0; j < res.getColumsNumber(); ++j)
				res.getMatrix()[i][j] = m1.getMatrix()[i][j] + m2.getMatrix()[i][j];
		}
		return res;
	}
	else
		throw std::length_error("Matrices sizes mismatch");
}

template<typename T, typename U, size_t N, size_t M, size_t K, size_t L>
auto operator-(const Matrix<T, N, M>& m1, const Matrix<U, K, L>& m2)
{
	if (m1.getRowsNumber() == m2.getRowsNumber() && m1.getColumsNumber() == m2.getColumsNumber())
	{
		auto t = m1.getValueAt(1, 1) - m2.getValueAt(1, 1);
		Matrix<decltype(t), N, M> res;
		for (int i = 0; i < res.getRowsNumber(); ++i)
		{
			for (int j = 0; j < res.getColumsNumber(); ++j)
				res.getMatrix()[i][j] = m1.getMatrix()[i][j] - m2.getMatrix()[i][j];
		}
		return res;
	}
	else
		throw std::length_error("Matrices sizes mismatch");
}

template<typename T, typename U,  size_t N, size_t M>
auto operator*(const Matrix<T, N, M>& m, U scalar)
{
	auto t = m.getValueAt(1, 1) * scalar;
	Matrix<decltype(t), N, M> res;
	for (int i = 0; i < res.getRowsNumber(); ++i)
	{	
			for (int j = 0; j < res.getColumsNumber(); ++j)				
				res.getMatrix()[i][j] = m.getMatrix()[i][j] * scalar;
	}
	return res;
}

template<typename T, typename U, size_t N, size_t M, size_t K, size_t L>
auto operator*(const Matrix<T, N, M>& m1, const Matrix<U, K, L>& m2)
{
	if (m1.getColumsNumber() == m2.getRowsNumber())
	{
		auto t = m1.getValueAt(1, 1) * m2.getValueAt(1, 1);
		Matrix<decltype(t), N, L> res;
		for (int i = 0; i < m1.getRowsNumber(); ++i)
		{
			for (int j = 0; j < m2.getColumsNumber(); ++j)
			{
				decltype(t) sum = 0;
				for (int k = 0; k < m2.getRowsNumber(); ++k)
					sum += m1.getMatrix()[i][k] * m2.getMatrix()[k][j];
				res.getMatrix()[i][j] = sum;
			}
		}
		return res;
	}
	else
		throw std::length_error("Matrices sizes mismatch");
}
