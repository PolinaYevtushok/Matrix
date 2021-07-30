#pragma once
#include <iostream>
#include <type_traits>

template<typename T, size_t N, size_t M>
class Matrix
{
public:
	Matrix();
	Matrix(const Matrix& m);
	Matrix(Matrix&& m);
	Matrix& operator=(const Matrix& m);
	~Matrix();
	T& operator()(size_t row, size_t colum) const;
	Matrix& operator++();
	Matrix& operator--();
	Matrix& operator-();
	void fillValues (T value);
	void print() const;
	void setValueAt(size_t row, size_t colum, T value);
	T getValueAt(size_t row, size_t colum) const;
	size_t GetRowsNumber() const { return m_rows; }
	size_t GetColumsNumber() const { return m_colums; }
	T** GetMatrix() const { return m_matrix; }

private:
	const size_t m_rows = N;
	const size_t m_colums = M;
	T** m_matrix;
};

template<typename T, size_t N, size_t M>
Matrix<T, N, M>::Matrix()
{
	m_matrix = new T* [m_rows];
	for (int i = 0; i < m_colums; ++i)
		m_matrix[i] = new T[m_colums];
	fillValues(0);
}

template<typename T, size_t N, size_t M>
Matrix<T, N, M>::Matrix(const Matrix& m)
{
	m_matrix = new T * [m_rows];
	for (int i = 0; i < m_colums; ++i)
		m_matrix[i] = new T[m_colums];
	for (int i = 0; i < m_rows; ++i)
	{
		for (int j = 0; j < m_colums; ++j)
			m_matrix[i][j] = m.GetMatrix()[i][j];
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
Matrix<T, N, M>& Matrix<T, N, M>::operator=(const Matrix& m)
{
	if (&m == this)
		return *this;
	else
	{
		//delete
		for (int i = 0; i < m_colums; ++i)
			delete[] m_matrix[i];
		delete[] m_matrix;
	}
	//create new with new parametrs
	m_matrix = new T * [m_rows];
	for (int i = 0; i < m_colums; ++i)
		m_matrix[i] = new T[m_colums];
	// fill new matrix
	for (int i = 0; i < m_rows; ++i)
	{
		for (int j = 0; j < m_colums; ++j)
			m_matrix[i][j] = m.GetMatrix()[i][j];
	}
	return *this;
}

template<typename T, size_t N, size_t M>
Matrix<T, N, M>::~Matrix()
{
	if (m_matrix)
	{
		for (int i = 0; i < m_colums; ++i)
		{
			if (m_matrix[i])
				delete[] m_matrix[i];
		}			
		delete[] m_matrix;
	}
}

template<typename T, size_t N, size_t M>
T& Matrix<T, N, M>::operator()(size_t row, size_t colum) const
{	
	if (m_rows >= row && m_colums >= colum )
		return GetMatrix()[row - 1][colum - 1];
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

template<typename T, typename U, size_t N, size_t M, size_t K, size_t L>
auto operator+(const Matrix<T, N, M>& m1, const Matrix<U, K, L>& m2)
{	
	if (m1.GetRowsNumber() == m2.GetRowsNumber() && m1.GetColumsNumber() == m2.GetColumsNumber())
	{
		auto t = m1.getValueAt(1, 1) + m2.getValueAt(1, 1);
		Matrix<decltype(t), N, M> res;
		for (int i = 0; i < res.GetRowsNumber(); ++i)
		{
			for (int j = 0; j < res.GetColumsNumber(); ++j)
				res.GetMatrix()[i][j] = m1.GetMatrix()[i][j] + m2.GetMatrix()[i][j];
		}
		return res;
	}
	else
		throw std::length_error("Matrices sizes mismatch");
}

template<typename T, typename U, size_t N, size_t M, size_t K, size_t L>
auto operator-(const Matrix<T, N, M>& m1, const Matrix<U, K, L>& m2)
{
	if (m1.GetRowsNumber() == m2.GetRowsNumber() && m1.GetColumsNumber() == m2.GetColumsNumber())
	{
		auto t = m1.getValueAt(1, 1) - m2.getValueAt(1, 1);
		Matrix<decltype(t), N, M> res;
		for (int i = 0; i < res.GetRowsNumber(); ++i)
		{
			for (int j = 0; j < res.GetColumsNumber(); ++j)
				res.GetMatrix()[i][j] = m1.GetMatrix()[i][j] - m2.GetMatrix()[i][j];
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
	for (int i = 0; i < res.GetRowsNumber(); ++i)
	{	
			for (int j = 0; j < res.GetColumsNumber(); ++j)				
				res.GetMatrix()[i][j] = m.GetMatrix()[i][j] * scalar;
	}
	return res;
}

template<typename T, typename U, size_t N, size_t M, size_t K, size_t L>
auto operator*(const Matrix<T, N, M>& m1, const Matrix<U, K, L>& m2)
{
	if (m1.GetColumsNumber() == m2.GetRowsNumber())
	{
		auto t = m1.getValueAt(1, 1) * m2.getValueAt(1, 1);
		Matrix<decltype(t), N, L> res;
		for (int i = 0; i < m1.GetRowsNumber(); ++i)
		{
			for (int j = 0; j < m2.GetColumsNumber(); ++j)
			{
				decltype(t) sum = 0;
				for (int k = 0; k < m2.GetRowsNumber(); ++k)
					sum += m1.GetMatrix()[i][k] * m2.GetMatrix()[k][j];
				res.GetMatrix()[i][j] = sum;
			}
		}
		return res;
	}
	else
		throw std::length_error("Matrices sizes mismatch");
}
