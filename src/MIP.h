#ifndef LP_HPP
#define LP_HPP

#include <vector>
#include <bitset>

#include "SparseMatrix.h"
#include "Common.h"

#ifndef BOOST_FOUND
	using bitset = std::vector<bool>;
#else
	using bitset = boost::dynamic_bitset;
#endif

template <typename REAL>
struct VectorView
{
	VectorView(const REAL* _array, const size_t* _indices, size_t _size) : coefs(_array), indices(_indices), size(_size) {}

	VectorView(const VectorView<REAL>&) = default;

	const REAL* coefs;
	// use int?
	const size_t* indices;
	size_t size;
};

template<typename REAL>
class MIP
{
public:
	MIP() {};

	MIP(MIP<REAL>&&) = default;

	MIP<REAL>& operator=(MIP<REAL>&& other)
	{
		objective = std::move(other.objective);

		lhs = std::move(other.lhs);
		rhs = std::move(other.rhs);

		ub = std::move(other.ub);
		lb = std::move(other.lb);

		constMatrix = std::move(other.constMatrix);
		constMatrixT = std::move(other.constMatrixT);

		return *this;
	}

	size_t getNCols() const { return constMatrixT.nrows; }

	size_t getNRows() const { return constMatrix.nrows; }

	const std::vector<REAL>& getObj() const { return objective; }

	const std::vector<REAL>& getLB() const { return lb; }

	const std::vector<REAL>& getUB() const { return ub; }

	const std::vector<REAL>& getLHS() const { return lhs; }

	const std::vector<REAL>& getRHS() const { return rhs; }

	VectorView<REAL> getRow(size_t row) const {
		const REAL* coefBegin = constMatrix.coefficients.data() + constMatrix.rowStart[row];
		const size_t* indBegin = constMatrix.indices.data() + constMatrix.rowStart[row];
		size_t size = constMatrix.rowStart[row+1] - constMatrix.rowStart[row];
		VectorView<REAL> view(coefBegin, indBegin, size);
		return view;
	}

	VectorView<REAL> getCol(size_t col) const {
		const REAL* coefBegin = constMatrixT.coefficients.data() + constMatrixT.rowStart[col];
		const size_t* indBegin = constMatrixT.indices.data() + constMatrixT.rowStart[col];
		size_t size = constMatrixT.rowStart[col+1] - constMatrixT.rowStart[col];
		return VectorView<REAL>(coefBegin, indBegin, size);
	}

private:
	friend class mpsreader;

	// min {obj*x}
	std::vector<REAL> objective;

	// lhs <= Ax <= rhs
	std::vector<REAL> lhs;
	std::vector<REAL> rhs;

	// lb <= x <= ub
	std::vector<REAL> lb;
	std::vector<REAL> ub;

	// row-major sparse
	SparseMatrix<REAL> constMatrix;

	// column-major sparse
	SparseMatrix<REAL> constMatrixT;

	// TODO: add column flags
	bitset integer;
};

#endif
