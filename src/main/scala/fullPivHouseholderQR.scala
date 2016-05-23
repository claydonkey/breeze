import breeze.linalg._
import scala.io._
import scala.Console._
import breeze.numerics._
import scala.util.control._
import breeze.math._
import scala.util.control.Breaks._
import Helper._
import GlobalConsts._

class fullPivHouseholderQR(val QR :DenseMatrix[Complex]) {
  val  rows = QR.rows
  val  cols = QR.cols
  val  size = math.min(rows,cols);
  val hCoeffs=   DenseVector.zeros[Complex](size) //tau coeffs
  val rows_transpositions = Array[Int](size)
  val cols_transpositions  = Array[Int](size)

  def compute()={
    breakable
    {
      var k =0
      for ( k  <-  0 to  size- 1)
      {
	var row_of_biggest_in_corner , col_of_biggest_in_corner = 0
	val biggest_in_corner : Double= Helper.biggest(QR(QR.rows -1 - k to QR.rows -1, QR.cols -1- k  to QR.cols -1))
	val m_precision = EPSILON
	row_of_biggest_in_corner += k
	col_of_biggest_in_corner += k
	  var biggest  = 0.0

	if(k==0)
	    biggest = biggest_in_corner

	// if the corner is negligible, then we have less than full rank, and we can finish early
	if(isMuchSmallerThan(biggest_in_corner, biggest))
	{
	  var  m_nonzero_pivots = k
	  var i = 0
	  for(  i <- k to  size -1)
	  {
	    rows_transpositions(i) = i
	    cols_transpositions(i) = i
	    hCoeffs(i) =Complex(0.0,0.0)
	  }
	  break
	}
/*
	rows_transpositions(k) = row_of_biggest_in_corner
	cols_transpositions(k) = col_of_biggest_in_corner
	if(k != row_of_biggest_in_corner) {
	  val t = QR(::,row_of_biggest_in_corner).(cols-k)
	  QR(::,k)(cols-k)swap(QR(::,row_of_biggest_in_corner)(cols-k)
	   number_of_transpositions = number_of_transpositions + 1
	}
	if(k != col_of_biggest_in_corner) {
	  QR.col(k).swap(m_qr.col(col_of_biggest_in_corner));
	  number_of_transpositions = number_of_transpositions +1
	}

	val beta =0.0
	QR.col(k).tail(rows-k).makeHouseholderInPlace(hCoeffs(k), beta);
	QR.(k,k) = beta;

	// remember the maximum absolute value of diagonal coefficients
	if(abs(beta) > m_maxpivot) m_maxpivot = abs(beta);

	QR.bottomRightCorner(rows-k, cols-k-1)
	.applyHouseholderOnTheLeft(m_qr.col(k).tail(rows-k-1), m_hCoeffs.coeffRef(k), &m_temp.coeffRef(k+1));

      }

      m_cols_permutation.setIdentity(cols);
      for(Index k = 0; k < size; ++k)
	m_cols_permutation.applyTranspositionOnTheRight(k, m_cols_transpositions.coeff(k));

      m_det_pq = (number_of_transpositions%2) ? -1 : 1;
      m_isInitialized = true;
*/

    }

  }
}
object fullPivHouseholderQR {

  def apply(M :DenseMatrix[Complex]) =
  {
    new fullPivHouseholderQR(M.copy).compute()
  }
  def   solve(Rhs :  DenseMatrix[Complex]) =
  {

  }
}
}
