import breeze.linalg._
import scala.io._
import scala.Console._
import breeze.numerics._
import scala.util.control._
import breeze.math._
import scala.util.control.Breaks._
import Helper._
object matrixPow {

  implicit def enrichDenseMatrix(M: DenseMatrix[Complex]) =
    new {
      def sqrtFrac(pow: Double = 3.43): DenseMatrix[Complex] = fract(pow, M)
    }
  /* no checks... */
  def sqrtTriangular(m_A: DenseMatrix[Complex]) = {

    val result = DenseMatrix.zeros[Complex](m_A.cols, m_A.rows)
    var i = 0
    for (i <- 0 to m_A.cols - 1)
      result(i, i) = breeze.numerics.pow(m_A(i, i), 0.5)
    var j = 1
    for (j <- 1 to m_A.cols - 1) {
      for (i <- (j - 1) to 0 by -1) {
        val tmp = result(i, (i + 1) to (j - i)) * result((i + 1) to (j - i), j)
        result(i, j) = (m_A(i, j) - tmp) / (result(i, i) + result(j, j))
      }
    }
    result
  }
  val maxNormForPade: Double = 2.789358995219730e-1; // double precision
  var numberOfSquareRoots = 0

  def getIMinusT(T: DenseMatrix[Complex], numSquareRoots: Int = 0, deg1: Double = 10.0, deg2: Double = 0.0): DenseMatrix[Complex] = {
    val IminusT = DenseMatrix.eye[Complex](T.rows) - T
    val normIminusT = norm1(sum(IminusT(::, *)).t.reduceLeft((x, y) => if (norm1(x) > norm1(y)) x else y))
    numberOfSquareRoots = numSquareRoots
    debugPrint(s"normIminusT:\n$normIminusT\nIminusT:\n$IminusT\n", "", 2)
    if (normIminusT < maxNormForPade) {
      val rdeg1 = padePower.degree(normIminusT)
      val rdeg2 = padePower.degree(normIminusT / 2)
      if (rdeg1 - rdeg2 <= 1.0) return T

      return getIMinusT(upperTriangular(sqrtTriangular(T)), numSquareRoots + 1, rdeg1, rdeg2)
    }
    return getIMinusT(upperTriangular(sqrtTriangular(T)), numSquareRoots + 1, deg1, deg2)

  }
  def fract(pow: Double = 3.43, M: DenseMatrix[Complex] = DenseMatrix((1, 2, 4, 4), (5, 6, 7, 9), (9, 10, 11, 12), (13, 14, 15, 16)).mapValues(Complex(_, 0.0))): DenseMatrix[Complex] =
    {
      debugPrint(s"M:\n$M\n", "", 1)
      val T = upperTriangular(M)
      debugPrint(s"T:\n$T", "", 2)
      val mySchur = complexSchur(M) //side effects hmmm
      val m_T = mySchur.matT
      val m_U = mySchur.matQ

      val MatrixH = mySchur.hess.MatrixH
      val MatrixP = mySchur.hess.MatrixP


      debugPrint(s"MatrixH:\n$MatrixH\n", "", 1)
      debugPrint(s"MatrixQ:\n$MatrixP\n", "", 1)
      debugPrint(s"M:\n$M\n", "", 1)
      debugPrint(mySchur.hess.hCoeffs, "mySchur.hess.hCoeffs", 1)

      debugPrint(s"m_T:\n$m_T\n", "", 1)
      debugPrint(s"m_U:\n$m_U\n", "", 1)
      debugPrint(s"M:\n$M\n", "", 2)
      val frac_power = pow % 1.0

      debugPrint(s"T:\n$T", "", 1)

      numberOfSquareRoots = 0

      val IT = getIMinusT(T)
      debugPrint(numberOfSquareRoots, "number of Square Roots", 1)
      debugPrint(s"IT:\n$IT", "", 1)
      val r = padePower(IT, frac_power)

      debugPrint(s"M:\n$M\n", "", 1)
      r
      /*
	     float frac_power = fmod(power, 1.0f);
	 pade(degree, IminusT, frac_power, res);
	 for (; numberOfSquareRoots; --numberOfSquareRoots) {
	 compute2x2(res, std::ldexp(frac_power, -numberOfSquareRoots), m_T);
	 res = res.triangularView<Upper>();
	 }
	 compute2x2(res, frac_power, m_T);
	 revertSchur(m_tmp, res, m_U);
	 int ipower = floor(power);
	 computeIntPower(res, ipower, M);
	 res = m_tmp * res;

	 */

    }
}

