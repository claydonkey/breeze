import breeze.linalg._
import scala.io._
import scala.Console._
import breeze.numerics._
import scala.util.control._
import breeze.math._
import scala.util.control.Breaks._
import Helper._
import Schur._

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
    //debugPrint(normIminusT, "normIminusT", 2)
    //debugPrint(IminusT, "IminusT", 2)
    if (normIminusT < maxNormForPade) {
      val rdeg1 = padePower.degree(normIminusT)
      val rdeg2 = padePower.degree(normIminusT / 2)
      if (rdeg1 - rdeg2 <= 1.0) return IminusT

      return getIMinusT(upperTriangular(sqrtTriangular(T)), numSquareRoots + 1, rdeg1, rdeg2)
    }
    return getIMinusT(upperTriangular(sqrtTriangular(T)), numSquareRoots + 1, deg1, deg2)

  }
  // atan(z) = (i/2) log((i + z)/(i - z))
  def atan2h(x: Complex, y: Complex): Complex = {
    val z = x / y
    if ((y == 0) || (abs2(z) > pow(GlobalConsts.EPSILON, 0.5)))
      (0.5) * log((y + x) / (y - x))
    else
      z + z * z * z / 3
  }
  val M_PI = 3.14159265358979323846
  def computeSuperDiag(curr: Complex, prev: Complex, p: Double) = {

    import breeze.numerics._
    val logCurr: Complex = log(curr)
    val logPrev: Complex = log(prev)
    val unwindingNumber = ceil(((logCurr - logPrev).imag - M_PI) / (2 * M_PI))
    val w = atan2h(curr - prev, curr + prev) + Complex(0.0, M_PI * unwindingNumber)
    (2.0 * exp(Complex(0.5, 0.0) * p * (logCurr + logPrev)) * scala.math.sinh(p * w.real) / (curr - prev))

  }

  def compute2x2(res: DenseMatrix[Complex], p: Double, m_A: DenseMatrix[Complex]) =
    {
      debugPrint(res, "  Result zero", 2)
      debugPrint(p, "  Result zero", 2)
      debugPrint(m_A, "  Result zero", 2)
      res(0, 0) = pow(m_A(0, 0), p)
      debugPrint(res(0, 0), "  Result zero", 2)
      var i = 0
      for (i <- 1 to m_A.cols - 1) {
        res(i, i) = pow(m_A(i, i), p);
        debugPrint(res(i, i), "  Result zero loop", 2)
        if (m_A(i - 1, i - 1) == m_A(i, i)) {
          res(i - 1, i) = p * pow(m_A(i, i), p - 1)
          debugPrint(res(i - 1, i), "  Result 1", 2)
        } else if (2 * abs(m_A(i - 1, i - 1)) < abs(m_A(i, i)) || 2 * abs(m_A(i, i)) < abs(m_A(i - 1, i - 1))) {
          res(i - 1, i) = (res(i, i) - res(i - 1, i - 1)) / (m_A(i, i) - m_A(i - 1, i - 1))
          debugPrint(res(i - 1, i), "  Result 2", 2)
        } else {
          res(i - 1, i) = computeSuperDiag(m_A(i, i), m_A(i - 1, i - 1), p)
          debugPrint(res(i - 1, i), "  Result 3", 2)
        }
        res(i - 1, i) *= m_A(i - 1, i)
      }
      res
    }

  def computeIntPower(p: Double, in: DenseMatrix[Double]) = {
    var pp = abs(p);
    val m_tmp = if (p < 0.0) { inv(in) } else { in }

    var res = DenseMatrix.eye[Double](in.cols)
    while (pp >= 1) {
      if (pp % 2 >= 1)
        res = m_tmp * res;
      m_tmp *= m_tmp;
      pp = pp / 2
    }
    res
  }

  def fract(pow: Double = 3.43, M: DenseMatrix[Complex] = DenseMatrix((1, 2, 4, 4), (5, 6, 7, 9), (9, 10, 11, 12), (13, 14, 15, 16)).mapValues(Complex(_, 0.0))): DenseMatrix[Complex] =
    {
      debugPrint(M, "M", 2)
      val T = upperTriangular(M)

      val mySchur = complexSchur(M.copy) //side effects hmmm
      val m_T = mySchur.matT
      val m_U = mySchur.matQ

      val MatrixH = mySchur.hess.MatrixH
      val MatrixP = mySchur.hess.MatrixP

      debugPrint(m_T, "m_T", 2)
      debugPrint(m_U, "m_U", 2)
      debugPrint(MatrixH, "MatrixH", 2)
      debugPrint(MatrixP, "MatrixP", 2)

      debugPrint(mySchur.hess.hCoeffs, "mySchur.hess.hCoeffs", 2)

      debugPrint(M, "M", 2)
      debugPrint(T, "T", 2)

      numberOfSquareRoots = 0
      val frac_power = pow % 1.0
      val IT = getIMinusT(T)
      debugPrint(numberOfSquareRoots, "number of Square Roots", 2)
      debugPrint(IT, "IT", 2)
      debugPrint(frac_power, "frac_power", 2)

      var res = padePower(IT, frac_power)

      debugPrint(res, "Starting Result", 2)

      var i = numberOfSquareRoots

      for (i <- numberOfSquareRoots to 0 by -1) {
        res = upperTriangular(compute2x2(res, frac_power * scala.math.pow(2, -i), m_T))
        debugPrint(res, "  Result zero", 2)
      }
      debugPrint(res, "res", 2)

      val m_tmp =(res.revertSchur(m_U.mapValues(_.real))).mapValues(Complex(_, 0.0))
      debugPrint(m_tmp, "m_tmp", 2)
      val ipower = floor(frac_power)
      res = computeIntPower(ipower, M.mapValues(_.real)).mapValues(Complex(_, 0.0))
      res = m_tmp * res
      debugPrint(res, "RESULT", 2)
      res

    }

}