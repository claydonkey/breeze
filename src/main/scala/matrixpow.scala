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
      def mPow(pow: Double): DenseMatrix[Complex] = fract(pow, M)
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

  def getIMinusT(T: DenseMatrix[Complex], numSquareRoots: Int = 0, deg1: Double = 10.0, deg2: Double = 0.0): (DenseMatrix[Complex], Int) = {
    val IminusT = DenseMatrix.eye[Complex](T.rows) - T
    val normIminusT = norm1(sum(IminusT(::, *)).t.reduceLeft((x, y) => if (norm1(x) > norm1(y)) x else y))
    if (normIminusT < maxNormForPade) {
      val rdeg1 = padePower.degree(normIminusT)
      val rdeg2 = padePower.degree(normIminusT / 2)
      if (rdeg1 - rdeg2 <= 1.0) return (IminusT, numSquareRoots) else getIMinusT(upperTriangular(sqrtTriangular(T)), numSquareRoots + 1, rdeg1, rdeg2)
    }
    getIMinusT(upperTriangular(sqrtTriangular(T)), numSquareRoots + 1, deg1, deg2)
  }
  // atan(z) = (i/2) log((i + z)/(i - z))
  def atan2h(x: Complex, y: Complex): Complex = {
    val z = x / y
    if ((y == 0) || (abs2(z) > pow(GlobalConsts.EPSILON, 0.5)))
      (0.5) * log((y + x) / (y - x))
    else
      z + z * z * z / 3
  }

  def computeSuperDiag(curr: Complex, prev: Complex, p: Double): Complex = {

    import breeze.numerics._
    val logCurr: Complex = log(curr)
    debugPrint(logCurr, "  logCurr", 2)
    val logPrev: Complex = log(prev)
    debugPrint(logPrev, "  logPrev", 2)
    val unwindingNumber = ceil(((logCurr - logPrev).imag - M_PI) / (2 * M_PI))
    debugPrint(unwindingNumber, "  unwindingNumber", 2)
    val w = atan2h(curr - prev, curr + prev) + Complex(0.0, M_PI * unwindingNumber)
    debugPrint(w, "  w", 2)
    val res = (2.0 * exp(0.5 * p * (logCurr + logPrev)) * sinh2(p * w) / (curr - prev))
    debugPrint(res, "  res", 2)
    res
  }

  def compute2x2(res: DenseMatrix[Complex], p: Double, m_A: DenseMatrix[Complex]) =
    {
      debugPrint(res, "  Result res", 2)
      debugPrint(p, "  Result p", 2)
      debugPrint(m_A, "  Result m_A", 2)
      res(0, 0) = breeze.numerics.pow(m_A(0, 0), p)
      debugPrint(res(0, 0), "  Result res(0,0)", 2)
      var i = 1
      for (i <- 1 to m_A.cols - 1) {
        res(i, i) = breeze.numerics.pow(m_A(i, i), p);
        debugPrint(res(i, i), "  Result  loop 1", 2)
        if (m_A(i - 1, i - 1) == m_A(i, i)) {
          res(i - 1, i) = p * breeze.numerics.pow(m_A(i, i), p - 1)
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

  def computeIntPower2(p: Double, in: DenseMatrix[Double]) = {
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

  def computeIntPower(exp: Double, value: DenseMatrix[Double]): DenseMatrix[Double] = {
    if (exp == 1) value
    else if (exp % 2 == 1) value * (pow(exp - 1, value))
    else {
      val half = pow(exp / 2, value)
      half * half
    }
  }

  def fract(pow: Double = 3.43, M: DenseMatrix[Complex] = DenseMatrix((1, 2, 4, 4), (5, 6, 7, 9), (9, 10, 11, 12), (13, 14, 15, 16)).mapValues(Complex(_, 0.0))): DenseMatrix[Complex] =
    {
      debugPrint(M, "M", 1)
      val T = upperTriangular(M)

      val (m_T, m_U, m_P, m_H, hCoeffs) = M getSchur

      debugPrint(m_T, "m_T", 1)
      debugPrint(m_U, "m_U", 1)
      debugPrint(m_H, "MatrixH", 1)
      debugPrint(m_P, "MatrixP", 1)

      debugPrint(hCoeffs, "mySchur.hess.hCoeffs", 1)

      debugPrint(M, "M", 1)
      debugPrint(T, "T", 1)

      val frac_power = pow % 1.0
      val (iminusT: DenseMatrix[Complex], noOfSqRts: Int) = getIMinusT(T)

      debugPrint(noOfSqRts, "  noOfSqRts", 1)
      debugPrint(iminusT, "iminusT", 1)
      debugPrint(frac_power, "frac_power", 1)



      /*
       *      //var res = padePower(iminusT, frac_power)
       *
             def sum2x2(noOfSqRts: Int, value: DenseMatrix[Complex], m_T: DenseMatrix[Complex]): DenseMatrix[Complex] = {
        if (noOfSqRts == 0) upperTriangular(compute2x2(value, frac_power * scala.math.pow(2, -noOfSqRts), m_T))
        else
          sum2x2(noOfSqRts - 1, upperTriangular(compute2x2(value, frac_power * scala.math.pow(2, -noOfSqRts), m_T)), m_T)
      }
    var res  = new Function3[Int, DenseMatrix[Complex], DenseMatrix[Complex], DenseMatrix[Complex]] {
      def apply(noOfSqRts: Int, value: DenseMatrix[Complex], m_T: DenseMatrix[Complex]): DenseMatrix[Complex] =
	if (noOfSqRts == 0) upperTriangular(compute2x2(value, frac_power * scala.math.pow(2, -noOfSqRts), m_T))
      else
	apply(noOfSqRts - 1, upperTriangular(compute2x2(value, frac_power * scala.math.pow(2, -noOfSqRts), m_T)), m_T)
    }.apply(noOfSqRts, padePower(iminusT, frac_power), m_T)
    debugPrint(res, "Starting Result", 2)
*/

      /* if tail recursion is a problem
      var i = 0

      for (i <- noOfSqRts to 0 by -1) {
        res = upperTriangular(compute2x2(res, frac_power * scala.math.pow(2, -i), m_T))
      }
*/

      val res = Y3[Int, DenseMatrix[Complex], DenseMatrix[Complex], DenseMatrix[Complex]](f => (noOfSqRts, value, m_T) =>
        if (noOfSqRts == 0) upperTriangular(compute2x2(value, frac_power * scala.math.pow(2, -noOfSqRts), m_T))
        else
          f(noOfSqRts - 1, upperTriangular(compute2x2(value, frac_power * scala.math.pow(2, -noOfSqRts), m_T)), m_T))(noOfSqRts, padePower(iminusT, frac_power), m_T)

      debugPrint(res, "res", 1)

      val m_tmp = (res.revertSchur(m_U))
      debugPrint(m_tmp, "revertSchur", 1)
      val ipower = floor(frac_power)
      val res2 = m_tmp * computeIntPower(ipower, M.mapValues(_.real)).mapValues(Complex(_, 0.0))

      debugPrint(res2, "RESULT", 1)
      res2

    }

}