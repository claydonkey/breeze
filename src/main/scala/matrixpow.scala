import breeze.linalg._
import scala.io._
import scala.Console._
import breeze.numerics._
import scala.util.control._
import breeze.math._
import scala.util.control.Breaks._
import Helper._
import Schur._
import scala.annotation.tailrec

object matrixPow {

  implicit def IMPL_mpowD(M: DenseMatrix[Double]) =
    new {
      def mPow(pow: Double) = fractD(pow, M)
    }

  implicit def IMPL_mpowI(M: DenseMatrix[Int]) =
    new {
      def mPow(pow: Double) = fractI(pow, M)
    }

  implicit def IMPL_mpowC(M: DenseMatrix[Complex]) =
    new {
      def mPow(pow: Double) = fractC(pow, M)
    }

  def degree = (normIminusT: Double) => {
    val maxNormForPade = Array(2.8064004e-1f /* degree = 3 */ , 4.3386528e-1f)
    var degree = 3
    for (degree <- 3 to 4)
      if (normIminusT <= maxNormForPade(degree - 3))
        degree
    degree
  }

  def padePower(IminusT: DenseMatrix[Complex], m_p: Double) = {

    val _degree = degree(m_p)
    val i = _degree << 1
    val res = IminusT.map(_ * (m_p - _degree.toDouble) / ((i - 1) << 1))
    val index = 0
    val M: DenseMatrix[Complex] = DenseMatrix.tabulate[Complex](res.rows, res.cols) { (x, y) => if (x == y) Complex(1.0, 0) else res(x, y) }
    val T1 = -1.0 * Complex(m_p, 0)
    val T = IminusT * T1
    (M.mapValues(_.real) \ T.mapValues(_.real)).mapValues(Complex(_, 0.0)) :+ DenseMatrix.eye[Complex](IminusT.rows) // BIG PROBLEMMMO

  }

  /* no checks... */
  def sqrtTriangular(m_A: DenseMatrix[Complex]) = {

    val result = DenseMatrix.zeros[Complex](m_A.cols, m_A.rows)
    var i = 0
    for (i <- 0 to m_A.cols - 1)
      result(i, i) = breeze.numerics.pow(m_A(i, i), 0.5)
    var j = 1
    for (j <- 1 until m_A.cols) {
      for (i <- (j - 1) to 0 by -1) {
        val tmp = result(i, (i + 1) to (j - i)) * result((i + 1) to (j - i), j)
        result(i, j) = (m_A(i, j) - tmp) / (result(i, i) + result(j, j))
      }
    }
    result
  }

  val maxNormForPade: Double = 2.789358995219730e-1;

  // double precision

  def getIMinusT(T: DenseMatrix[Complex], numSquareRoots: Int = 0, deg1: Double = 10.0, deg2: Double = 0.0): (DenseMatrix[Complex], Int) = {
    val IminusT = DenseMatrix.eye[Complex](T.rows) - T
    val normIminusT = norm1(sum(IminusT(::, *)).t.reduceLeft((x, y) => if (norm1(x) > norm1(y)) x else y))
    debugPrint(IminusT, "IminusT", 5)
    debugPrint(normIminusT, "normIminusT", 5)
    debugPrint(numSquareRoots, "numSquareRoots", 5)
    if (normIminusT < maxNormForPade) {
      val rdeg1 = degree(normIminusT)
      val rdeg2 = degree(normIminusT / 2)
      if (rdeg1 - rdeg2 <= 1.0) return (IminusT, numSquareRoots) else getIMinusT(upperTriangular(sqrtTriangular(T)), numSquareRoots + 1, rdeg1, rdeg2)
    }
    getIMinusT(upperTriangular(sqrtTriangular(T)), numSquareRoots + 1, deg1, deg2)
  }

  def computeSuperDiag(curr: Complex, prev: Complex, p: Double): Complex = {

    import breeze.numerics._
    val logCurr: Complex = log(curr)
    debugPrint(logCurr, "  logCurr", 3)
    val logPrev: Complex = log(prev)
    debugPrint(logPrev, "  logPrev", 3)
    val unwindingNumber = ceil(((logCurr - logPrev).imag - M_PI) / (2 * M_PI))
    debugPrint(unwindingNumber, "  unwindingNumber", 3)
    val w = atan2h(curr - prev, curr + prev) + Complex(0.0, M_PI * unwindingNumber)
    debugPrint(w, "  w", 3)
    val res = (2.0 * exp(0.5 * p * (logCurr + logPrev)) * sinh2(p * w) / (curr - prev))
    debugPrint(res, "  res", 3)
    res
  }

  def compute2x2(r: DenseMatrix[Complex], m_A: DenseMatrix[Complex], p: Double) = {
    import breeze.numerics._
    val res = r

    debugPrint(p, "  Result p", 3)
    debugPrint(m_A, "  Result m_A", 3)
    res(0, 0) = pow(m_A(0, 0), p)
    debugPrint(res(0, 0), "  Result res(0,0)", 3)
    var i = 1
    for (i <- 1 until m_A.cols) {
      debugPrint(p, "  p", 3)
      debugPrint(res(i, i), "  res(i, i) ", 3)
      debugPrint(m_A(i, i), "  m_A.coeff(i, i)", 3)

      res(i, i) = pow(m_A(i, i), p)

      debugPrint(res(i, i), "  Result  loop 1", 3)
      if (m_A(i - 1, i - 1) == m_A(i, i)) {
        res(i - 1, i) = p * pow(m_A(i, i), p - 1)
        debugPrint(res(i - 1, i), "  Result 1", 3)
      } else if (2 * abs(m_A(i - 1, i - 1)) < abs(m_A(i, i)) || 2 * abs(m_A(i, i)) < abs(m_A(i - 1, i - 1))) {

        debugPrint(res(i, i), "  res(i , i)", 3)
        debugPrint(res(i - 1, i), "  res(i - 1, i)", 3)
        debugPrint(res(i - 1, i - 1), "  res(i - 1, i-1)", 3)
        debugPrint(m_A(i, i), "  m_A(i, i)", 3)
        debugPrint(m_A(i - 1, i - 1), "  m_A(i-1, i-1)", 3)
        res(i - 1, i) = (res(i, i) - res(i - 1, i - 1)) / (m_A(i, i) - m_A(i - 1, i - 1))

        debugPrint(res(i - 1, i), "  Result 2", 3)
      } else {
        res(i - 1, i) = computeSuperDiag(m_A(i, i), m_A(i - 1, i - 1), p)
        debugPrint(res(i - 1, i), "  Result 3", 3)
      }
      res(i - 1, i) *= m_A(i - 1, i)
    }
    res
  }

  def compute2x2b(pade: DenseMatrix[Complex], sT: DenseMatrix[Complex], p: Double) = {
    import breeze.numerics._

    val R = pade
    R(0, 0) = pow(sT(0, 0), p)
    var i = 1
    for (i <- 1 until sT.cols) {
      R(i, i) = pow(sT(i, i), p)

      if (sT(i - 1, i - 1) == sT(i, i)) {
        R(i - 1, i) = p * pow(sT(i, i), p - 1)

      } else if (2 * abs(sT(i - 1, i - 1)) < abs(sT(i, i)) || 2 * abs(sT(i, i)) < abs(sT(i - 1, i - 1))) {
        R(i - 1, i) = (pade(i, i) - pade(i - 1, i - 1)) / (sT(i, i) - sT(i - 1, i - 1))

      } else {
        R(i - 1, i) = computeSuperDiag(sT(i, i), sT(i - 1, i - 1), p)
      }
      R(i - 1, i) *= sT(i - 1, i)
    }
    R
  }

  def compute2x2c(pade: DenseMatrix[Complex], sT: DenseMatrix[Complex], p: Double) = {
    import breeze.numerics._

    val R = pade
    R(0, 0) = pow(sT(0, 0), p)
    var i = 1

    for (i <- 1 until sT.cols) {

      val prev = sT(i - 1, i - 1)
      val curr = sT(i, i)
      val pPade = pade(i - 1, i - 1)
      val cPade = pade(i, i)

      R(i, i) = pow(curr, p)

      R(i - 1, i) = if (prev == curr) {
        p * pow(curr, p - 1)
      } else if (2 * abs(prev) < abs(curr) || 2 * abs(curr) < abs(prev)) {
        (cPade - pPade) / (curr - prev)
      } else {
        computeSuperDiag(sT(i, i), sT(i - 1, i - 1), p)
      } * sT(i - 1, i)

    }
    R
  }

  def computeIntPowerD(exp: Double, value: DenseMatrix[Double]): DenseMatrix[Double] = {
    if (exp == 1) value
    else if (exp % 2 == 1) value * (computeIntPowerD(exp - 1, value))
    else {
      val half = computeIntPowerD(exp / 2, value)
      half * half
    }
  }

  def computeIntPowerC(exp: Double, value: DenseMatrix[Complex]): DenseMatrix[Complex] = {
    if (exp == 1) value
    else if (exp % 2 == 1) value * (computeIntPowerC(exp - 1, value))
    else {
      val half = computeIntPowerC(exp / 2, value)
      half * half
    }
  }

  @tailrec
  def computeFracPower(value: DenseMatrix[Complex], sT: DenseMatrix[Complex], frac_power: Double, noOfSqRts: Int): DenseMatrix[Complex] = {
    if (noOfSqRts == 0) return upperTriangular(compute2x2(value, sT, frac_power * scala.math.pow(2, -noOfSqRts)))
    else
      computeFracPower(upperTriangular(compute2x2(value, sT, frac_power * scala.math.pow(2, -noOfSqRts))), sT, frac_power, noOfSqRts - 1)
  }

  def cFracPart(M: DenseMatrix[Complex], pow: Double): (Option[DenseMatrix[Complex]], DenseMatrix[Complex]) =
    {
      val (sT, sQ, tau, house) = M schurDecomp

      debugPrint(sT, "schur T", 1)
      debugPrint(sQ, "schur Q", 1)
      debugPrint(house, "house Q", 1)
      val fpow = pow % 1

      if (abs(fpow) > 0) {
        val (iminusT, noOfSqRts) = getIMinusT(upperTriangular(M))
        val pP = padePower(iminusT, fpow)

        debugPrint(iminusT, "iminusT", 1)
        debugPrint(fpow, "frac_power", 1)
        debugPrint(noOfSqRts, "  noOfSqRts", 1)
        debugPrint(pP, "padePower", 1)

        val pT = computeFracPower(pP, sT, fpow, noOfSqRts)

        debugPrint(sT, "schur T", 1)
        debugPrint(sQ, "schur Q", 1)
        debugPrint(pT, "pade T", 1)
        (Some(pT), sQ)
      } else
        (None, sQ)

    }

  def dFracPart(MD: DenseMatrix[Double], pow: Double): (Option[DenseMatrix[Complex]], DenseMatrix[Complex]) =
    {
      val (sT, sQ, tau, house) = MD schurDecomp
      val M = MD.mapValues(Complex(_, 0.0))
      debugPrint(sT, "schur T", 1)
      debugPrint(sQ, "schur Q", 1)
      debugPrint(house, "house Q", 1)
      val fpow = pow % 1

      if (abs(fpow) > 0) {
        val (iminusT, noOfSqRts) = getIMinusT(upperTriangular(M))
        val pP = padePower(iminusT, fpow)

        debugPrint(iminusT, "iminusT", 1)
        debugPrint(fpow, "frac_power", 1)
        debugPrint(noOfSqRts, "  noOfSqRts", 1)
        debugPrint(pP, "padePower", 1)

        val pT = computeFracPower(pP, sT, fpow, noOfSqRts)

        debugPrint(sT, "schur T", 1)
        debugPrint(sQ, "schur Q", 1)
        debugPrint(pT, "pade T", 1)
        (Some(pT), sQ)
      } else
        (None, sQ)

    }

  def fractC(pow: Double, M: DenseMatrix[Complex]): DenseMatrix[Complex] = {
    val ipower = abs(floor(pow))
    return cFracPart(M, pow) match {
      case (None, _) => computeIntPowerC(ipower, M)
      case (pT, sQ) => if (ipower > 0) (pT.get revertSchur sQ) * computeIntPowerC(ipower, M) else (pT.get revertSchur sQ)
    }
  }

  def fractI(pow: Double, M: DenseMatrix[Int]): DenseMatrix[Complex] = fractD(pow, M.mapValues(_.toDouble))

  def fractD(pow: Double, M: DenseMatrix[Double]): DenseMatrix[Complex] = {



    val ipower = abs(floor(pow))
    val intPow = if (ipower < 0.0) inv(M) else M
    return dFracPart(M, pow) match {
          //return cFracPart(M.mapValues(Complex(_, 0.0)), pow) match {
      case (None, _) => computeIntPowerD(ipower, intPow).mapValues(Complex(_, 0.0))
      case (pT, sQ) => if (ipower > 0) (pT.get revertSchur sQ) * computeIntPowerD(ipower, intPow).mapValues(Complex(_, 0.0)) else (pT.get revertSchur sQ)
    }

  }

}