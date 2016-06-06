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

/*
 Copyright 2016 Anthony Campbelll

 Licensed under the Apache License, Version 2.0 (the "License")
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
*/

/*
 * Based on:
 *A SCHUR–PADÉ ALGORITHM FOR
* FRACTIONAL POWERS OF A MATRIX*
* NICHOLAS J. HIGHAM† AND LIJING LIN†

 * The algorithm starts with a Schur decomposition, takes kk square roots of the triangular factor TT,  * evaluates
 *  an [m/mm/m] Padé approximant of (1−x)p(1-x)p at I−T1/2kI-T1/2k, and squares the result
 * kk times. The parameters kk and mm are chosen to minimize the cost subject to achieving double precision
 * accuracy in the evaluation of the Padé approximant, making use of a result that bounds the error in the
 *  matrix Padé approximant by the error in the scalar Padé approximant with argument the norm of the matrix.
 * The Padé approximant is evaluated from the continued fraction representation in bottom-up fashion,
 *  which is shown to be numerically stable. In the squaring phase the diagonal and first superdiagonal are
 *  computed from explicit formulae for Tp/2jTp/2j, yielding increased accuracy. Since the basic algorithm is
 * designed for p∈(−1,1)p∈(-1,1), a criterion for reducing an arbitrary real pp to this range is developed,
 * making use of bounds for the condition number of the ApAp problem. How best to compute AkAk for a
 * negative integer kk is also investigated. In numerical experiments the new algorithm is found to be
 * superior in accuracy and stability to several alternatives, including the use of an eigendecomposition
 * and approaches based on the formula Ap=exp(plog(A))Ap=exp(plog(A)).

Read More: http://epubs.siam.org/doi/abs/10.1137/10081232X?journalCode=sjmael
 */

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

  def getIMinusT(T: DenseMatrix[Complex], numSquareRoots: Int = 0, deg1: Double = 10.0, deg2: Double = 0.0): (DenseMatrix[Complex], Int) = {

    val IminusT = DenseMatrix.eye[Complex](T.rows) - T
    val normIminusT = norm1(sum(IminusT(::, *)).t.reduceLeft((x, y) => if (norm1(x) > norm1(y)) x else y))

    if (normIminusT < maxNormForPade) {
      val rdeg1 = degree(normIminusT)
      val rdeg2 = degree(normIminusT / 2)
      if (rdeg1 - rdeg2 <= 1.0) return (IminusT, numSquareRoots) else getIMinusT(upperTriangular(sqrtTriangular(T)), numSquareRoots + 1, rdeg1, rdeg2)
    }
    getIMinusT(upperTriangular(sqrtTriangular(T)), numSquareRoots + 1, deg1, deg2)
  }

  def computeSuperDiag(curr: Complex, prev: Complex, p: Double): Complex = {

    val logCurr: Complex = log(curr)
    val logPrev: Complex = log(prev)
    val unwindingNumber = ceil(((logCurr - logPrev).imag - M_PI) / (2 * M_PI))
    val w = atan2h(curr - prev, curr + prev) + Complex(0.0, M_PI * unwindingNumber)
    val res = (2.0 * exp(0.5 * p * (logCurr + logPrev)) * sinh2(p * w) / (curr - prev))
    res
  }

  def compute2x2(r: DenseMatrix[Complex], m_A: DenseMatrix[Complex], p: Double) = {

    val res = r
    res(0, 0) = pow(m_A(0, 0), p)
    var i = 1
    for (i <- 1 until m_A.cols) {
      res(i, i) = pow(m_A(i, i), p)

      if (m_A(i - 1, i - 1) == m_A(i, i)) {
        res(i - 1, i) = p * pow(m_A(i, i), p - 1)
      } else if (2 * abs(m_A(i - 1, i - 1)) < abs(m_A(i, i)) || 2 * abs(m_A(i, i)) < abs(m_A(i - 1, i - 1))) {
        res(i - 1, i) = (res(i, i) - res(i - 1, i - 1)) / (m_A(i, i) - m_A(i - 1, i - 1))
      } else {
        res(i - 1, i) = computeSuperDiag(m_A(i, i), m_A(i - 1, i - 1), p)
      }
      res(i - 1, i) *= m_A(i - 1, i)
    }
    res
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
      val fpow = pow % 1

      if (abs(fpow) > 0) {
        val (iminusT, noOfSqRts) = getIMinusT(upperTriangular(M))
        val pP = padePower(iminusT, fpow)
        val pT = computeFracPower(pP, sT, fpow, noOfSqRts)
        (Some(pT), sQ)
      } else {
        (None, sQ)
      }
    }

  def dFracPart(MD: DenseMatrix[Double], pow: Double): (Option[DenseMatrix[Complex]], DenseMatrix[Complex]) =
    {
      val (sT, sQ, tau, house) = MD schurDecomp
      val M = MD.mapValues(Complex(_, 0.0))
      val fpow = pow % 1

      if (abs(fpow) > 0) {
        val (iminusT, noOfSqRts) = getIMinusT(upperTriangular(M))
        val pP = padePower(iminusT, fpow)
        val pT = computeFracPower(pP, sT, fpow, noOfSqRts)

        (Some(pT), sQ)
      } else
        (None, sQ)

    }

  def fractC(pow: Double, M: DenseMatrix[Complex]): DenseMatrix[Complex] = {
    val ipower = abs(floor(pow))
    val intPow = if (pow < 0.0)
      throw new IllegalArgumentException("Cannot currently invert complex matrices.")

    return cFracPart(M, pow) match {
      case (None, _) => computeIntPowerC(ipower, M)
      case (pT, sQ) => if (ipower > 0) (pT.get revertSchur sQ) * computeIntPowerC(ipower, M) else (pT.get revertSchur sQ)
    }
  }

  def fractI(pow: Double, M: DenseMatrix[Int]): DenseMatrix[Complex] = fractD(pow, M.mapValues(_.toDouble))
  def fractD(pow: Double, M: DenseMatrix[Double]): DenseMatrix[Complex] = {

    val ipower = abs(floor(pow))
    val intPow = if (pow < 0.0) inv(M) else M
    return dFracPart(M, pow) match {
      case (None, _) => computeIntPowerD(ipower, intPow).mapValues(Complex(_, 0.0))
      case (pT, sQ) => if (ipower > 0) (pT.get revertSchur sQ) * computeIntPowerD(ipower, intPow).mapValues(Complex(_, 0.0)) else (pT.get revertSchur sQ)
    }
  }
}