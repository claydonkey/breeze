import breeze.linalg._
import scala.io._
import scala.Console._
import breeze.numerics._
import scala.util.control._
import breeze.math._
import scala.util.control.Breaks._

object Main {
  def promptEnterKey(): Option[Unit] = if (Console.in.read > 10) None else promptEnterKey

  def applySchur() = {
    val M = DenseMatrix((1, 2, 4, 4), (5, 6, 7, 9), (9, 10, 11, 12), (13, 14, 15, 16)).mapValues(Complex(_, 0))
    val S = schur(M)
    val matT = S.matT
    val matU = S.matU
    println(s"matU:\n $matU   \n matT: \n$matT\n")
    S.reduceToTriangularForm

  }

  def main(args: Array[String]): Unit = {
    promptEnterKey();
      applySchur();
/*
    val H = householder(DenseMatrix((1, 2, 3, 4), (5, 6, 7, 8), (9, 10, 11, 12), (13, 14, 15, 16)).mapValues(Complex(_, 0)), DenseVector(1, 2, 3).mapValues(Complex(_, 0)))
    val HD = hessenbergDecomposition(H)
    val matH = HD.MatrixQ
    val test = householder.householderSequence(DenseMatrix((1, 2, 3, 4), (5, 6, 7, 8), (9, 10, 11, 12), (13, 14, 15, 16)).mapValues(_.toDouble), DenseVector(1, 2, 3).mapValues(_.toDouble), 1)

    println("test\n" + test.mapValues(x => x - (x % 0.0001)))*/
  }
}

object Helper {
  def abs2(n: Complex) = { (n.real * n.real) + (n.imag * n.imag) }
  def conj(n: Complex) = { Complex(n.real, -n.imag) }
  def norm1(n: Complex): Double = { abs(n.real) + abs(n.imag) }
}

import Helper._

class jacobiRotation(val m_c: Complex, val m_s: Complex) {
  def transpose() = jacobiRotation(m_c, -conj(m_s))
  def adjoint() = jacobiRotation(conj(m_c), -m_s)
  def applyRotationinPlane(_x: DenseVector[Complex], _y: DenseVector[Complex]) = {
    jacobiRotation.applyRotationinPlane(_x, _y, this)
  }
}

object jacobiRotation {
  import Helper._
  def apply(t: (Complex, Complex)) = { new jacobiRotation(t._1, t._2) }
  def apply(m_c: Complex, m_s: Complex) = { new jacobiRotation(m_c, m_s) }

  def applyRotationinPlane(_x: DenseVector[Complex], _y: DenseVector[Complex], j: jacobiRotation) = {

    println(_x.mapValues(x => x.real - (x.real % 0.0001)))
    println(_y.mapValues(x => x.real - (x.real % 0.0001)))
    assert(_x.length == _y.length)
    if (j.m_c == 1 && j.m_s == 0)
      DenseMatrix.zeros[Complex](_x.size, 2)
    DenseMatrix.vertcat((DenseVector.tabulate(_x.size)(i => j.m_c * _x(i) + conj(j.m_s) * _y(i))).toDenseMatrix, (DenseVector.tabulate(_x.size)(i => -j.m_s * _x(i) + conj(j.m_c) * _y(i))).toDenseMatrix)
  }

  /*This function implements the continuous Givens rotation
   *generation algorithm found in Anderson (2000),
   *Discontinuous Plane Rotations and the Symmetric Eigenvalue Problem.
   *LAPACK Working Note 150, University of Tennessee, UT-CS-00-454, December 4, 2000. */
  def makeGivens(p: Complex, q: Complex) = {

    (p, q) match {
      case (Complex(0.0, 0.0), _) =>
        val m_c = if (p.real < 0) Complex(-1, 0) else Complex(1, 0)
        val m_s = Complex(0, 0)
        //if(r) *r = m_c * p;
        new jacobiRotation(m_c, m_s)
      case (_, Complex(0.0, 0)) =>
        val m_c = Complex(0, 0)
        val m_s = -q / abs(q)
        new jacobiRotation(m_c, m_s)
      case _ =>
        var p1 = norm1(p)
        var q1 = norm1(q)
        if (p1 >= q1) {
          val ps = p / p1
          val p2 = abs2(ps)
          val qs = q / p1
          val q2 = abs2(qs)

          var u = breeze.numerics.pow(Complex(1.0, 0.0) + (q2 / p2), 0.5)
          if (p.real < 0)
            u = -u;
          val m_c = 1 / u
          val m_s = -qs * conj(ps) * (m_c / p2)
          new jacobiRotation(m_c, m_s)
        } else {
          var ps = p / q1
          val p2 = abs2(ps)
          val qs = q / q1
          val q2 = abs2(qs);

          var u = q1 * breeze.numerics.pow((p2 + q2), 0.5)
          if (p.real < 0)
            u = -u
          p1 = abs(p)
          ps = p / p1
          val m_c = p1 / u
          val m_s = -conj(ps) * (q / u)
          new jacobiRotation(Complex(m_c, 0.0), m_s)
        }
    }
  }
}

/*
 *  4 x 4 example
 *    x    x    x    x
 *    x    [ Row0 ]
 *    e    [  bott   ]
 *    e    [   bott  ]
 */
class householder(val matrixH: DenseMatrix[Complex], val tau: DenseVector[Complex]) {

  def this(matrixH: DenseMatrix[Complex]) = this(matrixH, DenseVector.zeros[Complex](matrixH.cols - 1))
  val size = matrixH.cols - 1
  val beta = Array.ofDim[Double](size)
  val essential = Array.ofDim[DenseVector[Complex]](size)

  def applyHouseholder(cnt: Int) =
    {
      essential(cnt) = matrixH((cnt + 2) to matrixH.rows - 1, cnt)
      val eNorm = if (essential(cnt).length == 0) 0.0 else sum(essential(cnt).map(x => scala.math.pow(x.real, 2) + scala.math.pow(x.imag, 2))) // Does Complex component need squaring?
      val c0 = matrixH(cnt + 1, cnt);
      (eNorm, c0.imag) match {
        case (0, 0) =>
          beta(cnt) = c0.real
          tau(cnt) = Complex(0, 0)
        //   essential(cnt) = DenseVector.zeros[Complex](essential(cnt).length)
        case _ =>
          val c0norm = scala.math.pow(c0.real, 2) + scala.math.pow(c0.imag, 2)
          beta(cnt) = if (c0.real >= 0) -Math.sqrt(c0norm + eNorm) else Math.sqrt(c0norm + eNorm)
          tau(cnt) = ((beta(cnt) - c0) / beta(cnt))
          essential(cnt) = (essential(cnt) / (c0 - beta(cnt)))
      }
      matrixH((cnt + 1), cnt) = Complex(beta(cnt), 0)
      matrixH((cnt + 2) to matrixH.rows - 1, cnt) := essential(cnt)
      val matH2 = matrixH.mapValues(_.real)
      //     println(s"applyHouseholder matrixH $cnt \n  $matH2\n")

      this

    }

  /*
   *  4 x 4 example
   *    x    c0    Right
   *    x    c0    Right
   *    e    c0    Right
   *    e    c0    Right
   */
  def applyHouseholderRight(cnt: Int) = {

    if (matrixH.cols == 1) {
      matrixH *= 1 - tau(cnt)
      this
    } else {
      var c0 = matrixH(::, (cnt + 1) to (cnt + 1))
      var right = matrixH(::, (cnt + 2) to matrixH.cols - 1)
      val tmp = (right * (essential(cnt).toDenseMatrix.map(x => Complex(x.real, -x.imag))).t) + c0
      c0 -= tmp * tau(cnt)
      right -= tmp * tau(cnt) * essential(cnt).toDenseMatrix
      val matH2 = matrixH.mapValues(_.real)
      //    println(s"applyHouseholderRight matrixH $cnt \n $matH2\n")
      this
    }
  }
  /*
   *  4 x 4 example
   *    x    x    x    x
   *    x    r0  r0  r0
   *    e     Bottom
   *    e     Bottom
   */
  def applyHouseholderBottom(cnt: Int): householder = {

    if (matrixH.cols == 1) {
      matrixH *= 1 - tau(cnt)
      this
    } else {
      var r0 = matrixH((cnt + 1), (cnt + 1) to matrixH.cols - 1).t
      var bottom = matrixH((cnt + 2) to matrixH.rows - 1, (cnt + 1) to matrixH.cols - 1)
      val tmp = (bottom.t * essential(cnt).map(x => Complex(x.real, -x.imag))) + r0
      r0 -= (tmp.toDenseMatrix * tau(cnt)).toDenseVector
      bottom -= (essential(cnt).toDenseMatrix.t * tmp.toDenseMatrix) * tau(cnt)
      val matH2 = matrixH.mapValues(_.real)
      //     println(s"applyHouseholderBottom matrixH $cnt \n $matH2\n")
      this
    }
  }

}
/*
 *  4 x 4 example
 *    x    x    x    x
 *    c0  x    x    x
 *    e    x    x    x
 *    e    x    x    x
 */

object householder {

  def apply(m_matrix: DenseMatrix[Complex], tau: DenseVector[Complex]): householder = new householder(m_matrix, tau)
  def apply(m_matrix: DenseMatrix[Complex]): householder = new householder(m_matrix)

  def householderSequence(hMatrix: DenseMatrix[Double], hCoeffs: DenseVector[Double], order: Int) = {
    //HouseholderSequence shifted 1 with size -1
    /*
     *  4 x 4 example
     *    x    x    x     x
     *   H    H    H    x
     *   H    H    H    x
     *   H    H    H    x
     */

    if ((hMatrix.rows - order) != hCoeffs.length)
      throw new MatrixNotSquareException // change to correct exception

    val matHS = hMatrix(order to (hMatrix.cols - 1), 0 to (hMatrix.rows - order - 1))

    val I = DenseMatrix.eye[Double](matHS.cols)
    val hhMv = (upperTriangular(matHS.t) *:* -(I - 1.0)) +:+ I
    var sum = I
    var cnt = 0

    for (cnt <- 0 to matHS.cols - 1) {
      val hHvect = hhMv(cnt, ::).t
      //    val adj = det(hHvect.toDenseMatrix) * inv(hHvect.toDenseMatrix)  //  adjoint of vector is vector itself??
      val adj = hHvect.toDenseMatrix
      val newa = -(hCoeffs(cnt) * hHvect * adj - I)
      sum = sum * newa

    }
    cnt = 0

    println("sum\n" + sum.mapValues(x => x - (x % 0.0001)))
    //  println("sum\n" + sum.mapValues(x => x - (x % 0.0001)))
    if (order != 0) {
      val wrapper = DenseMatrix.zeros[Double](hMatrix.cols, hMatrix.cols)
      println("wrapper\n" + wrapper.mapValues(x => x - (x % 0.0001)))
      wrapper(order to wrapper.cols - 1, order to wrapper.cols - 1) := sum
      for (cnt <- 0 to order - 1) {
        wrapper(cnt, cnt) = 1.0;

      }
      wrapper
    } else
      sum
  }
}

class schur(val hess: hessenbergDecomposition) {

  val matT = hess.MatrixH()
  val matU = hess.MatrixQ()
  val m_maxIterationsPerRow = 30
  val maxIters = m_maxIterationsPerRow * hess.matH.rows

  /**
   * If matT(i+1,i) is neglegible in floating point arithmetic
   * compared to matT(i,i) and matT(j,j), then set it to zero and
   * return true, else return false.
   */
  val EPSILON: Double = 0.00001

  def isMuchSmallerThan(x: Double, y: Double) = { abs(x) <= abs(y) * EPSILON }

  def subdiagonalEntryIsNeglegible(i: Int) =
    {
      val d = norm1(matT(i, i)) + norm1(matT(i + 1, i + 1))
      val sd = norm1(matT(i + 1, i))
      if (isMuchSmallerThan(sd, d)) {
        matT(i + 1, i) = Complex(0, 0)
        true
      }
      false
    }

  /** Compute the shift in the current QR iteration. */
  def computeShift(iu: Int, iter: Int) = {
    if (iter == 10 || iter == 20) {
      // exceptional shift, taken from http://www.netlib.org/eispack/comqr.f
      abs(matT(iu, iu - 1).real) + abs(matT(iu - 1, iu - 2).real)
    }
    // compute the shift as one of the eigenvalues of t, the 2x2
    // diagonal block on the bottom of the active submatrix
    val NormT = sum(matT((iu - 1) to matT.cols - 1, (iu - 1) to matT.cols - 1).mapValues(abs(_)))
    val T = matT((iu - 1) to matT.cols - 1, (iu - 1) to matT.cols - 1) / Complex(NormT, 0)
    /*
       (-0.464689,0)  (0.421814,0)
       (-0.120907,0)  (0.141053,0)
       T
       (-0.404618,0)  (0.367286,0)
       (-0.105277,0)  (0.122819,0)
       */

    //  println(s"T\n $T\n")
    val b = T(0, 1) * T(1, 0)
    val c = T(0, 0) - T(1, 1)
    val disc = breeze.numerics.pow((c * c + 4.0 * b), 0.5)
    val det = T(0, 0) * T(1, 1) - b
    val trace = T(0, 0) + T(1, 1)
    var eival1 = (trace + disc) / 2.0
    var eival2 = (trace - disc) / 2.0;

    if (norm1(eival1) > norm1(eival2))
      eival2 = det / eival1;
    else
      eival1 = det / eival2;
    //  println(s"eival  $eival1 \n $eival2 \n")
    // choose the eigenvalue closest to the bottom entry of the diagonal
    if (norm1(eival1 - T(1, 1)) < norm1(eival2 - T(1, 1)))
      NormT * eival1;
    else
      NormT * eival2;

  }

  def reduceToTriangularForm() = {
    /**A horrbly disFunctional implementation using breaks, iterations  and mutations BLEUURRGGH will change... **/
    // The matrix matT is divided in three parts.
    // Rows 0,...,il-1 are decoupled from the rest because matT(il,il-1) is zero.
    // Rows il,...,iu is the part we are working on (the active submatrix).
    // Rows iu+1,...,end are already brought in triangular form.
    var iu = matT.cols - 1
    var il: Int = 0
    var iter = 0 // number of iterations we are working on the (iu,iu) element
    var totalIter = 0 // number of iterations for whole matrix
    breakable {
      while (true) {
        // find iu, the bottom row of the active submatrix
        breakable {
          while (iu > 0) {
            if (!subdiagonalEntryIsNeglegible(iu - 1)) break
            iter = 0
            iu = iu - 1
          }
        }

        // if iu is zero then we are done; the whole matrix is triangularized
        if (iu == 0) break

        // if we spent too many iterations, we give up
        iter = iter + 1
        totalIter = totalIter + 1
        if (totalIter > maxIters) break

        // find il, the top row of the active submatrix
        il = iu - 1
        while (il > 0 && !subdiagonalEntryIsNeglegible(il - 1)) {
          il = il - 1
        }
        /* perform the QR step using Givens rotations. The first rotation
	   *creates a bulge; the (il+2,il) element becomes nonzero. This
	   *bulge is chased down to the bottom of the active submatrix.
	   */
        val shift = computeShift(iu, iter)
        val rot = jacobiRotation.makeGivens(matT(il, il) - shift, matT(il + 1, il))

        val m_ca = rot.adjoint().m_c
        val m_sa = rot.adjoint().m_s
        val m_c = rot.m_c
        val m_s = rot.m_s

        matT(0 to 1, il to matT.cols - 1) := jacobiRotation.applyRotationinPlane(matT(il, ::).t, matT(il + 1, ::).t, rot.adjoint)
        matT(il to matT.rows - 2, 0 to 1) := jacobiRotation.applyRotationinPlane(matT(0 to matT.cols - 2 - il, il), matT(0 to matT.cols - 2 - il, il + 1), rot.transpose).t

        //       println(matT.mapValues(x => x.real - (x.real % 0.0001)))

        val matUC = matU.mapValues(Complex(_, 0.0))
        val temp = jacobiRotation.applyRotationinPlane(matUC(0 to matUC.cols - 1 - il, il), matUC(0 to matUC.cols - 1 - il, il + 1), rot.transpose).t
        //       println("TEMP\n" + temp.mapValues(x => x.real - (x.real % 0.0001)))
        // m_matT.topRows((std::min)(il + 2, iu) + 1).applyOnTheRight(il, il + 1, rot);
        /*
	   if (computeU) m_matU.applyOnTheRight(il, il + 1, rot);

	   for (Index i = il + 1; i < iu; i++) {
	   rot.makeGivens(m_matT.coeffRef(i, i - 1), m_matT.coeffRef(i + 1, i - 1), &m_matT.coeffRef(i, i - 1));
	   m_matT.coeffRef(i + 1, i - 1) = ComplexScalar(0);
	   m_matT.rightCols(m_matT.cols() - i).applyOnTheLeft(i, i + 1, rot.adjoint());
	   m_matT.topRows((std::min)(i + 2, iu) + 1).applyOnTheRight(i, i + 1, rot);
	   if (computeU) m_matU.applyOnTheRight(i, i + 1, rot);
	   }*/

      }
    }
}
  }

object schur {

  def apply(M: DenseMatrix[Complex]): schur = {

    val S = new schur(hessenbergDecomposition(M))
    if (M.rows != M.cols)
      throw new MatrixNotSquareException

    if (M.cols == 1) {
      S //   return (DenseMatrix.eye[Complex](1), DenseMatrix.eye[Complex](1))
    }

    S
    //  (M, M)
  }

  /**
   *  Hessenberg decomposition of given matrix.
   *
   * The Hessenberg decomposition is computed by bringing the columns of the
   * matrix successively in the required form using Householder reflections
   * (see, e.g., Algorithm 7.4.2 in Golub \& Van Loan, <i>%Matrix
   * Computations</i>).
   */
}
class hessenbergDecomposition(val H: householder) {
  val hCoeffs = H.tau
  val matH = H.matrixH

  def reduceToHessenberg() = {
    val icnt = 0
    for (icnt <- 0 to matH.rows - 2)
      H.applyHouseholder(icnt).applyHouseholderRight(icnt).applyHouseholderBottom(icnt)
    H
  }

  def MatrixQ() = householder.householderSequence(matH.mapValues(_.real), hCoeffs.mapValues(_.real), 1)

  /*  4 x 4 example
     *   x    x    x     x
     *   x    x    x    x
     *   0    x     x    x
     *   0    0     x    x
     */
  def MatrixH() = DenseMatrix.tabulate(matH.rows, matH.rows)((i, j) => if (j >= i - 1) matH(i, j) else Complex(0, 0))
}

/**
 * Computes the elementary reflector H such that:
 * \f$ H *this = [ beta 0 ... 0]^T \f$
 * where the transformation H is:
 * \f$ H = I - tau v v^*\f$
 * and the vector v is:
 * \f$ v^T = [1 essential^T] \f$
 *
 * On output:
 * \param essential the essential part of the vector \c v
 * \param tau the scaling factor of the Householder transformation
 * \param beta the result of H * \c *this
 *
 * sa MatrixBase::makeHouseholderInPlace(), MatrixBase::applyHouseholderOnTheLeft(),
 *     MatrixBase::applyHouseholderOnTheRight()
 */

object hessenbergDecomposition {

  def apply(H: householder) = new hessenbergDecomposition(H)
  def apply(M: DenseMatrix[Complex]) =
    {
      //    if (M.rows < 2)
      //       return Nil
      new hessenbergDecomposition(reduceToHessenberg(M))
    }

  def reduceToHessenberg(matA: DenseMatrix[Complex]) = {

    val H: householder = householder(matA)
    val icnt = 0
    for (icnt <- 0 to matA.rows - 2)
      H.applyHouseholder(icnt).applyHouseholderRight(icnt).applyHouseholderBottom(icnt)
    H
}
  }
