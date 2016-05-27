import Householder._
import breeze.linalg._
import scala.io._
import scala.Console._
import breeze.numerics._
import scala.util.control._
import breeze.math._
import scala.util.control.Breaks._

object Hessenberg {
  //PHP^(H)=A,
  /*
A Hessenberg decomposition is a matrix decomposition of a matrix A into a unitary matrix P and a Hessenberg matrix H such that
 PHP^(H)=A,
where P^(H) denotes the conjugate transpose.
Hessenberg decomposition is implemented in the Wolfram Language as HessenbergDecomposition[m].
Hessenberg decomposition is the first step in Schur decomposition. Hessenberg decomposition on an n×n matrix requires 14n^3/3 arithmetic operations.
*/

  implicit class hessenbergDecomposition(M: DenseMatrix[Complex]) {
val H = hessenbergDecomposition.reduceToHessenberg2(M.copy)

    val hCoeffs = H.coeffs
    val matH = H.matrixH

    def reduceToHessenberg() = {
      val icnt = 0
      for (icnt <- 0 to matH.rows - 2)
        H.applyHouseholder(icnt).applyHouseholderRight(icnt).applyHouseholderBottom(icnt)
    this
    }

    /*  4 x 4 example of the form
     *  1    0    0     0   ^  ------- order
     *   0    1    0    0   v
     *   0    0     x    x
     *   0    0     x    x
     */
    def MatrixP() = householder.householderSequence(matH.mapValues(_.real), hCoeffs.mapValues(_.real), 1)

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

    // def apply(H: householder) = new hessenbergDecomposition(H)
    /* def apply(M: DenseMatrix[Complex]) =
    {

      //    if (M.rows < 2)
      //       return Nil
      new hessenbergDecomposition(reduceToHessenberg(M.copy))

    }
*/
    def reduceToHessenberg2(M: DenseMatrix[Complex]) = {

      //  val H: householder = householder(M.copy)
      val M2 = M.copy
      val icnt = 0
      for (icnt <- 0 to M.rows - 2)
        M2.applyHouseholder(icnt).applyHouseholderRight(icnt).applyHouseholderBottom(icnt)
      M2

    }
  }

}