import Householder._
import breeze.linalg._
import scala.io._
import scala.Console._
import breeze.numerics._
import scala.util.control._
import breeze.math._
import scala.util.control.Breaks._
import Helper._

object Hessenberg {

  //PHP^(H)=A,
  /*
A Hessenberg decomposition is a matrix decomposition of a matrix A into a unitary matrix P and a Hessenberg matrix H such that
 PHP^(H)=A,
where P^(H) denotes the conjugate transpose.
Hessenberg decomposition is implemented in the Wolfram Language as HessenbergDecomposition[m].
Hessenberg decomposition is the first step in Schur decomposition. Hessenberg decomposition on an n×n matrix requires 14n^3/3 arithmetic operations.
*/

  implicit class IMPL_hessenbergDecomposition(M: DenseMatrix[Complex]) {
    val House = hessenbergDecomposition.reduceToHessenberg(M)

    val tau = House.tau
    val matH = House.matrixH

    def getHessenberg() = this

    /*  4 x 4 example of the form
     *  1    0    0     0   ^  ------- order
     *   0    1    0    0   v
     *   0    0     x    x
     *   0    0     x    x
     */
    def MatrixP() = householder.householderSequence(matH.mapValues(_.real), tau.mapValues(_.real), 1)

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
    * H * M = [ beta 0 ... 0]^T
    * where the transformation H is:
    * H = I - tau v v^*
    * and the vector v is:
    * v^T = [1 essential^T]
    *
    * Householder obj variables:
    * \ essential the essential part of the vector  v
    * \ tau the scaling factor of the Householder transformation
    * \ beta the result of H * M
    */

  object hessenbergDecomposition {

    def reduceToHessenberg(M: DenseMatrix[Complex]) = {

      val H: IMPL_householder = IMPL_householder(M.copy)
      val icnt = 0
      for (icnt <- 0 to M.rows - 2)
        H.makeHouseholder(icnt).applyHouseholderRight(icnt).applyHouseholderBottom(icnt)
    
      H

    }
  }

}