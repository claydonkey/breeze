import Householder._
import breeze.linalg._
import scala.io._
import scala.Console._
import breeze.numerics._
import scala.util.control._
import breeze.math._
import scala.util.control.Breaks._
import Helper._

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


object Hessenberg {

  /*
  PHP^(H)=A,

  A Hessenberg decomposition is a matrix decomposition of a matrix A into a unitary matrix P and a Hessenberg matrix H such that
  PHP^(H)=A,
  where P^(H) denotes the conjugate transpose.
  Hessenberg decomposition is implemented in the Wolfram Language as HessenbergDecomposition[m].
  Hessenberg decomposition is the first step in Schur decomposition. Hessenberg decomposition on an n×n matrix requires 14n^3/3 arithmetic operations.

  Computes the elementary reflector H such that:
  H * M = [ beta 0 ... 0]^T
  where the transformation H is:
  H = I - tau v v^*
  and the vector v is:
  v^T = [1 essential^T]

  Householder obj variables:
  essential the essential part of the vector  v
  House.tau  the scaling factor of the Householder transformation
  House.beta the result of H * M
 */

  implicit class hessenbergDecompositionC(M: DenseMatrix[Complex]) extends hessenberg[Complex](M) {
    override val House = {
      val H = new Householder(M.copy, DenseVector.zeros[Complex](M.cols - 1))
      val icnt = 0
      //LAPACK misses this step
      for (icnt <- 0 to M.rows - 2)
        H.makeHouseholder(icnt).applyHouseholderRight(icnt).applyHouseholderBottom(icnt)
      H
    }
  }

  implicit class hessenbergDecompositionD(M: DenseMatrix[Double]) extends hessenberg[Double](M) {
    override val House = {
      val (h, hLO, hHI, tau) = hessenberg(M)
      new Householder(h.mapValues(Complex(_, 0.0)), DenseVector.tabulate[Complex](M.cols - 1)(i => Complex(tau(i), 0.0)))
    }
  }

  implicit class hessenbergDecompositionI(M: DenseMatrix[Int]) extends hessenberg[Int](M) {
    override val House = {
      val (h, hLO, hHI, tau) = hessenberg(M.mapValues(_.toDouble))
      new Householder(h.mapValues(Complex(_, 0.0)), DenseVector.tabulate[Complex](M.cols - 1)(i => Complex(tau(i), 0.0)))
    }
  }

  abstract class hessenberg[T](M: DenseMatrix[T]) {

    val House: Householder

    def hessenbergDecomp() = (this, House)

    /*  4 x 4 example of the form
     *  1    0    0     0   ^  ------- order
     *   0    1    0    0   v
     *   0    0     x    x
     *   0    0     x    x
     */
    def MatrixP(): DenseMatrix[Complex] = householder.householderTransformation(House.matrixH.mapValues(_.real), House.tau.mapValues(_.real), 1).mapValues(Complex(_, 0.0))

    /*  4 x 4 example
     *   x    x    x     x
     *   x    x    x    x
     *   0    x     x    x
     *   0    0     x    x
     */
    def MatrixH(): DenseMatrix[Complex] = DenseMatrix.tabulate(House.matrixH.rows, House.matrixH.rows)((i, j) => if (j >= i - 1) House.matrixH(i, j) else Complex(0, 0))
  }

}