package breeze.linalg
/*
 Copyright 2016 @author claydonkey (Anthony Campbell)

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
import breeze.math._
import scala.annotation.tailrec
import breeze.generic.UFunc
import DenseMatrix.canMapValues
import breeze.linalg.Helper._
/*
 * The Householder linear transformation. Describes a reflection about a plane or hyperplane containing the origin.
 * The reflection hyperplane can be defined by a unit vector v (a vector with length 1) which is orthogonal to the hyperplane. The reflection of a point x about this hyperplane is:
* x-2<x,v> v=x-2v( v^H * x)
* where v is given as a column unit vector with Hermitian transpose vH. This is a linear transformation given by the Householder matrix:
* P=I-2vv^H, where I is the identity matrix.
*
* Widely used for tridiagonalization of symmetric matrices and for transforming non-symmetric matrices to a Hessenberg form.
* (https://en.wikipedia.org/wiki/Householder_transformation)
 */
object householder extends UFunc {

  /* takes a Matrix and optional computes a householder Transform. Creates householder class for further householder transformations. The result is stored in matrixH */
  implicit object DMC_DVC_IMPL_H extends Impl3[DenseMatrix[Complex], DenseVector[Complex], Boolean, Householder] {
    def apply(M: DenseMatrix[Complex], tau: DenseVector[Complex], generate: Boolean = true): Householder = {
      new Householder(M, tau, generate)
    }
  }

  /* takes a Matrix and optional computes a householder Transform. Creates householder class for further householder transformations. The result is stored in matrixH */
  implicit object DMC_IMPL_H extends Impl2[DenseMatrix[Complex], Boolean, Householder] {
    def apply(M: DenseMatrix[Complex], generate: Boolean = true): Householder = {
      new Householder(M, generate)
    }
  }

  class Householder(val M: DenseMatrix[Complex], val tau: DenseVector[Complex], val essential: Array[DenseVector[Complex]], generate: Boolean = false) {

    def this(M: DenseMatrix[Complex], tau: DenseVector[Complex], gen: Boolean) = this(M, tau, Array.ofDim[DenseVector[Complex]](M.cols - 1), gen)
    def this(M: DenseMatrix[Complex], gen: Boolean) = this(M, DenseVector.zeros[Complex](M.cols - 1), Array.ofDim[DenseVector[Complex]](M.cols - 1), gen)

    val size: Int = M.cols - 1
    val beta: Array[Double] = Array.ofDim[Double](size)
    val matrixH = M.copy

    if (generate) {
      for (icnt <- 0 to matrixH.rows - 2) {
        makeHouseholder(icnt)
        applyHouseholderOnTheLeft(icnt)
        applyHouseholderOnTheRight(icnt)
    }
      }

    /*
     *  4 x 4 example
     *    x    x    x    x
     *    x    [ Row0 ]
     *    [e]    [  bott   ]
     *    [e]    [   bott  ]
     *  compute A = H A H'
     */
    def makeHouseholder(shift: Int) = {
      essential(shift) = matrixH((shift + 2) to matrixH.rows - 1, shift)
      val eNorm = if (essential(shift).length == 0) 0.0 else sum(essential(shift).map(x => scala.math.pow(x.real, 2) + scala.math.pow(x.imag, 2))) // Does Complex component need squaring?

      val c0 = matrixH(shift + 1, shift);
      (eNorm, c0.imag) match {
        case (0, 0) =>
          beta(shift) = c0.real
          tau(shift) = Complex(0, 0)
        case _ =>
          val c0norm = scala.math.pow(c0.real, 2) + scala.math.pow(c0.imag, 2)
          beta(shift) = if (c0.real >= 0) -Math.sqrt(c0norm + eNorm) else Math.sqrt(c0norm + eNorm)
          essential(shift) = (essential(shift) / (c0 - beta(shift)))
          tau(shift) = conj((beta(shift) - c0) / beta(shift))
      }
      matrixH((shift + 1), shift) = Complex(beta(shift), 0)
      matrixH((shift + 2) until matrixH.rows, shift) := essential(shift)

    }
    /**
     * Apply the elementary reflector H given by
     *  H = I - tau v v^*
     * with
     * v^T = [1 essential^T]
     * from the left to a vector or matrix.
     *  4 x 4 example
     *    x    x    x    x
     *    x    r0  r0  r0
     *    e     Bottom
     *    e     Bottom
     */
    def applyHouseholderOnTheLeft(shift: Int, matShift: Boolean = true) = {

      val matshift = if (matShift) shift else -1

      if (matrixH.rows == 1) {
        matrixH(0, 0) = Complex(1.0, 0.0) - tau(shift)
      } else {
        try {
          val ess = essential(shift)
          var r0 = matrixH((matshift + 1), (matshift + 1) to matrixH.cols - 1).t
          var bottom = matrixH((matshift + 2) to matrixH.rows - 1, (matshift + 1) to matrixH.cols - 1)
          val tmp = (ess.t * bottom) + conj(r0.t)
          r0 -= (tmp * tau(shift)).toDenseVector
          bottom -= (ess * tmp.toDenseMatrix) * tau(shift)
        } catch { case e: Exception => }
      }

    }
    /*
     * Apply the elementary reflector H given by
     *  H = I - tau v v^*
     * with
     *  v^T = [1 essential^T]
     * from the right to a vector or matrix.
     *  4 x 4 example
     *    x    c0    Right
     *    x    c0    Right
     *    e    c0    Right
     *    e    c0    Right
     */
    def applyHouseholderOnTheRight(shift: Int) = {
      try {
        val ess = essential(shift)
        //    if (ess.length != 0) {
        var c1 = matrixH(::, shift + 1).toDenseMatrix
        var right = matrixH(::, (shift + 2) to matrixH.cols - 1)
        val essTrans = essential(shift).t
        val tmp2 = conj((essTrans * right.t) + conj(c1))
        matrixH(0 until matrixH.cols, (shift + 1)) -= (tmp2.toDenseMatrix * conj(tau(shift))).toDenseVector
        right -= conj(tmp2.t * conj(essTrans) * tau(shift))

        //   } else { return this }
      } catch { case e: Exception => }

    }
  }

  //householderTransformation shifted 1 with size -1
  /*  4 x 4 example of the form
     *  1    0    0     0   ^  ------- order (wrapper)
     *   0    1    0    0   v
     *   0    0     x    x
     *   0    0     x    x
     */
  def householderTransformationD(House: Householder, order: Int): DenseMatrix[Double] = {

    val hMatrix = House.matrixH.mapValues(_.real)
    val tau = House.tau.mapValues(_.real)

    if ((hMatrix.rows - order) != tau.length)
      throw new MatrixNotSquareException // change to correct exception

    val matHS = hMatrix(order to (hMatrix.cols - 1), 0 to (hMatrix.rows - order - 1))

    val I = DenseMatrix.eye[Double](matHS.cols)
    val hhMv = (upperTriangular(matHS.t) *:* -(I - 1.0)) +:+ I
    var cnt = 0

    val sum2 = (new ((Int, Int, DenseMatrix[Double]) => DenseMatrix[Double]) {
      @tailrec def apply(from: Int, to: Int, s: DenseMatrix[Double]): DenseMatrix[Double] = {
        if (from == to) return s;
        apply(from + 1, to, s * -(tau(from) * hhMv(from, ::).t * hhMv(from, ::).t.toDenseMatrix - I))
      }
    })(0, matHS.cols - 1, I)

    if (order != 0) {
      val wrapper = DenseMatrix.zeros[Double](hMatrix.cols, hMatrix.cols)
      wrapper(order to wrapper.cols - 1, order to wrapper.cols - 1) := sum2
      for (cnt <- 0 to order - 1) {
        wrapper(cnt, cnt) = 1.0;
      }
      wrapper
    } else
      sum2
  }
  //householderTransformation shifted 1 with size -1
  /*  4 x 4 example of the form
     *  1    0    0     0   ^  ------- order (wrapper)
     *   0    1    0    0   v
     *   0    0     x    x
     *   0    0     x    x
     */
  def householderTransformationC(house: Householder, shift: Int): DenseMatrix[Complex] = {

    var M = DenseMatrix.eye[Complex](house.matrixH.cols)
    var k = 0
    for (k <- (house.matrixH.rows - 2) to 0 by -1) {
      val kshift = house.matrixH.rows - k - shift
      val corner = M((k + 1) to (M.cols) - 1, (k + 1) to (M.cols - 1))
      try {
        corner :=
          {
            new Householder(corner, house.tau.mapValues(i => Complex(i.real, -1 * i.imag)), house.essential, false) {
              applyHouseholderOnTheLeft(k, false)
            }.matrixH
          }
      } catch { case e: Exception => }
    }
    M
  }
}

