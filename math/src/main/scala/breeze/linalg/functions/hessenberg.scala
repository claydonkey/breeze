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
import org.netlib.util.intW
import com.github.fommil.netlib.LAPACK.{ getInstance => lapack }
import breeze.generic.UFunc
import breeze.math._
import breeze.linalg.householder._
import DenseMatrix.canMapValues
import reflect.runtime.universe._

/*
 * A Hessenberg decomposition is a matrix decomposition of a matrix A into a unitary matrix P and a Hessenberg matrix H
 * such that PHP^H=A,  where P^H denotes the conjugate transpose.
 * Hessenberg decomposition is the first step in Schur decomposition.
 * the factory method "hessenberg" returns a tuple from hessenberg  decomposition where
* p is a unitary matrix
* h is an Hessenberg  matrix,
*/
object hessenberg extends UFunc {

  implicit object DMC_IMPL_HeCHo extends Impl[DenseMatrix[Complex], (DenseMatrix[Complex], DenseMatrix[Complex])] {
    def apply(M: DenseMatrix[Complex]): (DenseMatrix[Complex], DenseMatrix[Complex]) = {
      new Hessenberg[Complex](M, new Householder(M.copy).generate()).decompose
    }
  }

  implicit object DMD_IMPL_HeDHo extends Impl[DenseMatrix[Double], (DenseMatrix[Complex], DenseMatrix[Complex])] {
    def apply(M: DenseMatrix[Double]): (DenseMatrix[Complex], DenseMatrix[Complex]) = {
      val (h, hLO, hHI, tau) = getHessenbergLAPACK(M)
      new Hessenberg[Double](M, new Householder(h.mapValues(Complex(_, 0.0)), DenseVector.tabulate[Complex](M.cols - 1)(i => Complex(tau(i), 0.0)))).decompose()
    }
  }

  implicit object DMD_IMPL_DMI_DMI_Ho extends Impl[DenseMatrix[Int], (DenseMatrix[Complex], DenseMatrix[Complex])] {
    def apply(M: DenseMatrix[Int]): (DenseMatrix[Complex], DenseMatrix[Complex]) = {
      val (h, hLO, hHI, tau) = getHessenbergLAPACK(M.mapValues(_.toDouble))
      new Hessenberg[Int](M, new Householder(h.mapValues(Complex(_, 0.0)), DenseVector.tabulate[Complex](M.cols - 1)(i => Complex(tau(i), 0.0)))).decompose()
    }
  }
  class Hessenberg[T: TypeTag](M: DenseMatrix[T], val House: Householder[T]) {

    def decompose() = (P, H)

    /*  4 x 4 example of the form
     *  1    0    0     0   ^  ------- order
     *   0    1    0    0   v
     *   0    0     x    x
     *   0    0     x    x
     */
    def P = House.P(1)

    /*  4 x 4 example tridiagonal matrix
     *   x    x    0     0
     *   x    x    x    0
     *   0    x     x    x
     *   0    0     x    x
     */
    //This does nothing except zero the corners     i.e H = householder.triDiagonalize(M)
    def H = DenseMatrix.tabulate(House.matrixH.rows, House.matrixH.rows)((i, j) => if (j >= i - 1) House.matrixH(i, j) else Complex(0, 0))

  }

  def getHessenbergLAPACK(X: DenseMatrix[Double]): (DenseMatrix[Double], Int, Int, Array[Double]) = {

    val M = X.rows
    val N = X.cols
    val y = X.copy
    var iLO = new intW(0)
    var iHI = new intW(0)
    var iScale = Array.ofDim[Double](N)
    var tau = Array.ofDim[Double](N)
    var info = new intW(0)

    lapack.dgebal(
      "P",
      N,
      y.data,
      scala.math.max(1, M),
      iLO,
      iHI,
      iScale,
      info
    )

    if (info.`val` > 0)
      throw new NotConvergedException(NotConvergedException.Iterations)
    else if (info.`val` < 0)
      throw new IllegalArgumentException()

    var workcnt = Array.ofDim[Double](1)

    lapack.dgehrd(N, iLO.`val`, iHI.`val`, y.data, scala.math.max(1, N), tau, workcnt, -1, info)
    lapack.dgehrd(N, iLO.`val`, iHI.`val`, y.data, scala.math.max(1, N), tau, Array.ofDim[Double](workcnt(0).toInt), workcnt(0).toInt, info)

    //Error check
    if (info.`val` > 0)
      throw new NotConvergedException(NotConvergedException.Iterations)
    else if (info.`val` < 0)
      throw new IllegalArgumentException()

    (y, iLO.`val`, iHI.`val`, tau)
  }

}