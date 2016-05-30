import breeze.linalg._
import scala.io._
import scala.Console._
import breeze.numerics._
import scala.util.control._
import breeze.math._
import scala.util.control.Breaks._
import Helper._

object Householder {
  /*
   *  4 x 4 example
   *    x    x    x    x
   *    x    [ Row0 ]
   *    e    [  bott   ]
   *    e    [   bott  ]
   */
  implicit class IMPL_householder(val matrixH: DenseMatrix[Complex]) {

    val size = matrixH.cols - 1
    val beta = Array.ofDim[Double](size)
    val coeffs = DenseVector.zeros[Complex](size) //tau coeffs
    val essential = Array.ofDim[DenseVector[Complex]](size)

    //implicit def enrichArray[T](xs: Array[T]) = new RichArray[T]
    implicit def enrichDenseMatrix(i: DenseMatrix[Complex]) = new IMPL_householder(i)

    def applyHouseholder(cnt: Int) =
      {
        debugPrint(matrixH, "applyHouseholder start matrixH", 6)
        essential(cnt) = matrixH((cnt + 2) to matrixH.rows - 1, cnt)
        val eNorm = if (essential(cnt).length == 0) 0.0 else sum(essential(cnt).map(x => scala.math.pow(x.real, 2) + scala.math.pow(x.imag, 2))) // Does Complex component need squaring?
        val c0 = matrixH(cnt + 1, cnt);
        (eNorm, c0.imag) match {
          case (0, 0) =>
            beta(cnt) = c0.real
            coeffs(cnt) = Complex(0, 0)
          case _ =>
            val c0norm = scala.math.pow(c0.real, 2) + scala.math.pow(c0.imag, 2)
            beta(cnt) = if (c0.real >= 0) -Math.sqrt(c0norm + eNorm) else Math.sqrt(c0norm + eNorm)
            coeffs(cnt) = ((beta(cnt) - c0) / beta(cnt))
            essential(cnt) = (essential(cnt) / (c0 - beta(cnt)))
        }
        matrixH((cnt + 1), cnt) = Complex(beta(cnt), 0)
        matrixH((cnt + 2) to matrixH.rows - 1, cnt) := essential(cnt)
        val matH2 = matrixH.mapValues(_.real)
        debugPrint(essential, "applyHouseholder essential", 6)
        debugPrint(beta, "applyHouseholder beta", 6)
        debugPrint(coeffs, "applyHouseholder coeffs", 6)
        debugPrint(matrixH, "applyHouseholder matrixH", 6)
        debugPrint(matH2, "applyHouseholder matH2", 6)
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
        matrixH *= 1 - coeffs(cnt)
        this
      } else {
        var c0 = matrixH(::, (cnt + 1) to (cnt + 1))
        var right = matrixH(::, (cnt + 2) to matrixH.cols - 1)
        val tmp = (right * (essential(cnt).toDenseMatrix.map(x => Complex(x.real, -x.imag))).t) + c0
        c0 -= tmp * coeffs(cnt)
        right -= tmp * coeffs(cnt) * essential(cnt).toDenseMatrix
        val matH2 = matrixH.mapValues(_.real)

        debugPrint(matrixH, "applyHouseholderRight matrixH", 6)
        debugPrint(matH2, "applyHouseholderRight matH2", 6)
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
    def applyHouseholderBottom(cnt: Int): IMPL_householder = {

      if (matrixH.cols == 1) {
        matrixH *= 1 - coeffs(cnt)
        this
      } else {
        var r0 = matrixH((cnt + 1), (cnt + 1) to matrixH.cols - 1).t
        var bottom = matrixH((cnt + 2) to matrixH.rows - 1, (cnt + 1) to matrixH.cols - 1)
        val tmp = (bottom.t * essential(cnt).map(x => Complex(x.real, -x.imag))) + r0
        r0 -= (tmp.toDenseMatrix * coeffs(cnt)).toDenseVector
        bottom -= (essential(cnt).toDenseMatrix.t * tmp.toDenseMatrix) * coeffs(cnt)
        val matH2 = matrixH.mapValues(_.real)
        debugPrint(matrixH, "applyHouseholderBottom matrixH", 6)
        debugPrint(matH2, "applyHouseholderBottom matH2", 6)
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

    //def apply(m_matrix: DenseMatrix[Complex]): householder = new householder(m_matrix)

    def householderSequence(hMatrix: DenseMatrix[Double], hCoeffs: DenseVector[Double], order: Int): DenseMatrix[Double] = {
      //HouseholderSequence shifted 1 with size -1
      /*  4 x 4 example of the form
     *  1    0    0     0   ^  ------- order (wrapper)
     *   0    1    0    0   v
     *   0    0     x    x
     *   0    0     x    x
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

      if (order != 0) {
        val wrapper = DenseMatrix.zeros[Double](hMatrix.cols, hMatrix.cols)
        wrapper(order to wrapper.cols - 1, order to wrapper.cols - 1) := sum
        for (cnt <- 0 to order - 1) {
          wrapper(cnt, cnt) = 1.0;

        }
        wrapper
      } else
        sum
    }
  }
}