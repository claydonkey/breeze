import breeze.linalg._
import java.text.DecimalFormat
import scala.Console._
import breeze.numerics._
import scala.util.control._
import breeze.math._
import scala.util.control.Breaks._
import Helper._
import jacobiRotation._
import Hessenberg._
/* The Schur decomposition is computed by first reducing the
 * matrix to Hessenberg form using the class
 * HessenbergDecomposition. The Hessenberg matrix is then reduced
 * to triangular form by performing QR iterations with a single
 * shift. The cost of computing the Schur decomposition depends
 * on the number of iterations; as a rough guide, it may be taken
 * on the number of iterations; as a rough guide, it may be taken
 * to be \f$25n^3\f$ complex flops, or \f$10n^3\f$ complex flops
 * if \a computeU is false.*/

object Schur {

  implicit class complexSchur(val M: DenseMatrix[Complex]) {
    if (M.rows != M.cols)
      throw new MatrixNotSquareException
    /*
     *  Hessenberg decomposition of given matrix.
     *
     * The Hessenberg decomposition is computed by bringing the columns of the
     * matrix successively in the required form using Householder reflections
     * (see, e.g., Algorithm 7.4.2 in Golub \& Van Loan, <i>%Matrix
     * Computations</i>).
     */
    val hess = M.getHessenberg()
    val matU = hess.MatrixH
    val matQ = hess.MatrixP.mapValues(Complex(_, 0.0))
    val m_maxIterationsPerRow = 30
    val maxIters = m_maxIterationsPerRow * hess.matH.rows
    reduceToTriangularForm()
    /**
     * If matT(i+1,i) is neglegible in floating point arithmetic
     * compared to matT(i,i) and matT(j,j), then set it to zero and
     * return true, else return false.
     */


    def getSchur() =
      {
        (matU, matQ, hess.MatrixP, hess.MatrixH, hess.hCoeffs)
      }

    def subdiagonalEntryIsNeglegible(i: Int) =
      {
        val d = norm1(matU(i, i)) + norm1(matU(i + 1, i + 1))
        val sd = norm1(matU(i + 1, i))
        if (isMuchSmallerThan(sd, d)) {
          val t = matU(i + 1, i)
          val M = matU.mapValues(_.real)
          matU(i + 1, i) = Complex(0, 0)

          true
        } else
          false
      }
    /*def revertSchur( T: DenseMatrix[Complex]) {
     M * (upperTriangular(T) * adj(M))
     }
     */
    def revertSchur(T: DenseMatrix[Complex]) = {

      M * upperTriangular(T) * Helper.adj(M)

    }

    /** Compute the shift in the current QR iteration. */
    def computeShift(iu: Int, iter: Int) = {
      if (iter == 10 || iter == 20) {
        // exceptional shift, taken from http://www.netlib.org/eispack/comqr.f
        abs(matU(iu, iu - 1).real) + abs(matU(iu - 1, iu - 2).real)
      }
      // compute the shift as one of the eigenvalues of t, the 2x2
      // diagonal block on the bottom of the active submatrix

      var t = matU((iu - 1) to matU.cols - 1, (iu - 1) to matU.cols - 1) // Complex(NormT, 0)
      val normt = Complex(sum(t.mapValues(abs(_))), 0)

      t = t / normt

      val b = t(0, 1) * t(1, 0)
      val c = t(0, 0) - t(1, 1)
      val disc = breeze.numerics.pow((c * c + 4.0 * b), 0.5)
      val det = t(0, 0) * t(1, 1) - b
      val trace = t(0, 0) + t(1, 1)
      var eival1 = (trace + disc) / 2.0
      var eival2 = (trace - disc) / 2.0;

      if (norm1(eival1) > norm1(eival2))
        eival2 = det / eival1;
      else
        eival1 = det / eival2;

      // choose the eigenvalue closest to the bottom entry of the diagonal
      if (norm1(eival1 - t(1, 1)) < norm1(eival2 - t(1, 1)))
        normt * eival1;
      else
        normt * eival2;

    }

    def reduceToTriangularForm() = {
      /**A horrbly disFunctional implementation using breaks, iterations  and mutations BLEUURRGGH will change... **/
      // The matrix matT is divided in three parts.
      // Rows 0,...,il-1 are decoupled from the rest because matT(il,il-1) is zero.
      // Rows il,...,iu is the part we are working on (the active submatrix).
      // Rows iu+1,...,end are already brought in triangular form.

      var matnum = 0
      var matnum2 = 0
      var newrot = false

      var iu = matU.cols - 1
      val maxIters = m_maxIterationsPerRow * matU.rows
      var il = 0
      var iter = 0 // number of iterations we are working on the (iu,iu) element
      var totalIter = 0 // number of iterations for whole matrix

      debugPrint(matU, "Start MatU", 4)

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

          debugPrint(il, "il ", 3)
          debugPrint(iu, "iu ", 3)

          val shift = computeShift(iu, iter)
          debugPrint(shift, "A shift", 3)

         implicit  val rot = makeGivens(matU(il, il) - shift, matU(il + 1, il))

          rot.applyOnTheLeft(matU(il, ::), matU(il + 1, ::))
          debugPrint(matU, "A MatT", 4)
          rot.applyOnTheRight(matU(0 to matU.cols - il - 2, il), matU(0 to matU.cols - il - 2, il + 1))
          debugPrint(matU, "B MatT", 4)
          rot.applyOnTheRight(matQ(::, il), matQ(::, il + 1))
          debugPrint(matQ, "AQ", 0)

          val idx: Int = 0

          for (idx <- ((il + 1) to iu - 1)) {

            val rot2 = makeGivens(matU(idx, idx - 1), matU(idx + 1, idx - 1))
            matU(idx, idx - 1) = rot2.rot
            matU(idx + 1, idx - 1) = Complex(0.0, 0.0)
            rot2.applyOnTheLeft(matU(idx, idx to matU.cols - 1), matU(idx + 1, idx to matU.cols - 1))

            debugPrint(matU, "1 MatT", 4)
            debugPrint(idx, "idx", 3)

            rot2.applyOnTheRight(matU(::, idx), matU(::, idx + 1))
            //CHECK FOR 3 needing to be equal to cell size..
            debugPrint(matU, "2 MatT", 4)
            rot2.applyOnTheRight(matQ(::, idx), matQ(::, idx + 1))
            debugPrint(matQ, "1 MatQ", 0)

          }

        }
      }
    }
  }

}