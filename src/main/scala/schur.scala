import breeze.linalg._
import scala.Console._
import breeze.numerics._
import scala.util.control._
import breeze.math._
import scala.util.control.Breaks._
import Helper._
import jacobi._
import Hessenberg._
import Householder._

/* Schur decomposition  -computed by first reducing the
 * matrix to Hessenberg form using the class
 * HessenbergDecomposition. The Hessenberg matrix is then reduced
 * to triangular form by performing QR iterations with a single
 * shift. The cost of computing the Schur decomposition depends
 * on the number of iterations; as a rough guide, it may be taken
 * on the number of iterations; as a rough guide, it may be taken
 * to be \f$25n^3\f$ complex flops, or \f$10n^3\f$ complex flops
 * if \a computeU is false.*/

object Schur {

  abstract class Schur[T](val M: DenseMatrix[T]) {
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
    val hessenberg: Tuple2[hessenberg[T], Householder]

    lazy val matT = hessenberg._1.MatrixH() //matT is  is an Upper Triangle
    lazy val matQ = hessenberg._1.MatrixP()

    def decomposition() = (matT, matQ, hessenberg._2.tau, hessenberg._2.matrixH)
  }

  class SchurLAPACKD(M: DenseMatrix[Double]) extends Schur[Double](M) {
    lazy override val hessenberg = M hessenbergDecomp
    val (s, z, y2, wR, wI, sLO, sHI) = schur(M)

    debugPrint(y2, "y2", 1)
            debugPrint(hessenberg._1.MatrixP(), "hessenberg._1.MatrixP()", 1)
    lazy override val matT = DenseMatrix.tabulate[Complex](M.cols, M.rows)((i, j) => Complex(s(i, j), 0.0)) //z(i, j)))
    lazy override val matQ = hessenberg._1.MatrixP() * DenseMatrix.tabulate[Complex](M.cols, M.rows)((i, j) => Complex(z(i, j), 0.0)) 
  }

  class SchurLAPACKI(M: DenseMatrix[Int]) extends Schur[Double](M.mapValues(_.toDouble)) {
    val MD = M.mapValues(_.toDouble)
    lazy override val hessenberg = MD hessenbergDecomp
    val (s, z, y2, wR, wI, sLO, sHI) = schur(MD)
    debugPrint(y2, "y2", 1)
        debugPrint(hessenberg._1.MatrixP(), "hessenberg._1.MatrixP()", 1)
    lazy override val matT = DenseMatrix.tabulate[Complex](M.cols, M.rows)((i, j) => Complex(s(i, j), 0.0))
    lazy override val matQ = hessenberg._1.MatrixP() * DenseMatrix.tabulate[Complex](M.cols, M.rows)((i, j) => Complex(z(i, j), 0.0))
  }

  class SchurC(M: DenseMatrix[Complex]) extends Schur[Complex](M) {
    override val hessenberg = M hessenbergDecomp

    val m_maxIterationsPerRow = 30
    val maxIters = m_maxIterationsPerRow * hessenberg._2.matrixH.rows

    reduceToTriangularForm()
    /**
     * If matT(i+1,i) is neglegible in floating point arithmetic
     * compared to matT(i,i) and matT(j,j), then set it to zero and
     * return true, else return false.
     */
    def subdiagonalEntryIsNeglegible(i: Int) =
      {
        val d = norm1(matT(i, i)) + norm1(matT(i + 1, i + 1))
        val sd = norm1(matT(i + 1, i))
        if (isMuchSmallerThan(sd, d)) {
          matT(i + 1, i) = Complex(0.0, 0.0)
          debugPrint("SMALL VALUE", "SMALL VALUE", 0)
          true
        } else
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

      var t = matT((iu - 1) to iu, (iu - 1) to iu) // Complex(NormT, 0)
      val normt = Complex(sum(t.mapValues(abs(_))), 0)

      debugPrint(t, "t", 3)
      debugPrint(t.mapValues(abs(_)), "t.mapValues(abs(_))", 3)

      debugPrint(normt, "normt", 3)
      t = t / normt

      val b = t(0, 1) * t(1, 0)

      val c = t(0, 0) - t(1, 1)
      val disc = breeze.numerics.pow((c * c + 4.0 * b), 0.5)
      val det = (t(0, 0) * t(1, 1)) - b

      debugPrint(t(0, 0) * t(1, 1), "t(0, 0) * t(1, 1)", 3)
      debugPrint(b, "b", 3)
      debugPrint(det, "det", 3)

      val trace = t(0, 0) + t(1, 1)
      var eival1 = (trace + disc) / 2.0
      var eival2 = (trace - disc) / 2.0

      if (norm1(eival1) > norm1(eival2))
        eival2 = det / eival1
      else
        eival1 = det / eival2

      debugPrint(eival1, "eival1", 3)
      debugPrint(eival2, "eival2", 3)
      // choose the eigenvalue closest to the bottom entry of the diagonal
      if (norm1(eival1 - t(1, 1)) < norm1(eival2 - t(1, 1)))
        normt * eival1
      else
        normt * eival2

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

      var iu = matT.cols - 1
      val maxIters = m_maxIterationsPerRow * matT.rows
      var il = 0
      var iter = 0 // number of iterations we are working on the (iu,iu) element
      var totalIter = 0 // number of iterations for whole matrix

      debugPrint(matT, "Start MatT", 0)
      debugPrint(matQ, "Start MatQ", 0)
      breakable {
        while (true) {
          // find iu, the bottom row of the active submatrix
          breakable {
            while (iu > 0) {
              if (subdiagonalEntryIsNeglegible(iu - 1)) {
                debugPrint(matT, "Negligible MatT", 0)
              }
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
          debugPrint(shift, " shift ", 3)
          val rot = makeGivens(matT(il, il) - shift, matT(il + 1, il))

          matT(il to il + 1, ::) rotateoL (rot)
          debugPrint(matT, "MatT A", 0)
          matT(0 to (min(il + 2, iu)), il to il + 1) rotateoR (rot)
          if (il + 2 > iu) debugPrint(iu, "lesser", 0);
          debugPrint(matT, "MatT B", 0)
          matQ(::, il to il + 1) rotateoR (rot)

          debugPrint(matQ, "MatQ A", 0)

          val idx: Int = 0

          for (idx <- ((il + 1) to iu - 1)) {

            val rot2 = makeGivens(matT(idx, idx - 1), matT(idx + 1, idx - 1))
            matT(idx, idx - 1) = rot2.rot
            matT(idx + 1, idx - 1) = Complex(0.0, 0.0)
            matT(idx to idx + 1, idx to matT.cols - 1) rotateoL (rot2)
            debugPrint(matT, "MatT 1", 0)
            if (idx + 2 > iu) debugPrint(iu, "lesser", 0);
            matT(0 to (min(idx + 2, iu)), idx to idx + 1) rotateoR (rot2)
            //CHECK FOR 3 needing to be equal to cell size..
            debugPrint(matT, "MatT 2", 0)
            matQ(::, idx to idx + 1) rotateoR (rot2)
            debugPrint(matQ, "MatQ 1", 0)

          }

        }
      }
    }
  }
  implicit def schurDecompositionM(M: DenseMatrix[Complex]) = new { def schurDecomp = { new SchurC(M).decomposition() } }

  implicit def schurDecompositionD(M: DenseMatrix[Double]) = new { def schurDecomp = { new SchurLAPACKD(M).decomposition() } }

  implicit def schurDecomposI(M: DenseMatrix[Int]) = new { def schurDecomp = { new SchurLAPACKI(M).decomposition() } }

  implicit def IMPL_Schur2(M: DenseMatrix[Complex]) = new {
    def revertSchur(T: DenseMatrix[Complex]): DenseMatrix[Complex] = {

      debugPrint(Tadj(T), "Tadj(T)", 1)
      debugPrint(upperTriangular(M), "upperTriangular(M) ", 1)
      debugPrint(T, "T ", 1)
      return T * upperTriangular(M) * Tadj(T)
  }
    }

}
