import breeze.linalg._
import scala.io._
import java.text.DecimalFormat
import scala.Console._
import breeze.numerics._
import scala.util.control._
import breeze.math._
import scala.util.control.Breaks._
import Helper._
/* The Schur decomposition is computed by first reducing the
   * matrix to Hessenberg form using the class
   * HessenbergDecomposition. The Hessenberg matrix is then reduced
   * to triangular form by performing QR iterations with a single
   * shift. The cost of computing the Schur decomposition depends
   * on the number of iterations; as a rough guide, it may be taken
   * on the number of iterations; as a rough guide, it may be taken
   * to be \f$25n^3\f$ complex flops, or \f$10n^3\f$ complex flops
   * if \a computeU is false.*/
class complexSchur(val hess: hessenbergDecomposition) {

  val matT = hess.MatrixH
 private  val matU = hess.MatrixQ
  val matUC = matU.mapValues(Complex(_, 0.0))
  val m_maxIterationsPerRow = 30
  val maxIters = m_maxIterationsPerRow * hess.matH.rows
  // reduceToTriangularForm()
  /**
   * If matT(i+1,i) is neglegible in floating point arithmetic
   * compared to matT(i,i) and matT(j,j), then set it to zero and
   * return true, else return false.
   */
  val EPSILON: Double = 2.22045e-016

  def isMuchSmallerThan(x: Double, y: Double) = { abs(x) <= abs(y) * EPSILON }

  def subdiagonalEntryIsNeglegible(i: Int) =
    {
      val d = norm1(matT(i, i)) + norm1(matT(i + 1, i + 1))
      val sd = norm1(matT(i + 1, i))
      if (isMuchSmallerThan(sd, d)) {
        val t = matT(i + 1, i)
        val M = matT.mapValues(_.real)
        debugPrint(s"SMALLER   i = $i matT($i + 1, $i) :  $t")
        matT(i + 1, i) = Complex(0, 0)

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
    val NormT = sum(matT((iu - 1) to matT.cols - 1, (iu - 1) to matT.cols - 1).mapValues(abs(_)))
    val T = matT((iu - 1) to matT.cols - 1, (iu - 1) to matT.cols - 1) / Complex(NormT, 0)

    //  debugPrint(s"T\n $T\n")
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
    //  debugPrint(s"eival  $eival1 \n $eival2 \n")
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
    //
    //
    //
    //

    debugPrint("reduceToTriangularForm" , "")

    debugPrint( matT,"Start MatT")
    debugPrint(matU, "Start MatU")

    var matnum = 0
    var matnum2 = 0

    var iu = matT.cols - 1
    val maxIters = m_maxIterationsPerRow * matT.rows
    var il = 0
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
        //       debugPrint("Prior or with Zeroing:\n\n" + matT.mapValues(x => x.real - (x.real % 0.0001)))
        val shift = computeShift(iu, iter)
        var rot = jacobiRotation.makeGivens(matT(il, il) - shift, matT(il + 1, il))

        val A = jacobiRotation.applyRotationinPlane(matT(il, il to matT.cols - 1).t, matT(il + 1, il to matT.cols - 1).t, rot.adjoint)
        //      debugPrint("J1\n" + J1.mapValues(x => x.real - (x.real % 0.0001)))
        matT(il to il + 1, il to matT.cols - 1) := A
        debugPrint(matT, "A MatT")

        val B = jacobiRotation.applyRotationinPlane(matT(0 to matT.cols - 2 - il, il), matT(0 to matT.cols - 2 - il, il + 1), rot.transpose).t
        //    debugPrint("J2\n" + J2.mapValues(x => x.real - (x.real % 0.0001)))
        matT(il to matT.rows - 2, 0 to 1) := B
        debugPrint(matT, "B MatT")

        val AmatUC = jacobiRotation.applyRotationinPlane(matUC(0 to matUC.cols - 1 - il, il), matUC(0 to matUC.cols - 1 - il, il + 1), rot.transpose).t
        // debugPrint("matUc J3\n" + J3.mapValues(x => x.real - (x.real % 0.0001)))

        matUC(il to matUC.cols - 1, 0 to 1) := AmatUC
        //    debugPrint("matUC\n" + matUC.mapValues(x => x.real - (x.real % 0.0001)))

        var idx: Int = 0
        for (idx <- ((il + 1) to iu - 1)) {

          rot = rot.makeGivens(matT(idx, idx - 1), matT(idx + 1, idx - 1))
          matT(idx, idx - 1) = rot.r
          debugPrint(s"rot.r  into $idx,$idx -1\n" + rot.r + "\n")
          matT(idx + 1, idx - 1) = Complex(0.0, 0.0)
          //      debugPrint(s"$idx\n")

          val J1 = jacobiRotation.applyRotationinPlane(matT(idx, idx to matT.cols - 1).t, matT(idx + 1, idx to matT.cols - 1).t, rot.adjoint)
          //      debugPrint("J4\n" + J4.mapValues(x => x.real - (x.real % 0.0001)))
          matT(idx to idx + 1, idx to matT.cols - 1) := J1
          debugPrint( matT,"1 MatT")

          val J2 = jacobiRotation.applyRotationinPlane(matT(0 to matT.cols - 1, idx), matT(0 to matT.cols - 1, idx + 1), rot.transpose).t
          //     debugPrint("J5\n" + J5.mapValues(x => x.real - (x.real % 0.0001)))
          matT(0 to matT.rows - 1, idx to idx + 1) := J2 //CHECK FOR 3 needing to be equal to cell size..
          debugPrint(matT,"2 MatT")

          val J3 = jacobiRotation.applyRotationinPlane(matUC(0 to matUC.cols - 1, idx), matUC(0 to matUC.cols - 1, idx + 1), rot.transpose).t
          //     debugPrint("J6\n" + J6.mapValues(x => x.real - (x.real % 0.0001)))
          matUC(0 to matUC.rows - 1, idx to idx + 1) := J3
          debugPrint(matU,"1 MatU")
        }

      }
    }
}
  }

object complexSchur {

  def apply(M: DenseMatrix[Complex]): complexSchur = {

    val S = new complexSchur(hessenbergDecomposition(M.copy))
    S.reduceToTriangularForm
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