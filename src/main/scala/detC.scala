/**
 * Computes the determinant of the given real matrix.
 */
import breeze.math.Complex
import spire.implicits.cforRange
import breeze.generic._
import breeze.linalg._

object detC extends UFunc {
  implicit def canDetUsingLU2[T](implicit luImpl: LU.Impl[T, (DenseMatrix[Complex], Array[Int])]):Impl[T, Complex] = {
    new Impl[T, Complex] {
      def apply(X: T): Complex = {

        // For triangular N-by-N matrices X, the determinant of X equals the product
        // of the diagonal elements X(i,i) where 0 <= i < N.
        // Since det(AB) = det(A) * det(B), the LU factorization is well-suited for
        // the computation of the determinant of general N-by-N matrices.
        val (m:DenseMatrix[Complex], ipiv:Array[Int]) = LU(X)

        // Count the number of exchanged rows.  ipiv contains an array of swapped indices,
        //  but it also contains indices that weren't swapped.  To count the swapped
        //  indices, we have to compare them against their position within the array.  A
        //  final complication is that the array indices are 1-based, due to the LU call
        //  into LAPACK.
        val numExchangedRows = ipiv.map(_ - 1).zipWithIndex.count { piv => piv._1 != piv._2 }

        var acc = if (numExchangedRows % 2 == 1) Complex(-1.0,0) else Complex(1.0,0)
        cforRange(0 until m.rows){ i =>
          acc *= m(i,i)
        }

        acc
      }

    }
  }
}