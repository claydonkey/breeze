import breeze.linalg._
import scala.io._
import scala.Console._
import breeze.numerics._
import scala.util.control._
import breeze.math._
import scala.util.control.Breaks._
import Helper._

/* This class represents a Jacobi or Givens rotation.
   * This is a 2D rotation in the plane \c J of angle \f$ \theta \f$ defined by
   * its cosine \c c and sine \c s as follow:
   * \f$ J = \left ( \begin{array}{cc} c & \overline s \\ -s  & \overline c \end{array} \right ) \f$
   *
   * You can apply the respective counter-clockwise rotation to a column vector \c v by
   * applying its adjoint on the left: \f$ v = J^* v \f$
   */
class jacobiRotation(val m_c: Complex, val m_s: Complex, var r: Complex) {

  def this(m_c: Complex, m_s: Complex) = this(m_c, m_s, Complex(0.0, 0.0))

  def transpose() = jacobiRotation(m_c, -conj(m_s))
  def adjoint() = jacobiRotation(conj(m_c), -m_s)
  def applyRotationinPlane(_x: DenseVector[Complex], _y: DenseVector[Complex]) = {
    jacobiRotation.applyRotationinPlane(_x, _y, this)
  }

  /*This function implements the continuous Givens rotation
     *generation algorithm found in Anderson (2000),
     *Discontinuous Plane Rotations and the Symmetric Eigenvalue Problem.
     *LAPACK Working Note 150, University of Tennessee, UT-CS-00-454, December 4, 2000. */
  def makeGivens(p: Complex, q: Complex) = {

    (p, q) match {
      case (Complex(0.0, 0.0), _) =>
        val m_c = if (p.real < 0) Complex(-1, 0) else Complex(1, 0)
        val m_s = Complex(0, 0)
        r = m_c * p;
        new jacobiRotation(m_c, m_s, r)
      case (_, Complex(0.0, 0)) =>
        val m_c = Complex(0, 0)
        val m_s = -q / abs(q)
        r = Complex(abs(q), 0)
        new jacobiRotation(m_c, m_s, r)
      case _ =>
        var p1 = norm1(p)
        var q1 = norm1(q)
        if (p1 >= q1) {
          val ps = p / p1
          val p2 = abs2(ps)
          val qs = q / p1
          val q2 = abs2(qs)

          var u = breeze.numerics.pow(Complex(1.0, 0.0) + (q2 / p2), 0.5)
          if (p.real < 0)
            u = -u;
          val m_c = 1 / u
          val m_s = -qs * conj(ps) * (m_c / p2)
          r = p * u;
          new jacobiRotation(m_c, m_s, r)
        } else {
          var ps = p / q1
          val p2 = abs2(ps)
          val qs = q / q1
          val q2 = abs2(qs);

          var u = q1 * breeze.numerics.pow((p2 + q2), 0.5)
          if (p.real < 0)
            u = -u
          p1 = abs(p)
          ps = p / p1
          val m_c = p1 / u
          val m_s = -conj(ps) * (q / u)
          r = ps * u;
          new jacobiRotation(Complex(m_c, 0.0), m_s, r)
        }
    }
}
  }

object jacobiRotation {
  import Helper._

  def apply(m_c: Complex, m_s: Complex) = { new jacobiRotation(m_c, m_s) }
  def applyRotationinPlane(_x: DenseVector[Complex], _y: DenseVector[Complex], j: jacobiRotation) = {

    //    debugPrint(_x.mapValues(x => x.real - (x.real % 0.0001)))
    //   debugPrint(_y.mapValues(x => x.real - (x.real % 0.0001)))
    assert(_x.length == _y.length)
    if (j.m_c == 1 && j.m_s == 0)
      DenseMatrix.zeros[Complex](_x.size, 2)
    DenseMatrix.vertcat((DenseVector.tabulate(_x.size)(i => j.m_c * _x(i) + conj(j.m_s) * _y(i))).toDenseMatrix, (DenseVector.tabulate(_x.size)(i => -j.m_s * _x(i) + conj(j.m_c) * _y(i))).toDenseMatrix)
  }

  /*This function implements the continuous Givens rotation
     *generation algorithm found in Anderson (2000),
     *Discontinuous Plane Rotations and the Symmetric Eigenvalue Problem.
     *LAPACK Working Note 150, University of Tennessee, UT-CS-00-454, December 4, 2000. */
  def makeGivens(p: Complex, q: Complex, r1: Complex = null) = {
    var r: Complex = r1
    (p, q) match {
      case (Complex(0.0, 0.0), _) =>
        val m_c = if (p.real < 0) Complex(-1, 0) else Complex(1, 0)
        val m_s = Complex(0, 0)
        if (r != null) r = m_c * p;
        new jacobiRotation(m_c, m_s)
      case (_, Complex(0.0, 0)) =>
        val m_c = Complex(0, 0)
        val m_s = -q / abs(q)
        if (r != null) r = Complex(abs(q), 0)
        new jacobiRotation(m_c, m_s)
      case _ =>
        var p1 = norm1(p)
        var q1 = norm1(q)
        if (p1 >= q1) {
          val ps = p / p1
          val p2 = abs2(ps)
          val qs = q / p1
          val q2 = abs2(qs)

          var u = breeze.numerics.pow(Complex(1.0, 0.0) + (q2 / p2), 0.5)
          if (p.real < 0)
            u = -u;
          val m_c = 1 / u
          val m_s = -qs * conj(ps) * (m_c / p2)
          if (r != null) r = p * u;
          new jacobiRotation(m_c, m_s)
        } else {
          var ps = p / q1
          val p2 = abs2(ps)
          val qs = q / q1
          val q2 = abs2(qs);

          var u = q1 * breeze.numerics.pow((p2 + q2), 0.5)
          if (p.real < 0)
            u = -u
          p1 = abs(p)
          ps = p / p1
          val m_c = p1 / u
          val m_s = -conj(ps) * (q / u)
          if (r != null) r = ps * u;
          new jacobiRotation(Complex(m_c, 0.0), m_s)
        }
    }
}
  }