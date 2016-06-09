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


/*
This function implements the continuous Givens rotation generation algorithm found
in Anderson (2000), Discontinuous Plane Rotations and the Symmetric Eigenvalue Problem.
LAPACK Working Note 150, University of Tennessee, UT-CS-00-454, December 4, 2000.

 This class represents a Jacobi or Givens rotation.
 * This is a 2D rotation in the plane \c J of angle \f$ \theta \f$ defined by
 * its cosine \c c and sine \c s as follow:
 * \f$ J = \left ( \begin{array}{cc} c & \overline s \\ -s  & \overline c \end{array} \right ) \f$
 *
 * You can apply the respective counter-clockwise rotation to a column vector \c v by
 * applying its adjoint on the left: \f$ v = J^* v \f$
 */

object jacobi {
  implicit class IMPL_Jacobi(val M: DenseMatrix[Complex]) {

    def rotateL(jrot: jRotation) = jacobi.applyRotationinPlane(M(0, ::).t, M(1, ::).t, new jRotation(conj(jrot.c), -(jrot.s), jrot.rotation))
    def rotateR(jrot: jRotation) = jacobi.applyRotationinPlane(M(::, 0), M(::, 1), new jRotation((jrot.c), -conj(jrot.s), jrot.rotation))

    //Need to fix these...
    //def rotateL(jrot: jRotation) =  jacobi.getRotationLeft(M(0, ::).t, M(1, ::).t, new jRotation(conj(jrot.c), -(jrot.s), jrot.rotation))
    //def rotateR(jrot: jRotation) = jacobi.getRotationRight(M(::, 0), M(::, 1), new jRotation((jrot.c), -conj(jrot.s), jrot.rotation))
  }
/*
  def getRotationLeft(_x: DenseVector[Complex], _y: DenseVector[Complex], j: jRotation) = {
    val x1 = DenseVector.tabulate[Complex](_x.length) { (i) => j.c * _x(i) + conj(j.s) * _y(i) }
    val y1 = DenseVector.tabulate[Complex](_y.length) { (i) => -j.s * _x(i) + conj(j.c) * _y(i) }
    DenseMatrix.vertcat(x1.t, y1.t)
  }

  def getRotationRight(_x: DenseVector[Complex], _y: DenseVector[Complex], j: jRotation) = {
    val x1 = DenseVector.tabulate[Complex](_x.length) { (i) => j.c * _x(i) + conj(j.s) * _y(i) }
    val y1 = DenseVector.tabulate[Complex](_y.length) { (i) => -j.s * _x(i) + conj(j.c) * _y(i) }
    DenseVector.horzcat(x1, y1)
  }
*/
private  def applyRotationinPlane(_x: DenseVector[Complex], _y: DenseVector[Complex], j: jRotation) = {

    if (j.c == 1 && j.s == 0)
      DenseMatrix.zeros[Complex](_x.size, 2)

    val x1 = DenseVector.tabulate[Complex](_x.length) { (i) => j.c * _x(i) + conj(j.s) * _y(i) }
    val y1 = DenseVector.tabulate[Complex](_y.length) { (i) => -j.s * _x(i) + conj(j.c) * _y(i) }
    val res = DenseMatrix.vertcat(x1.t, y1.t)
    val res1 = DenseVector.horzcat(x1, y1)

    val i = 0

    for (i <- 0 to _x.length - 1) {
      val xi = _x(i)
      val yi = _y(i)
      _x(i) = j.c * xi + conj(j.s) * yi
      _y(i) = -j.s * xi + conj(j.c) * yi
    }
  }

/*This function implements the continuous Givens rotation
 *generation algorithm found in Anderson (2000),
 *Discontinuous Plane Rotations and the Symmetric Eigenvalue Problem.
 *LAPACK Working Note 150, University of Tennessee, UT-CS-00-454, December 4, 2000. */
/*
 * Makes *this as a Givens rotation G such that applying $ G^* $ to the left of the vector
 *  $ V = \left ( \begin{array}{c} p \\ q \end{array} \right )$ yields: $ G^* V = \left ( \begin{array}{c} r \\ 0 \end{array} \right )$.
The value of z is returned if z is not null (the default is null). Also note that G is built such that the cosine is always real.
 */

def makeGivens(p: Complex, q: Complex) = {

  (p, q) match {
    case (_, Complex(0.0, 0.0)) =>
      val c = if (p.real < 0) Complex(-1.0, 0.0) else Complex(1.0, 0.0)
      val s = Complex(0.0, 0.0)
      val r = c * p;

      new jRotation(c, s, r)
    case (Complex(0.0, 0.0), _) =>
      val c = Complex(0.0, 0.0)
      val s = -q / abs(q)
      val r = Complex(abs(q), 0.0)

      new jRotation(c, s, r)
    case _ =>
      val p1 = norm1(p)
      val q1 = norm1(q)
      if (p1 >= q1) {
	val ps = p / p1
	val p2 = abs2(ps)
	val qs = q / p1
	val q2 = abs2(qs)

	var u = pow(1.0 + (q2 / p2), 0.5)
	if (p.real < 0)
	  u = -u

	val c = Complex(1.0, 0) / u
	val s = -qs * conj(ps) * (c / p2)
	val r = p * u

	new jRotation(c, s, r)
      } else {

	val p2 = abs2(p / q1)
	val qs = q / q1
	val q2 = abs2(qs);

	var u = q1 * pow((p2 + q2), 0.5)

	if (p.real < 0)
	  u = -u

	val p1 = abs(p)
	val  ps2 = p / p1
	val c = Complex(p1 / u, 0.0)
	val s = -conj(ps2) * (q / u)
	val r = ps2 * u
	new jRotation(c, s, r)
      }
  }
}

  

  class jRotation(val c: Complex, val s: Complex, val rotation: Complex) {}
}