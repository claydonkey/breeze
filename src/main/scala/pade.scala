import breeze.linalg._
import scala.io._
import scala.Console._
import breeze.numerics._
import scala.util.control._
import breeze.math._
import scala.util.control.Breaks._
import Helper._

object padePower {

  def degree = (normIminusT: Double) => {
    val maxNormForPade = Array(2.8064004e-1f /* degree = 3 */ , 4.3386528e-1f)
    var degree = 3
    for (degree <- 3 to 4)
      if (normIminusT <= maxNormForPade(degree - 3))
        degree
    degree
  }
  def apply(IminusT: DenseMatrix[Complex], m_p: Double) = {

    val _degree = degree(m_p)
    val i = _degree << 1
    val r = IminusT.map(_ * (m_p - _degree.toDouble) / ((i - 1) << 1))
    val index = 0
    val M = upperTriangular(DenseMatrix.eye[Complex](IminusT.rows)) + r
    val T1 = -1.0 * m_p
    val T = IminusT.mapValues(_.real) * T1 // BIG PROBLEMMMO
    debugPrint(s"IminusT:\n$IminusT\nT:\n$T\nT1:\n$T1\nM:\n$M\n")

    val IminusR = IminusT.mapValues(_.real)
    debugPrint(s"IminusR:\n$IminusR\n")
    // else IminusT.map(_ * (m_p - (t >> 1)) / ((t - 1) << 1))
    val res = T \ M.mapValues(_.real) // BIG PROBLEMMMO
    res + DenseMatrix.eye[Double](IminusT.rows)
    res.mapValues(Complex(_, 0.0))
  }
}
