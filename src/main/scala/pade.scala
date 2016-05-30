import breeze.linalg._
import scala.io._
import scala.Console._
import breeze.numerics._
import scala.util.control._
import breeze.math._
import scala.util.control.Breaks._
import Helper._
import HouseholderQR._

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
    val res = IminusT.map(_ * (m_p - _degree.toDouble) / ((i - 1) << 1))
    val index = 0
    val M :DenseMatrix[Complex]= DenseMatrix.tabulate[Complex](res.rows, res.cols){(x, y) => if (x == y) Complex(1.0, 0) else res(x, y)}
    val T1 = -1.0 * Complex(m_p, 0)
    val T = IminusT * T1
   (M.mapValues(_.real) \ T.mapValues(_.real)).mapValues(Complex(_, 0.0)):+ DenseMatrix.eye[Complex](IminusT.rows)   // BIG PROBLEMMMO


  }
}
