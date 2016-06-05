import breeze.linalg._
import scala.io._
import java.text.DecimalFormat
import scala.Console._
import breeze.numerics._
import scala.util.control._
import breeze.math._
import reflect.runtime.universe._
import scala.util.control.Breaks._

object GlobalConsts {
  implicit class MutableInt(var value: Int) {
    def inc() = { value += 1 }
  }
  object printType extends Enumeration {
    type printType = Value
    val REAL, IMAG, BOTH = Value
  }
  def function(s: MutableInt): Boolean = {
    s.inc() // parentheses here to denote that method has side effects
    true
  }

  var showHouseholder = false
  var matnum1 = MutableInt(0)
  var matnum2 = MutableInt(0)
  val showComplex = false
  val showTitles = true
  val showLines = false
  val fileOutput = true
  val formatter = new DecimalFormat("#0.####E0")
  //val formatter = new DecimalFormat("#0.######")
  val printEnabled = Array(true, true, true, true, true, false, showHouseholder, true) //0,1 for debugging last part
  val EPSILON: Double = 2.22045e-016
  val currentPrintType = printType.BOTH
}

object Helper {

  import GlobalConsts._
  def adj(M: DenseMatrix[Complex]): DenseMatrix[Complex] = { (det(M.mapValues(_.real)) * inv(M.mapValues(_.real))).mapValues(Complex(_, 0.0)) } //  adjoint of vector is vector itself??
  //def adj(M: DenseMatrix[Complex]) = { inv(M) * det(M) }  //  adjoint of vector is vector itself??
  def abs2(n: Complex): Double = { (n.real * n.real) + (n.imag * n.imag) }
  def conj(n: Complex) = { Complex(n.real, -n.imag) }
  def norm1(n: Complex): Double = { abs(n.real) + abs(n.imag) }
  def biggest(M: DenseMatrix[Complex]) = norm1(sum(M(::, *)).t.reduceLeft((x, y) => if (norm1(x) > norm1(y)) x else y))
  val M_PI = 3.14159265358979323846
  def sinh2(c: Complex) = { Complex(sinh(c.real) * cos(c.imag), cosh(c.real) * sin(c.imag)) }
  def isMuchSmallerThan(x: Double, y: Double) = { abs(x) <= abs(y) * EPSILON }
  implicit def Y3[A1, A2, A3, B](f: ((A1, A2, A3) => B) => ((A1, A2, A3) => B)): (A1, A2, A3) => B = f(Y3(f))(_, _, _)
  type fType[A] = (Int, DenseMatrix[Complex]) => A
  def Y[A, B](f: (A => B) => (A => B)): A => B = f(Y(f))(_)
  def atan2h(x: Complex, y: Complex): Complex = { val z = x / y; if ((y == 0) || (abs2(z) > pow(GlobalConsts.EPSILON, 0.5))) (0.5) * log((y + x) / (y - x)) else z + z * z * z / 3 }
  //implicit def Y1[Int, DenseMatrix[Complex] ,  DenseMatrix[Complex]](f: (fType[DenseMatrix[Complex]]) => ((Int, DenseMatrix[Complex]) => A)): (Int, DenseMatrix[Complex]) => A = f(Y1(f))(_, _)

  def printcount2(name: String) = {

    function(matnum1)
    val count = matnum1.value
    "************************************************************** " + count + " *** " + name + "**************************************************************\n"
  }
  implicit def enrichString2(stuff: String) =
    new {
      def showTitle(name: String) = { if (showTitles) printcount2(name) + stuff else stuff }
    }

  implicit def enrichString(stuff: String) =
    new {
      def oneLiner = { if (showLines) stuff.filter(_ >= ' ') else stuff }
    }

  def output(str: String) = { if (fileOutput) Main.bw.get.write(str) else print(str) }
  def adder(x: Double) = { if (x > 0) { "+" } else { "" } }
  def debugPrint[T: TypeTag](M: T, name: String = "", loglevel: Int = 0) =

    typeTag[T].tpe match {
      case b if b =:= typeOf[DenseMatrix[Double]] => if (printEnabled(loglevel)) {
        currentPrintType match {

          case _ => output(("" + M.asInstanceOf[DenseMatrix[Double]].mapValues { (x) => formatter.format(x) }).oneLiner.showTitle(name) + "\n")
        }

      }
      case b if b =:= typeOf[DenseVector[Double]] => if (printEnabled(loglevel)) {
        currentPrintType match {

          case _ => output(("" + M.asInstanceOf[DenseVector[Double]].mapValues { (x) => formatter.format(x) }).oneLiner.showTitle(name) + "\n")
        }
      }

      case b if b =:= typeOf[DenseVector[Complex]] => if (printEnabled(loglevel)) {
        currentPrintType match {
          case printType.REAL => output(("" + M.asInstanceOf[DenseVector[Complex]].mapValues { (x) => formatter.format(x.real) }).oneLiner.showTitle(name) + "\n")
          case printType.IMAG => output(("" + M.asInstanceOf[DenseVector[Complex]].mapValues { (x) => formatter.format(x.imag) }).oneLiner.showTitle(name) + "\n")
          case _ => output(("" + M.asInstanceOf[DenseVector[Complex]].mapValues { (x) => formatter.format(x.real) + "," + formatter.format(x.imag) }).oneLiner.showTitle(name) + "\n")
        }
      }
      case b if b =:= typeOf[DenseMatrix[Complex]] => if (printEnabled(loglevel)) {
        currentPrintType match {
          case printType.REAL => output(("" + M.asInstanceOf[DenseMatrix[Complex]].mapValues { (x) => formatter.format(x.real) }).oneLiner.showTitle(name) + "\n")
          case printType.IMAG => output(("" + M.asInstanceOf[DenseMatrix[Complex]].mapValues { (x) => formatter.format(x.imag) }).oneLiner.showTitle(name) + "\n")
          case _ => output(("" + M.asInstanceOf[DenseMatrix[Complex]].mapValues { (x) => "(" + formatter.format(x.real) + "," + formatter.format(x.imag) + ")" }).oneLiner.showTitle(name) + "\n")
        }
      }

      case b if b =:= typeOf[Array[Double]] => if (printEnabled(loglevel)) {
        currentPrintType match {
          case _ => output("" + M.asInstanceOf[Array[Double]].deep.mkString("\n").oneLiner.showTitle(name) + "\n")
        }
      }

      case b if b =:= typeOf[Array[DenseVector[Complex]]] => if (printEnabled(loglevel)) {
        currentPrintType match {
          case _ => output("" + M.asInstanceOf[Array[DenseVector[Complex]]].deep.mkString("\n").oneLiner.showTitle(name) + "\n")
        }
      }
      case _ => if (printEnabled(loglevel)) {
        currentPrintType match {
          case printType.REAL => output(M.toString.oneLiner.showTitle(name) + "\n")
          case printType.IMAG => output(M.toString.oneLiner.showTitle(name) + "\n")
          case _ => output(M.toString.oneLiner.showTitle(name) + "\n")
        }
      }
    }
}
