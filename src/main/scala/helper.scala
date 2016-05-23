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

  def function(s: MutableInt): Boolean = {
    s.inc() // parentheses here to denote that method has side effects
    return true
  }
  var matnum1 = MutableInt(0)
  var matnum2 = MutableInt(0)
  val showComplex = false
  val formatter = new DecimalFormat("#0.000000")
  val printEnabled = Array(false, true, false, false, true)
  val EPSILON: Double = 2.22045e-016
}

object Helper {
  import GlobalConsts._

  def abs2(n: Complex) = { (n.real * n.real) + (n.imag * n.imag) }
  def conj(n: Complex) = { Complex(n.real, -n.imag) }
  def norm1(n: Complex): Double = { abs(n.real) + abs(n.imag) }
  def biggest(M: DenseMatrix[Complex]) = norm1(sum(M(::, *)).t.reduceLeft((x, y) => if (norm1(x) > norm1(y)) x else y))

  def isMuchSmallerThan(x: Double, y: Double) = { abs(x) <= abs(y) * EPSILON }

  def printcount() = {

    function(matnum1)
    val count = matnum1.value
    println(s"************************************************************** $count **************************************************************")
  }
  def debugPrint[T: TypeTag](M: T, name: String = "", loglevel: Int = 0) =

    typeTag[T].tpe match {
      case b if b =:= typeOf[DenseMatrix[Double]] => if (printEnabled(loglevel)) {
        println(s"$name\n\n" + M.asInstanceOf[DenseMatrix[Double]].mapValues(x => "( " + formatter.format(x) + " )"))
        printcount()
      }

      case b if b =:= typeOf[DenseVector[Double]] => if (printEnabled(loglevel)) {
        println(s"$name\n\n" + M.asInstanceOf[DenseVector[Double]].mapValues(x => "( " + formatter.format(x) + " )"))
        printcount()
      }

      case b if b =:= typeOf[DenseVector[Complex]] => if (printEnabled(loglevel)) {
        if (showComplex) println(s"$name\n\n" + M.asInstanceOf[DenseVector[Complex]].mapValues(x => "( " + formatter.format(x.real) + ", " + formatter.format(x.imag) + " )"))
        else println(s"$name\n\n" + M.asInstanceOf[DenseVector[Complex]].mapValues(x => "( " + formatter.format(x.real) + " )"))
        printcount()
      }

      case b if b =:= typeOf[DenseMatrix[Complex]] => if (printEnabled(loglevel)) {
        if (showComplex) println(s"$name\n\n" + M.asInstanceOf[DenseMatrix[Complex]].mapValues(x => "( " + formatter.format(x.real) + ", " + formatter.format(x.imag) + " )"))
        else println(s"$name\n\n" + M.asInstanceOf[DenseMatrix[Complex]].mapValues(x => "( " + formatter.format(x.real) + " )"))
        printcount()
      }
      case b if b =:= typeOf[String] =>
        if (printEnabled(loglevel)) {
          println(s"$name " + M)
          printcount()

        }
      case _ => if (printEnabled(loglevel)) {

        println(s"$name " + M)
        printcount()

      }
    }
  /*
   *
   *
   *
   *    case b :DenseMatrix[Double] => println(s"$name\n\n" + b.mapValues(x => "( " + formatter.format(x) + " )"))
   case b: DenseMatrix[Complex] => println(s"$name\n\n" + b.mapValues(x => "( " + formatter.format(x.real) + ", " + formatter.format(x.imag) + " )"))
   if (printEnabled(index)) {
   function(matnum1)
   val count = matnum1.value

   }
   *
   *
   (index: Int, m: T, name: String)

   if (printEnabled(index)) {
   function(matnum1)
   val count = matnum1.value
   println(s"************************************ $count ************************************")
   m match {

   */

}
