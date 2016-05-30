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
    val COMPLEX, COMPLEX2, REAL, GRAPH, GRAPH2, GRAPH3, GRAPHCOMPLEX2, GRAPHCOMPLEX3, GRAPHCOMPLEXONLY, GRAPHCOMPLEXONLY2 = Value
  }
  def function(s: MutableInt): Boolean = {
    s.inc() // parentheses here to denote that method has side effects
    return true
  }
  var matnum1 = MutableInt(0)
  var matnum2 = MutableInt(0)
  val showComplex = false
  //val formatter = new DecimalFormat("#0.###E0")
  val formatter = new DecimalFormat("#0.##########")
  val printEnabled = Array(true, false, false, false, false, false, false)
  val EPSILON: Double = 2.22045e-016
  val currentPrintType = printType.GRAPHCOMPLEXONLY2
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
  def printcount() = {

    function(matnum1)
    val count = matnum1.value
    println(s"************************************************************** $count **************************************************************")
  }

  def printcount2(name: String) = {

    function(matnum1)
    val count = matnum1.value
    "**************************************************************" + name + " " + count + "**************************************************************\n"
  }
  def adder(x: Double) = { if (x > 0) { "+" } else { "" } };
  def debugPrint[T: TypeTag](M: T, name: String = "", loglevel: Int = 0) =

    typeTag[T].tpe match {
      case b if b =:= typeOf[DenseMatrix[Double]] => if (printEnabled(loglevel)) {
        currentPrintType match {

          case printType.GRAPH =>

            println(s"$name\n\n" + M.asInstanceOf[DenseMatrix[Double]].mapValues {
              var i = 0; (x) =>
                i += 1;
                i + " " + formatter.format(x) + ", "
            })

          case printType.GRAPH2 => Main.bw.get.write(("" + M.asInstanceOf[DenseMatrix[Double]].mapValues { (x) => formatter.format(x) }).filter(_ >= ' ') + "\n")
          case printType.GRAPHCOMPLEX2 => Main.bw.get.write(("" + M.asInstanceOf[DenseMatrix[Double]].mapValues { (x) => formatter.format(x) }).filter(_ >= ' ') + "\n")
          case printType.GRAPH3 => Main.bw.get.write(printcount2(name) + M.asInstanceOf[DenseMatrix[Double]].mapValues { (x) => "( " + formatter.format(x) + " )" } + "\n")
          case printType.GRAPHCOMPLEX3 => Main.bw.get.write(printcount2(name) + M.asInstanceOf[DenseMatrix[Double]].mapValues { (x) => "( " + formatter.format(x) + " )" } + "\n")
          case printType.REAL => println(printcount2(name) + M.asInstanceOf[DenseMatrix[Double]].mapValues(x => "( " + formatter.format(x) + " )"))

          case printType.COMPLEX2 => println((printcount2(name) + M.asInstanceOf[DenseMatrix[Double]].mapValues { (x) => formatter.format(x) }))
          case _ => println(printcount2(name) + M.asInstanceOf[DenseMatrix[Double]].mapValues(x => "( " + formatter.format(x) + " )"))
        }

      }
      case b if b =:= typeOf[DenseVector[Double]] => if (printEnabled(loglevel)) {
        currentPrintType match {
          case printType.GRAPH2 => Main.bw.get.write(("" + M.asInstanceOf[DenseVector[Double]].mapValues { (x) => formatter.format(x) }).filter(_ >= ' ') + "\n")
          case printType.GRAPH3 => Main.bw.get.write(printcount2(name) + M.asInstanceOf[DenseVector[Double]].mapValues { (x) => "( " + formatter.format(x) + " )" } + "\n")
          case printType.REAL => println(printcount2(name) + M.asInstanceOf[DenseVector[Double]].mapValues(x => "( " + formatter.format(x) + " )"))
          case _ => println(printcount2(name) + M.asInstanceOf[DenseVector[Double]].mapValues(x => "( " + formatter.format(x) + " )"))

        }
      }

      case b if b =:= typeOf[DenseVector[Complex]] => if (printEnabled(loglevel)) {
        currentPrintType match {
          case printType.GRAPH2 => Main.bw.get.write(("" + M.asInstanceOf[DenseVector[Complex]].mapValues { (x) => formatter.format(x.real) }).filter(_ >= ' ') + "\n")
          case printType.GRAPH3 => Main.bw.get.write(printcount2(name) + M.asInstanceOf[DenseVector[Complex]].mapValues { (x) => "( " + formatter.format(x.real) + " )" } + "\n")
          case printType.GRAPHCOMPLEX3 => Main.bw.get.write(printcount2(name) + M.asInstanceOf[DenseVector[Complex]].mapValues { (x) => "( " + formatter.format(x.real) + ", " + formatter.format(x.imag) + " )" } + "\n")
          case printType.GRAPHCOMPLEXONLY => Main.bw.get.write(printcount2(name) + M.asInstanceOf[DenseVector[Complex]].mapValues { (x) => "( " + formatter.format(x.imag) + " )" } + "\n")
          case printType.REAL => println(printcount2(name) + M.asInstanceOf[DenseVector[Complex]].mapValues(x => "( " + formatter.format(x.real) + " )"))
          case _ => println(printcount2(name) + M.asInstanceOf[DenseVector[Complex]].mapValues(x => "( " + formatter.format(x.real) + ", " + formatter.format(x.imag) + " )"))
        }
      }

      case b if b =:= typeOf[DenseMatrix[Complex]] => if (printEnabled(loglevel)) {
        currentPrintType match {
          case printType.GRAPH2 => Main.bw.get.write(("" + M.asInstanceOf[DenseMatrix[Complex]].mapValues { (x) => formatter.format(x.real) }).filter(_ >= ' ') + "\n")
          case printType.GRAPH3 => Main.bw.get.write(printcount2(name) + M.asInstanceOf[DenseMatrix[Complex]].mapValues { (x) => "( " + formatter.format(x.real) + " )" } + "\n")
          case printType.GRAPHCOMPLEX3 => Main.bw.get.write(printcount2(name) + M.asInstanceOf[DenseMatrix[Complex]].mapValues { (x) => "( " + formatter.format(x.real) + ", " + formatter.format(x.imag) + " )" } + "\n")
          case printType.GRAPHCOMPLEXONLY => Main.bw.get.write(printcount2(name) + M.asInstanceOf[DenseMatrix[Complex]].mapValues { (x) => "( " + formatter.format(x.imag) + " )" } + "\n")
          case printType.GRAPHCOMPLEXONLY2 => Main.bw.get.write(("" + M.asInstanceOf[DenseMatrix[Complex]].mapValues { (x) => formatter.format(x.imag) }).filter(_ >= ' ') + "\n")
          case printType.REAL => println(printcount2(name) + M.asInstanceOf[DenseMatrix[Complex]].mapValues(x => "( " + formatter.format(x.real) + " )"))
          case _ => println(printcount2(name) + M.asInstanceOf[DenseMatrix[Complex]].mapValues(x => "( " + formatter.format(x.real) + ", " + formatter.format(x.imag) + " )"))
        }
      }

      case b if b =:= typeOf[Array[Double]] => if (printEnabled(loglevel)) {
        currentPrintType match {
          case _ => Main.bw.get.write(printcount2(name) + M.asInstanceOf[Array[Double]].deep.mkString("\n") + "\n")

        }

      }

      case _ => if (printEnabled(loglevel)) {
        currentPrintType match {
          case printType.GRAPH2 => Main.bw.get.write(M + "\n")
          case printType.GRAPH3 => Main.bw.get.write(printcount2(name) + M + "\n")
          case printType.GRAPHCOMPLEX3 => Main.bw.get.write(printcount2(name) + M + "\n")
          case printType.GRAPHCOMPLEX2 => Main.bw.get.write(M + "\n")
          case _ => println(printcount2(name) + M)
        }
        //   printcount()

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
