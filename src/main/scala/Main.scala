import breeze.linalg._
import scala.io._
import java.text.DecimalFormat
import scala.Console._
import breeze.numerics._
import scala.util.control._
import breeze.math._
import scala.util.control.Breaks._
import Helper._

object Main {
  def promptEnterKey(): Option[Unit] = if (Console.in.read > 10) None else promptEnterKey

  def applySchur() = {
    val formatter = new DecimalFormat("#0.0000")
    val M = DenseMatrix((1, 2, 4, 4), (5, 6, 7, 9), (9, 10, 11, 12), (13, 14, 15, 16)).mapValues(Complex(_, 0.0))
    val mySchur = complexSchur(M)
    val matT = mySchur.matT
   val matQ = mySchur.matQ

   debugPrint( matT, "mySchur.matT",1 )

 debugPrint( matQ, "mySchur.matQ" ,1)
  }

  def main(args: Array[String]): Unit = {
    promptEnterKey();
       applySchur() ;
    //matrixPow.fract(3.43, DenseMatrix((1, 2, 4, 4), (5, 6, 7, 9), (9, 10, 11, 12), (13, 14, 15, 16)).mapValues(Complex(_, 0.0)))



  }
}



