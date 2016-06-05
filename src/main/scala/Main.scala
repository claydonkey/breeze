import breeze.linalg._
import scala.io._
import java.io.BufferedWriter
import java.io.File
import java.io.FileWriter
import java.text.DecimalFormat
import scala.Console._
import breeze.numerics._
import scala.util.control._
import breeze.math._
import scala.util.control.Breaks._
import Helper._
import GlobalConsts._
import Householder._
import Hessenberg._
import Schur._
import matrixPow._

object Main {
  def promptEnterKey(): Option[Unit] = if (Console.in.read > 10) None else promptEnterKey
  val file = new File("schur.dat")
  val bw: Option[BufferedWriter] = if (fileOutput == true) { Some(new BufferedWriter(new FileWriter(file))) } else { None }

  def main(args: Array[String]): Unit = {
    // promptEnterKey();

    val M1 = DenseMatrix((1, 2, 4, 4), (5, 6, 7, 9), (9, 10, 11, 12), (13, 14, 15, 16))
    val M2 = DenseMatrix((1, 2, 4, 4, 5), (5, 6, 7, 9, 10), (9, 10, 11, 12, 12), (13, 14, 15, 16, 17), (18, 19, 20, 21, 22))
    debugPrint(M2 mPow 2.321, "RESULT", 1)

    if (bw != None) {
      bw.get.close()
    }
  }

}

