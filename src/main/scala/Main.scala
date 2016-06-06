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

  def getHessenberg(M: DenseMatrix[Int]) =
    {
      val MD = M.mapValues(_.toDouble)
      val MC = M.mapValues(i => Complex(i.toDouble, 0.0))

      val (hessC, houseC) = MC hessenbergDecomp
      val (hessI, houseI) = M hessenbergDecomp
      val (hessD, houseD) = MD hessenbergDecomp

      val (matTC, matQC, tauC, matrixHC) = MC schurDecomp
      val (matTI, matQI, tauI, matrixHI) = M schurDecomp
      val (matTD, matQD, tauD, matrixHD) = MD schurDecomp

      debugPrint(hessC.MatrixH, "hessHC", 1)
      debugPrint(hessC.MatrixP, "hessPC", 1)
      debugPrint(houseC.tau, "tauC", 1)
      debugPrint(houseC.matrixH, "houseHC", 1)

      debugPrint(matTC, "schurTC", 1)
      debugPrint(matQC, "schurQC", 1)

      debugPrint(hessI.MatrixH, "hessHI", 1)
      debugPrint(hessI.MatrixP, "hessPI", 1)
      debugPrint(houseI.tau, "tauI", 1)
      debugPrint(houseI.matrixH, "houseHI", 1)

      debugPrint(matTI, "schurTI", 1)
      debugPrint(matQI, "schurQI", 1)

      debugPrint(hessD.MatrixH, "hessHD", 1)
      debugPrint(hessD.MatrixP, "hessPD", 1)
      debugPrint(houseD.tau, "tauD", 1)
      debugPrint(houseD.matrixH, "houseHD", 1)

      debugPrint(matTD, "schurTD", 1)
      debugPrint(matQD, "schurQD", 1)

    }

  def main(args: Array[String]): Unit = {
    // promptEnterKey();

    val M1 = DenseMatrix((1, 2, 4, 4), (5, 6, 7, 9), (9, 10, 11, 12), (13, 14, 15, 16))
    val M2 = DenseMatrix((1, 2, 4, 4, 5), (5, 6, 7, 9, 10), (9, 10, 11, 12, 12), (13, 14, 15, 16, 17), (18, 19, 20, 21, 22))
    val M3 = DenseMatrix((1, 2, 3), (4, 5, 6), (4, 3, 2))
 //   getHessenberg(M3)
    debugPrint(M3 mPow 2.321, "RESULT", 1)
  //  debugPrint(M3.mapActiveValues(Complex(_, 0.0)) mPow 2.321, "RESULT", 1)
    if (bw != None) {
      bw.get.close()
    }
  }

}

