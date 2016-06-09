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
    }

  def main(args: Array[String]): Unit = {
    promptEnterKey();

    val M4 = DenseMatrix((1, 2, 4, 4), (5, 6, 7, 9), (9, 10, 11, 12), (13, 14, 15, 16))
    val M5 = DenseMatrix((1, 2, 4, 4, 5), (5, 6, 7, 9, 10), (9, 10, 11, 12, 12), (13, 14, 15, 16, 17), (18, 19, 20, 21, 22))
    val M3 = DenseMatrix((1, 4, 3), (4, 5, 6), (4, 3, 2))
    val M6 = DenseMatrix.rand(6, 6)

    val R = DenseMatrix((4, 1, -2, 2), (1, 2, 0, 1), (-2, 0, 3, -2), (2, 1, -2, -1))
    val I = DenseMatrix((1, 2, -3, 4), (1, 5, 0, 1), (-2, 6, 3, -2), (5, 1, -2, -1))

    val RtoC = DenseMatrix((4, 1, -2, 2), (1, 2, 0, 1), (-2, 0, 3, -2), (2, 1, -2, -1)).mapValues(Complex(_, 0.0))
    val C = DenseMatrix.tabulate[Complex](R.cols, R.rows)((i, j) => Complex(R(i, j), I(i, j)))

    //   getHessenberg(M3)
    debugPrint(M4, "MATRIX  1", 1)
    debugPrint(M4.mapActiveValues(Complex(_, 0.0)) mPow 3.21, "RESULT 1", 1)
    debugPrint(M5, "MATRIX 2", 1)
    debugPrint(M5 mPow 3.21, "RESULT 2", 1)

    debugPrint(C, "  C", 2)
    debugPrint(C mPow 2.321, "RESULT C", 2)

    debugPrint(RtoC, "  Check", 2)
    debugPrint(RtoC mPow 2.321, "RESULT RtoC", 2)

    if (bw != None) {
      bw.get.close()
    }
  }

}

