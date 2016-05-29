import sbt._
import sbt.Keys._
import com.typesafe.sbt.SbtScalariform
import com.typesafe.sbt.SbtScalariform.ScalariformKeys

/**
 * Please use plain text editor to edit this file instead of NetBeans (To be supported)
 */
object Build extends sbt.Build {

  lazy val root = Project("householder", file("."))
    .settings(basicSettings: _*)
    .settings(libraryDependencies ++= Dependencies.basic)

  lazy val basicSettings = Seq(
    organization := "claydonkey.com",
    version := "0.1.0",
    scalaVersion := "2.11.8",
    scalacOptions ++= Seq("-unchecked", "-deprecation"),
    javacOptions ++= Seq("-source", "1.6", "-target", "1.6"),
    resolvers ++= Seq(
   "My Scala Repo" at "http://think-station:8081/nexus/content/groups/public/",
      "My Scala Snapshots" at "http://think-station:8081/nexus/content/repositories/snapshots/",
      Resolver.mavenLocal,
      Resolver.sonatypeRepo("snapshots"),
      Resolver.sonatypeRepo("releases"),
      Resolver.typesafeRepo("releases")),
    initialCommands in console := """
	|import breeze._
	|import breeze.linalg._
	|import breeze.signal._
	|import com.quantifind.charts.Highcharts._
	|import Math._
	|import scala.util.Random
	|import breeze.math._
	""".stripMargin

  )

  // scalariform code format settings
  SbtScalariform.scalariformSettings // enable scalariformSettings
  import scalariform.formatter.preferences._
  ScalariformKeys.preferences := ScalariformKeys.preferences.value
    .setPreference(RewriteArrowSymbols, false)
    .setPreference(AlignParameters, true)
    .setPreference(AlignSingleLineCaseStatements, true)
    .setPreference(DoubleIndentClassDeclaration, true)
    .setPreference(IndentSpaces, 2)
}



object Dependencies {
  // ---- define dependencies libs
  var basic: Seq[ModuleID] = Seq( "org.scalanlp" %% "breeze" % "0.13-SNAPSHOT" changing () withJavadoc (),
    "org.scalanlp" %% "breeze-natives" % "0.13-SNAPSHOT",
    "org.scalanlp" %% "breeze-viz" % "0.13-SNAPSHOT",
    "com.quantifind" %% "wisp" % "0.0.4",
  "org.spire-math" %% "spire" % "0.11.0"
  )
}
