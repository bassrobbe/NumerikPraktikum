import Dependencies._

lazy val root = (project in file(".")).
  settings(
    inThisBuild(List(
      organization := "numerik_praktikum",
      scalaVersion := "2.12.4",
      version      := "0.1.0-SNAPSHOT"
    )),
    name := "SplineInterpolation",
    mainClass in Compile := Some("numerik_praktikum.SplineInterpolation"),
    libraryDependencies += scalaTest % Test
  )
