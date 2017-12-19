package numerik_praktikum

import Math._

object SplineInterpolation {

  type Tuple3 = (Double, Double, Double)
  type Tuple4 = (Double, Double, Double, Double)
  type Tuple5 = (Double, Double, Double, Double, Double)

  def flatten1(t: ((Int, Double, Double), Double)) : (Int, Double, Double, Double) = (t._1._1, t._1._2, t._1._3, t._2 )
  def flatten2(t: (Tuple3, Double)) : Tuple4 = (t._1._1, t._1._2, t._1._3, t._2 )

  val s : (Tuple5, Double) => Double = (c : Tuple5, t : Double) =>
    c._2 + c._3*(t-c._1) + c._4*(t-c._1)*(t-c._1)/2 + c._5*(t-c._1)*(t-c._1)*(t-c._1)/6

  val ds : (Tuple5, Double) => Double = (c : Tuple5, t : Double) => c._3 + c._4*(t-c._1) + c._5*(t-c._1)*(t-c._1)/2

  val dds : (Tuple5, Double) => Double = (c : Tuple5, t : Double) => c._4 + c._5*(t-c._1)

  val S = List(s, ds, dds)

  val testValues : List[(Double => Double, Double => Double, Double => Double, Double, Double, String)] = List(
    ((t : Double) => sin(t), (t : Double) => cos(t), (t : Double) => -sin(t), -PI, PI, "sin(t)"),
    ((t : Double) => cos(t), (t : Double) => -sin(t), (t : Double) => -cos(t), -PI, PI, "cos(t)"),
//    ((t : Double) => abs(t)*sin(t)*sin(t)/3, (t*sin(t)*sin(t)+2*t*t*cos(t))/3/abs(s), -PI, PI, "1/3*abs(t)*sin^2(t)"),
    ((t : Double) => sqrt(t), (t : Double) => 1/(2*sqrt(t)), (t : Double) => -pow(t, -3/2)/4, 1, 2, "t^(1/2)"),
    ((t : Double) => exp(t), (t : Double) => exp(t), (t : Double) => exp(t), -1, 1, "exp(t)"),
    ((t : Double) => 1f/(1+t*t), (t : Double) => -2*t/pow(t*t+1,2), (t : Double) => (6*t*t-2)/pow(t*t+1,3), -1, 1, "1/(1+t^2)"))

  def main(args : Array[String]) {

    def p_kToString(p : List[Tuple3]) : String = p map { case (h1, h2, p) => s"\t\th1 = ${h1}   h2 = ${h2}   p_k = ${p}\n"} reduce (_ + _)

    def printp(E : List[List[(Double, Double)]]) : Unit =

      (0 to 2 toList) foreach {k => println(s"--------------- ${k}-te Ableitung (k = ${k})\n${p_kToString(getp(E(k)))}")}

    def printf(in : (Double => Double, Double => Double, Double => Double, Double, Double, String)) : Unit = {

      val funcs = List((in._1,s), (in._2,ds), (in._3,dds))
      val h = (0 to 4 toList) map { i => (in._5-in._4) / pow(2, i) }

      println(s"\n===== f(t) = ${in._6}\n\n----------Interpolation Typ (a) - natuerlicher Spline\n")
      printp(funcs map { case (f,g) => h map { h_i => (getE(f, g, interpAC(in._1, in._4, in._5, h_i), in._4, in._5, h_i), h_i) } })

      println(s" ---------- Interpolation Typ (b) - eingespannter Spline\n")
      printp(funcs map { case (f,g) => h map { h_i =>
        (getE(f, g, interpB(in._1, in._4, in._5, h_i, in._2(in._4), in._2(in._5)), in._4, in._5, h_i), h_i) } })

      println(s"---------- Interpolation Typ (c)\n")
      printp(funcs map { case (f,g) => h map { h_i =>
        (getE(f, g, interpAC(in._1, in._4, in._5, h_i, in._3(in._4), in._3(in._5)), in._4, in._5, h_i), h_i) } })
    }

    //testValues foreach printf

    //only for testing purposes
    //show values step by step for f(t) = sin(t) with interpolation type (a)

    val coeffs = interpAC((t : Double) => sin(t), -PI, PI, PI/2)

    println(s"h = PI/2   sin(PI/6) ~= ${s(coeffs(2), PI/6)}")

    val E = getE((t : Double) => sin(t), s, coeffs, -PI, PI, PI/2)

    println(s"E_PI/2,0 = ${E}")

  }


  private def interpAC(f : Double => Double, a : Double, b : Double, h : Double, c_0 : Double = 0d, c_n : Double = 0d) : List[Tuple5] = {

    val values = (0 to ((b-a)/h).round.toInt).toList map { i => (i, a+i*h, f(a+i*h)) }
    val p = values.map{x => (x._1, x._3)}.toMap

    val coeffList = List.fill(p.size-2)(h/6,2*h/3,h/6) zip ((1 to p.size-2 toList) map { i => (p(i-1)+p(i+1)-2*p(i))/h })

    values zip (c_0 :: thomasAlgorithm(coeffList map flatten2)) :+ c_n map flatten1 sliding 2 map { case c :: n :: Nil =>
      (c._2, c._3, (n._3-c._3) / h - h * (2 * c._4 + n._4) / 6, c._4, (n._4 - c._4) / h) } toList
  }

  private def interpB(f : Double => Double, a : Double, b : Double, h : Double, df_0 : Double, df_n : Double) : List[Tuple5] = {

    val values = (0 to ((b-a)/h).round.toInt).toList map { i => (i, a+i*h, f(a+i*h)) }
    val p = values.map{x => (x._1, x._3)}.toMap

    val coeffList = ((0d, h/3, h/6) :: List.fill(p.size-2)(h/6,2*h/3,h/6)) :+ (h/6, h/3, 0d) zip
      ((p(1)-p(0))/h-df_0 :: ((1 to p.size-2 toList) map { i => (p(i-1)+p(i+1)-2*p(i))/h })) :+ df_n-(p(p.size-1)-p(p.size-2))/h

    values zip thomasAlgorithm(coeffList map flatten2) map flatten1 sliding 2 map { case c :: n :: Nil =>
      (c._2, c._3, (n._3-c._3) / h - h * (2 * c._4 + n._4) / 6, c._4, (n._4 - c._4) / h) } toList
  }

  private def getp(E : List[(Double, Double)]) : List[Tuple3] =

    (0 until E.length foldLeft List[Tuple3]()) { (p, i) => p ::: (E.drop(i + 1) map { e =>
      (E(i)._2, e._2, (log(E(i)._1) - log(e._1)) / (log(E(i)._2) - log(e._2))) }) }

  private def getE(f : Double => Double, g : (Tuple5, Double) => Double, coeffs : List[Tuple5], a : Double, b : Double, h : Double) : Double =

    (a to b by h/10 toList) map { tau => abs(g(coeffs(findInterval(coeffs map (_._1), tau)),tau) - f(tau)) } reduce (_ max _)

  private def thomasAlgorithm(A : List[Tuple4]) : List[Double] = {

    def forwardSub(i : Int, A : List[Tuple4]) : List[Tuple4] = { i match {

      case n if n == A.length => A

      case 0 => forwardSub(i+1, A updated(0, (0.0, A(0)._2, A(0)._3 / A(0)._2, A(0)._4 / A(0)._2)))

      case _ => forwardSub(i+1, A updated(i, (A(i)._1, A(i)._2, A(i)._3/(A(i)._2-A(i-1)._3*A(i)._1),
        (A(i)._4-A(i-1)._4*A(i)._1)/(A(i)._2-A(i-1)._3*A(i)._1))))
    }}

    def backwardSub(i : Int, A : List[Tuple4], x : List[Double] = List()) : List[Double] = { i match {

      case -1 => x

      case n if n == A.length-1 => backwardSub(i-1, A, A(i)._4 :: x)

      case _ => backwardSub(i-1, A, A(i)._4 - A(i)._3 * x(0) :: x)

    }}
    backwardSub(A.length-1, forwardSub(0, A))
  }

  def findInterval(borders : List[Double], t : Double) : Int =

    borders.indexOf(borders reduce { (x : Double, y : Double) => if ((x max y) < t) x max y else x min y })
}