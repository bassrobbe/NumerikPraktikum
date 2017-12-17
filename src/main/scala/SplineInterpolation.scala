package numerik_praktikum

import Math._

object SplineInterpolation {

  type Tuple5 = (Double, Double, Double, Double, Double)

  def flatten1(t: ((Int, Double, Double), Double)) : (Int, Double, Double, Double) = (t._1._1, t._1._2, t._1._3, t._2 )
  def flatten2(t: ((Double, Double, Double), Double)) : (Double, Double, Double, Double) = (t._1._1, t._1._2, t._1._3, t._2 )

  def tuple(t: List[Double]) : (Double, Double, Double) = t match { case a :: b ::c :: Nil => (a,b,c) }

  val s : (Tuple5, Double) => Double = (c : Tuple5, t : Double) =>
    c._2 + c._3*(t-c._1) + c._4*(t-c._1)*(t-c._1)/2 + c._5*(t-c._1)*(t-c._1)*(t-c._1)/6

  val ds : (Tuple5, Double) => Double = (c : Tuple5, t : Double) => c._3 + c._4*(t-c._1) + c._5*(t-c._1)*(t-c._1)/2

  val dds : (Tuple5, Double) => Double = (c : Tuple5, t : Double) => c._4 + c._5*(t-c._1)

  val S = List(s, ds, dds)

  //val epsilon = 0.001

  val testValues : List[(Double => Double,Double => Double,Double => Double, Double, Double, String)] = List(
    ((t : Double) => sin(t), (t : Double) => cos(t), (t : Double) => -sin(t), -PI, PI, "sin"))//,
//    ((t : Double) => cos(t), (t : Double) => -sin(t), (t : Double) => -cos(t), -PI, PI)),
//    ((t : Double) => 1d/3*abs(t)*sin(t)*sin(t), -PI, PI),
//    ((t : Double) => sqrt(t), 1, 2),
//    ((t : Double) => exp(t), -1, 1),
//    ((t : Double) => 1f/(1+t*t), -1, 1))

  def main(args : Array[String]) {

    println("Interpolation Typ (a) - natuerlicher Spline\n")

    val in = testValues(0)

    val a = interpolationAC(in._1, in._4, in._5, PI/2)
    val E2 = getE(in._1, s, a, in._4, in._5, PI/2)

    println(E2)





    val E1 = (0 to 3 toList) map { i => (in._5-in._4) / pow(2, i) } map
      { h => {/*println(h);*/(getE(in._1, s, interpolationAC(in._1, in._4, in._5, h), in._4, in._5, h), h) }}
    println(E1)

    val E = List((in._1,s), (in._2,ds), (in._3,dds)) map { case (f,g) => (0 to 4 toList) map { i =>
      (in._5-in._4) / pow(2, i) } map { h => (getE(f, g, interpolationAC(in._1, in._4, in._5, h), in._4, in._5, h), h) }}


    println(E)
    println

   println(s"f(t) = ${in._6}(t) => p_1 = ${getp(E(0))}   p_2 = ${getp(E(1))}   p_3 = ${getp(E(2))}")


  }


  private def interpolationAC(f : Double => Double, a : Double, b : Double, h : Double, c_0 : Double = 0d, c_n : Double = 0d) : List[Tuple5] = {

    val values = (0 to ((b-a)/h).round.toInt).toList map { i => (i, a+i*h, f(a+i*h)) }
    val p = values.map{x => (x._1, x._3)}.toMap

    val coeffList = List.fill(p.size-2)(1d/6,2d/3,1d/6) zip ((1 to p.size-2 toList) map { i => (p(i-1)+p(i+1)-2*p(i))/h })

    values zip (c_0 :: thomasAlgorithm(coeffList map flatten2)) :+ c_n map flatten1 sliding 2 map { case c :: n :: Nil =>
      (c._2, c._3, (n._3-c._3) / h - h * (2 * c._4 + n._4) / 6, c._4, (n._4 - c._4) / h) } toList
  }

  private def getp(E : List[(Double, Double)]) : Double ={

    println((0 until E.length foldLeft List[Double]()) { (p, i) => p ::: (E.drop(i + 1) map { e =>
      (log(E(i)._1) - log(e._1)) / (log(E(i)._2) - log(e._2)) })})

    ((0 until E.length foldLeft List[Double]()) { (p, i) => p ::: (E.drop(i + 1) map { e =>
      (log(E(i)._1) - log(e._1)) / (log(E(i)._2) - log(e._2)) })} reduce (_ + _)) / ((E.length * (E.length - 1)) / 2)}

  def getE(f : Double => Double, g : (Tuple5, Double) => Double, coeffs : List[Tuple5], a : Double, b : Double, h : Double) : Double =

    (a to b by h/10 toList) map { tau => abs(g(coeffs(findInterval(coeffs map (_._1), tau)),tau) - f(tau)) } reduce (_ max _)

  private def thomasAlgorithm(A : List[(Double,Double,Double,Double)]) : List[Double] = {

    def forwardSub(i : Int, A : List[(Double,Double,Double,Double)]) : List[(Double,Double,Double,Double)] = { i match {

      case n if n == A.length => A

      case 0 => forwardSub(i+1, A updated(0, (0.0, A(0)._2, A(0)._3 / A(0)._2, A(0)._4 / A(0)._2)))

      case _ => forwardSub(i+1, A updated(i, (A(i)._1, A(i)._2, A(i)._3/(A(i)._2-A(i-1)._3*A(i)._1),
        (A(i)._4-A(i-1)._4*A(i)._1)/(A(i)._2-A(i-1)._3*A(i)._1))))
    }}

    def backwardSub(i : Int, A : List[(Double,Double,Double,Double)], x : List[Double] = List()) : List[Double] = { i match {

      case -1 => x

      case n if n == A.length-1 => backwardSub(i-1, A, A(i)._4 :: x)

      case _ => backwardSub(i-1, A, A(i)._4 - A(i)._3 * x(0) :: x)

    }}
    backwardSub(A.length-1, forwardSub(0, A))
  }


  private def gaussWithPivot(A : List[List[Double]]) : List[Double] = {

    def forwardSub(A : List[List[Double]], next : (Int,Int)): List[List[Double]] = next match {

      case (i,j) if i==j && i==A.length-1 => A

      case (i,j) if i==j =>
        val pivot = (j + 1 until A.length foldLeft j) { (pivot,i) => if (A(i)(j).abs > A(pivot)(j).abs) i else pivot }
        forwardSub(A.updated(pivot, A(j)).updated(j, A(pivot)), (i+1,j))

      case (i,j) => val newRow = (0 to A.length toList) map { _ match {
            case k if k <= j => 0
            case k => A(i)(k) - A(j)(k)*A(i)(j)/A(j)(j)
        }}
        if (i+1<A.length) forwardSub(A.updated(i, newRow), (i+1,j)) else forwardSub(A.updated(i, newRow), (j+1,j+1))
    }

    def backwardSub(A : List[List[Double]], x : List[Double], next : (Int,Int)) : List[Double] = next match {

      case (i,j) if j==A.length => if (i==0) x else backwardSub(A, x, (i-1,i-1))

      case (i,j) if i==j => backwardSub(A, x.updated(i, x(i) / A(i)(i)), (i,j+1))

      case (i,j) => backwardSub(A, x.updated(i, x(i) - A(i)(j)*x(j)/A(i)(i)), (i,j+1))
    }

    val A1 = forwardSub(A, (0,0))
    backwardSub(A1, A1 map {_.last}, (A.length-1,A.length-1))
  }

  def findInterval(borders : List[Double], t : Double) : Int =

    borders.indexOf(borders reduce { (x : Double, y : Double) => if ((x max y) < t) x max y else x min y })
}