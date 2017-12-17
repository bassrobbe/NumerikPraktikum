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

  //val epsilon = 0.001

  val testValues : List[(Double => Double,Double => Double,Double => Double, Double, Double)] = List(
    ((t : Double) => sin(t), (t : Double) => cos(t), (t : Double) => -sin(t), -PI, PI))//,
//    ((t : Double) => cos(t), (t : Double) => -sin(t), (t : Double) => -cos(t), -PI, PI)),
//    ((t : Double) => 1d/3*abs(t)*sin(t)*sin(t), -PI, PI),
//    ((t : Double) => sqrt(t), 1, 2),
//    ((t : Double) => exp(t), -1, 1),
//    ((t : Double) => 1f/(1+t*t), -1, 1))

  def main(args : Array[String]) {

    println("Interpolation Typ (a) - natuerlicher Spline\n")

    val in = testValues(0)
    val E = getE(0,in._5-in._4, in)

  }


  private def interpolationA(f : Double => Double, a : Double, b : Double, h : Double) : List[Tuple5] = {

    val values = (0 to ((b-a)/h).round.toInt).toList map { i => (i, a+i*h, f(a+i*h)) }
    val p = values.map{x => (x._1, x._3)}.toMap

    val coeffList = List.fill(p.size-2)(1d/6,2d/3,1d/6) zip ((1 to p.size-2 toList) map { i => (p(i-1)+p(i+1)-2*p(i))/h })

    values zip (0d :: thomasAlgorithm(coeffList map flatten2)) :+ 0d map flatten1 sliding(2) map { case c :: n :: Nil =>
      (c._2, c._3, (n._3-c._3)/h-h*(2*c._4+n._4)/6, c._4, (n._4-c._4)/h)} toList
  }

  private def getp(E : List[(Double, Double)]) : List[(Double, Double, Double)] = {
//
//    def getNext(E : List[(Double, Double)], p : List[(Double, Double, Double)] = List()) : List[(Double, Double, Double)] = {
//      if (! E.isEmpty) getNext(E.tail, p ::: E.tail.map { e => (E.head._2, e._2, (log(E.head._1)-log(e._1))/(log(E.head._2)-log(e._2))) })
//      else p
//    }

    (0 until E.length foldLeft List[(Double, Double, Double)]()) { (p, i) => p ::: (E.drop(i + 1) map
      { e => (E(i)._2, e._2, (log(E(i)._1) - log(e._1)) / (log(E(i)._2) - log(e._2))) })}
  }

  private def getE(i : Int, h : Double, in : (Double => Double,Double => Double,Double => Double, Double, Double),
           E : List[(Double, Double, Double)] = List()) : List[(Double, Double, Double)] = {

    def getNext(f : Double => Double, df : Double => Double, ddf : Double => Double,
               a : Double, b : Double, h : Double, coeffs : List[Tuple5]) : (Double, Double, Double) = {

      val tau = (a to b by h/10 toList) map { tau => (tau, findInterval(coeffs map (_._1), tau)) }

      List((s,f), (ds,df), (dds,ddf)) map {case (g, f) => (tau map (t =>abs(g(coeffs(t._2),t._1)-f(t._1)))) reduce (_ max _) } match {
        case a :: b :: c :: Nil => (a, b, c)
      }
    }
    if (i==10) E
    else getE(i+1, h/2, in, E :+ getNext(in._1, in._2, in._3, in._4, in._5, h, interpolationA(in._1, in._4, in._5, h)))
  }

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