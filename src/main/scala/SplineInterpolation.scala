package numerik_praktikum

import Math._

object SplineInterpolation {

  type Tuple5 = (Double, Double, Double, Double, Double)

  //conversion method
  def flatten(t: ((Int, Double, Double), Double)) : (Int, Double, Double, Double) = (t._1._1, t._1._2, t._1._3, t._2 )

  val s : (Tuple5, Double) => Double = (c : Tuple5, t : Double) =>
    c._2 + c._3*(t-c._1) + c._4*(t-c._1)*(t-c._1)/2 + c._5*(t-c._1)*(t-c._1)*(t-c._1)/6

  val ds : (Tuple5, Double) => Double = (c : Tuple5, t : Double) => c._3 + c._4*(t-c._1) + c._5*(t-c._1)*(t-c._1)/2

  val dds : (Tuple5, Double) => Double = (c : Tuple5, t : Double) => c._4 + c._5*(t-c._1)

  val epsilon = 0.001

  var testValues = List(
    ((t : Double) => sin(t), -PI, PI),
    ((t : Double) => cos(t), -PI, PI),
    ((t : Double) => 1/3*abs(t)*sin(t)*sin(t), -PI, PI),
    ((t : Double) => sqrt(t), 1, 2),
    ((t : Double) => exp(t), -1, 1),
    ((t : Double) => 1/(1+t*t), -1, 1))

  def main(args : Array[String]) {
    //testValues = List((, , PI))

    println("Interpolation Typ (a) - natuerlicher Spline\n")
    val matrix = interpolationCaseA((t : Double) => sin(t), -PI, PI, 1)


    println(getE_h(testValues(0)._1, matrix, 1))

    //println(gaussWithPivot(List(List(4,-2,-2,10), List(8,-5,1,6), List(-2,7,-2,3))))
    //println(gaussWithPivot(List(List(2,-3,-5), List(-2,4,8))))
  }


  private def interpolationCaseA(f : Double => Double, a : Double, b : Double, h : Double) : List[Tuple5] = {

    var values = (a to b by h toList) zip (Stream from 0) map (_.swap) map { case (i, t) => (i,t,f(t)) }

    //only for testing purposes
  //  values = List((0, -2, -1), (1, -1, 0), (2, 0, 1), (3, 1, 0), (4, 2, -1))


    val p = values.map{x => (x._1, x._3)}.toMap
    val n = p.size

    val coeffMatrix = (List.fill(n+1)(0.0) updated(0, 1.0)) :: (List.fill(n+1)(0.0) updated(n-1, 1.0)) ::
      (1 to n-2 toList).map{ i => List.fill(n+1)(0.0) updated (i-1, h/6) updated (i, 2*h/3) updated (i+1, h/6) updated (n, (p(i-1)+p(i+1)-2*p(i))/h)}

    values zip gaussWithPivot(coeffMatrix) map flatten sliding(2) map { case c :: n :: Nil =>
      (c._2, c._3, (n._2-c._2)/h-h*(2*c._4+n._4)/6, c._4, (n._4-c._4)/h)} toList

  }

  private def getE_h(f : Double => Double, coeffs : List[Tuple5], h : Double) : (Double, Double, Double) = {

    val all = (0 until 10*coeffs.length toList) map { j => coeffs(0)._1 + h*j/10 } map { t =>
      (t, findInterval(coeffs map (_._1), t)) } map { case (t,i) =>
      (abs(s(coeffs(i), t)-f(t)), abs(ds(coeffs(i), t)-f(t)), abs(dds(coeffs(i), t)-f(t))) }

    val max = (y : List[Double]) => y reduce (_ max _)

    (max(all map (_._1)), max(all map (_._2)), max(all map (_._3)))
  }

  private def gaussWithPivot(A : List[List[Double]]) : List[Double] = {

    def forwardSub(A : List[List[Double]], nextCoordinate : (Int,Int)): List[List[Double]] = nextCoordinate match {

      case (i,j) if i==j && i==A.length-1 => A

      case (i,j) if i==j => {
        val pivot = (j + 1 until A.length foldLeft j) { (pivot,i) => if (A(i)(j).abs > A(pivot)(j).abs) i else pivot }
        forwardSub(A.updated(pivot, A(j)).updated(j, A(pivot)), (i+1,j))
      }

      case (i,j) => { val factor = A(i)(j)/A(j)(j)
        val newRow = (0 to A.length toList) map { _ match {
            case k if k <= j => 0
            case k => A(i)(k) - A(j)(k)*factor
          }
        }
        if (i+1<A.length) forwardSub(A.updated(i, newRow), (i+1,j)) else forwardSub(A.updated(i, newRow), (j+1,j+1))
      }
    }

    def backwardSub(A : List[List[Double]], x : List[Double], nextCoordinate : (Int,Int)) : List[Double] = nextCoordinate match {

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