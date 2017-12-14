val xxx = List((1,2), (5,4))
xxx.map {case(i, j) => (i+1, j+1)}


val fraction: PartialFunction[Int, Int] =
{ case d: Int if d != 0 ⇒ 42 / d }

fraction(5)