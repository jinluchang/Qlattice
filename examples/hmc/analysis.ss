#!/usr/bin/env scheme-script

(import (cslib all))

; http://arxiv.org/abs/hep-lat/9806005v2
; http://arxiv.org/abs/hep-lat/0309017v1

(define (get-a-wilson beta)
  (let* ([a1 -1.6805]
         [a2 -1.7139]
         [a3 0.8155]
         [a4 -0.6667]
         [b-6 (- beta 6)]
         [b-6^2 (sqr b-6)]
         [b-6^3 (* b-6 b-6^2)]
         [lna/r0 (+ a1 (* a2 b-6) (* a3 b-6^2) (* a4 b-6^3))])
    (* 0.5 (exp lna/r0))))

(define (get-a-iwasaki beta)
  (let* ([c1 -2.1281]
         [c2 -1.0056]
         [c3 0.6041]
         [b-3 (- beta 3)]
         [b-3^2 (sqr b-3)]
         [lna/r0 (+ c1 (* c2 b-3) (* c3 b-3^2))])
    (* 0.5 (exp lna/r0))))

(define (get-a-dbw2 beta)
  (let* ([d1 -1.6007]
         [d2 -2.3179]
         [d3 -0.8020]
         [d4 -19.8509]
         [b-1 (- beta 1)]
         [b-1^2 (sqr b-1)]
         [b-1^3 (* b-1 b-1^2)]
         [lna/r0 (+ d1 (* d2 b-1) (* d3 b-1^2) (* d4 b-1^3))])
    (* 0.5 (exp lna/r0))))

(define sqrt-t0/fm 0.1667) ; in unit of fm

(define (get-t0/a a/fm)
  (sqr (/ sqrt-t0/fm a/fm)))

(print "hello world")

(print (get-a-iwasaki 2.31))

(print (get-a-dbw2 0.7796))

(print "sqrt(t0)" (sqrt (* 2.7854 (sqr 0.0999))))
(print "sqrt(t0)" (sqrt (* 5.489 (sqr 0.0710))))
(print "sqrt(t0)" (sqrt (* 11.241 (sqr 0.0498))))

(print "sqrt(t0) = 0.1667")

(print (get-t0/a (get-a-dbw2 0.7796)))
(print (get-t0/a 0.2))
(print (get-t0/a 0.1))
(print (get-t0/a 0.125))

(print "a/fm")

(print (get-a-wilson 5.95935)) ; 0.1 fm
(print (get-a-iwasaki 2.5868)) ; 0.1 fm
(print (get-a-dbw2 1.0038)) ; 0.1 fm

(print (get-a-dbw2 0.91082)) ; 0.125 fm

(print (get-a-wilson 5.63594))
(print (get-a-dbw2 0.7796))

(print (get-a-wilson 5.8139))
(print (get-a-dbw2 0.89))
(print (get-t0/a (get-a-dbw2 0.89)))
